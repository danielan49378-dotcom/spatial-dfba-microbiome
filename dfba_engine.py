# Core dFBA math engine for cross-feeding
import numpy as np
import pandas as pd
from dataclasses import dataclass

@dataclass
class KineticParams:
    v_max: float
    k_m: float

class DFBASimulator:
    def __init__(self, m_bl, m_ah, init_x, init_s, kinetics, dt=0.1):
        self.m_bl = m_bl
        self.m_ah = m_ah
        self.dt = dt
        self.kinetics = kinetics
        
        self.time = 0.0
        self.x = init_x.copy()  # biomass (gDW/L)
        self.s = init_s.copy()  # substrates (mmol/L)
        self.history = []
        self._log_state()

    def _apply_mm_kinetics(self, model, name):
        if name not in self.kinetics:
            return
            
        for rxn_id, params in self.kinetics[name].items():
            # grab base met name (e.g., EX_glc_D_e_ -> glc)
            met = rxn_id.split('_')[1] 
            s_val = max(0.0, self.s.get(met, 0.0))
            
            v_uptake = params.v_max * (s_val / (params.k_m + s_val))
            x_val = self.x[name]
            
            # prevent negative concentrations on the next step
            v_max_phys = (s_val / (x_val * self.dt)) if x_val > 0 else 0.0
            
            try:
                model.reactions.get_by_id(rxn_id).lower_bound = -min(v_uptake, v_max_phys)
            except KeyError:
                pass # skip if model lacks the exchange rxn

    def step(self):
        self._apply_mm_kinetics(self.m_bl, 'b_longum')
        self._apply_mm_kinetics(self.m_ah, 'a_hallii')
        
        sol_bl = self.m_bl.optimize()
        sol_ah = self.m_ah.optimize()
        
        mu_bl = sol_bl.objective_value if sol_bl.status == 'optimal' else 0.0
        mu_ah = sol_ah.objective_value if sol_ah.status == 'optimal' else 0.0
        
        # hard cap growth at 0.5 to avoid infinite bg carbon breaking the physics
        MAX_MU = 0.5 
        scale_bl = MAX_MU / mu_bl if mu_bl > MAX_MU else 1.0
        scale_ah = MAX_MU / mu_ah if mu_ah > MAX_MU else 1.0

        mu_bl = min(mu_bl, MAX_MU)
        mu_ah = min(mu_ah, MAX_MU)
        
        fluxes = {
            'b_longum': {
                'glc': sol_bl.fluxes.get('EX_glc_D_e_', 0.0) * scale_bl,
                'ac': sol_bl.fluxes.get('EX_ac_e_', 0.0) * scale_bl, 
                'but': sol_bl.fluxes.get('EX_but_e_', 0.0) * scale_bl
            },
            'a_hallii': {
                'glc': sol_ah.fluxes.get('EX_glc_D_e_', 0.0) * scale_ah,
                'ac': sol_ah.fluxes.get('EX_ac_e_', 0.0) * scale_ah,
                'but': sol_ah.fluxes.get('EX_but_e_', 0.0) * scale_ah
            }
        }

        # euler update for biomass
        self.x['b_longum'] += mu_bl * self.x['b_longum'] * self.dt
        self.x['a_hallii'] += mu_ah * self.x['a_hallii'] * self.dt
        
        # euler update for metabolites
        for met in self.s.keys():
            ds_dt = (fluxes['b_longum'].get(met, 0.0) * self.x['b_longum']) + \
                    (fluxes['a_hallii'].get(met, 0.0) * self.x['a_hallii'])
            
            self.s[met] = max(0.0, self.s[met] + (ds_dt * self.dt))
            
        self.time += self.dt
        self._log_state()

    def _log_state(self):
        state = {'time': self.time}
        state.update({f"X_{k}": v for k, v in self.x.items()})
        state.update({f"S_{k}": v for k, v in self.s.items()})
        self.history.append(state)
        
    def simulate(self, t_max):
        for _ in range(int(t_max / self.dt)):
            self.step()
        return pd.DataFrame(self.history)