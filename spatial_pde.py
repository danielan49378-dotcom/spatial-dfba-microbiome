# spatial_pde.py
# 2D Reaction-Diffusion solver for dFBA

import numpy as np
from scipy.ndimage import laplace
from dfba_engine import DFBASimulator

class StatelessCellSimulator(DFBASimulator):
    """Prevents the DFBASimulator from eating memory when run thousands of times per step."""
    def _log_state(self):
        pass  

class SpatialPDEGrid:
    def __init__(self, grid_shape, dx, dt, diffusion_coeffs, model_b_longum, model_a_hallii, kinetics, thresh=1e-6):
        self.nx, self.ny = grid_shape
        self.dx = dx
        self.dt = dt
        self.D = diffusion_coeffs
        self.thresh = thresh
        self.t = 0.0
        
        self.x = {
            'b_longum': np.zeros(grid_shape, dtype=np.float64),
            'a_hallii': np.zeros(grid_shape, dtype=np.float64)
        }
        
        self.s = {
            'glc': np.zeros(grid_shape, dtype=np.float64),
            'ac':  np.zeros(grid_shape, dtype=np.float64),
            'but': np.zeros(grid_shape, dtype=np.float64)
        }
        
        self._check_cfl()
        
        # init dummy engine for cell-level reactions
        dummy_x = {k: 0.0 for k in self.x.keys()}
        dummy_s = {k: 0.0 for k in self.s.keys()}
        
        self.engine = StatelessCellSimulator(
            m_bl=model_b_longum,
            m_ah=model_a_hallii,
            init_x=dummy_x,
            init_s=dummy_s,
            kinetics=kinetics,
            dt=self.dt
        )

    def _check_cfl(self):
        """Warn if time steps are too big for the grid resolution."""
        max_d = max(self.D.values())
        cfl = (max_d * self.dt) / (self.dx ** 2)
        if cfl > 0.25:
            print(f"WARNING: CFL Condition violated ({cfl:.3f} > 0.25). Math might explode.")

    def initialize_state(self, x_init, s_init):
        for k, arr in x_init.items():
            self.x[k] = arr.copy()
        for k, arr in s_init.items():
            self.s[k] = arr.copy()

    def _diffuse(self):
        lap_factor = self.dt / (self.dx ** 2)
        
        for k, grid in self.x.items():
            if self.D.get(k, 0.0) > 0:
                delta = self.D[k] * laplace(grid, mode='reflect') * lap_factor
                self.x[k] = np.clip(grid + delta, a_min=0.0, a_max=None)

        for k, grid in self.s.items():
            if self.D.get(k, 0.0) > 0:
                delta = self.D[k] * laplace(grid, mode='reflect') * lap_factor
                self.s[k] = np.clip(grid + delta, a_min=0.0, a_max=None)

    def _react(self):
        # only solve LP if there are actually bacteria there
        active = (self.x['b_longum'] > self.thresh) | (self.x['a_hallii'] > self.thresh)
        coords = np.argwhere(active)
        
        for i, j in coords:
            for k in self.x.keys():
                self.engine.x[k] = self.x[k][i, j]
            for k in self.s.keys():
                self.engine.s[k] = self.s[k][i, j]
                
            self.engine.step()
            
            for k in self.x.keys():
                self.x[k][i, j] = self.engine.x[k]
            for k in self.s.keys():
                self.s[k][i, j] = self.engine.s[k]

    def step(self):
        self._diffuse()
        self._react()
        self.t += self.dt

    def simulate(self, t_max):
        steps = int(t_max / self.dt)
        history = []
        
        for step in range(steps):
            self.step()
            # save snapshot every 10 steps to save RAM
            if step % 10 == 0 or step == steps - 1:
                snap = {
                    'time': self.t,
                    'x': {k: v.copy() for k, v in self.x.items()},
                    's': {k: v.copy() for k, v in self.s.items()}
                }
                history.append(snap)
                
        return history