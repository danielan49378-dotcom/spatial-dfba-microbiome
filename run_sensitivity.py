# run_sensitivity.py
# GSA script for testing population density vs kinetic bottlenecks

import os
import warnings
import pandas as pd
from model_loader import load_agora_model
from dfba_engine import DFBASimulator, KineticParams
from prcc_sensitivity import run_global_sensitivity, plot_prcc

# mute cobra/glpk solver warnings
import logging
logging.getLogger('cobra').setLevel(logging.CRITICAL)
warnings.filterwarnings('ignore')

def simulate_yield(params):
    b_longum_path = os.path.join("models", "B_longum.xml")
    a_hallii_path = os.path.join("models", "A_hallii.xml")
    
    bounds = {
        'EX_h2o_e_': (-1000.0, 1000.0), 'EX_h_e_': (-1000.0, 1000.0),
        'EX_glc_D_e_': (0.0, 0.0), 'EX_ac_e_': (0.0, 1000.0), 'EX_but_e_': (0.0, 1000.0)
    }
    
    m_bl = load_agora_model(b_longum_path, bounds)
    m_ah = load_agora_model(a_hallii_path, bounds)

    kinetics = {
        'b_longum': {
            'EX_glc_D_e_': KineticParams(10.0, 0.5)
        },
        'a_hallii': {
            'EX_glc_D_e_': KineticParams(8.0, 0.6),
            'EX_ac_e_': KineticParams(params['Vmax_AH_Ac'], 1.2)  
        }
    }

    # sample the starting population densities
    init_x = {
        'b_longum': params['Init_BL_Biomass'], 
        'a_hallii': params['Init_AH_Biomass']
    }
    init_s = {'glc': 50.0, 'ac': 50.0, 'but': 0.0}

    sim = DFBASimulator(
        m_bl, m_ah,
        init_x, init_s,
        kinetics, dt=0.5 
    )
    
    try:
        res = sim.simulate(t_max=24.0)
        return float(res.iloc[-1]['S_but'])
    except Exception:
        return 0.0

if __name__ == "__main__":
    # LHS sampling bounds
    bounds = {
        'Init_AH_Biomass': (0.01, 0.20),
        'Init_BL_Biomass': (0.01, 0.20),
        'Vmax_AH_Ac': (5.0, 20.0)
    }

    print("Starting parallel GSA for butyrate drivers...")

    prcc_df = run_global_sensitivity(
        simulator_fn=simulate_yield,
        bounds=bounds,
        n_samples=200,      
        max_workers=30
    )

    if not prcc_df.empty:
        prcc_df.to_csv("prcc_results.csv")
        plot_prcc(prcc_df, title="PRCC: Population vs Kinetics")
        print("Saved GSA results to prcc_results.csv")