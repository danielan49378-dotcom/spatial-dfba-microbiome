# test_variance.py
# Diagnostic script to verify A. hallii butyrate production when artificially fed acetate

import os
import logging
from model_loader import load_agora_model
from dfba_engine import DFBASimulator, KineticParams

# mute solver logs so we can actually see the output
logging.getLogger('cobra').setLevel(logging.ERROR)
logging.getLogger('optlang').setLevel(logging.ERROR)
logging.getLogger('model_loader').setLevel(logging.ERROR)

def run_diagnostics():
    print("Loading models...")
    bounds = {
        'EX_h2o_e_': (-1000.0, 1000.0), 'EX_h_e_': (-1000.0, 1000.0),
        'EX_glc_D_e_': (0.0, 0.0), 'EX_ac_e_': (0.0, 1000.0), 'EX_but_e_': (0.0, 1000.0)
    }
    
    m_bl = load_agora_model(os.path.join("models", "B_longum.xml"), bounds)
    m_ah = load_agora_model(os.path.join("models", "A_hallii.xml"), bounds)

    print("\nTesting increasing glucose with 50.0 acetate seed...")
    
    for init_glc in [1.0, 5.0, 10.0, 15.0, 20.0]:
        kinetics = {
            'b_longum': {'EX_glc_D_e_': KineticParams(10.0, 0.5)},
            'a_hallii': {
                'EX_glc_D_e_': KineticParams(8.0, 0.6),
                'EX_ac_e_': KineticParams(12.0, 1.2)  
            }
        }
        
        # seed 50.0 acetate so A. hallii doesn't starve waiting for B. longum
        init_x = {'b_longum': 0.1, 'a_hallii': 0.1}
        init_s = {'glc': init_glc, 'ac': 50.0, 'but': 0.0}
        
        sim = DFBASimulator(
            m_bl, m_ah, 
            init_x, init_s, 
            kinetics, dt=0.25
        )
        
        df = sim.simulate(t_max=24.0)
        final_but = df.iloc[-1]['S_but']
        print(f"Glc: {init_glc:4.1f} | Final Butyrate: {final_but:.4f}")

if __name__ == "__main__":
    run_diagnostics()