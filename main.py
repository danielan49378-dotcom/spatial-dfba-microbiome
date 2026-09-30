# execution script for 2D spatial dFBA microbiome simulation
import os
import pickle
import warnings
import numpy as np
from model_loader import load_agora_model
from dfba_engine import KineticParams # Updated to match our new shorter dataclass name!
from spatial_pde import SpatialPDEGrid
from visualizer import plot_spatial_snapshot

# mute cobra solver warnings
import logging
logging.getLogger('cobra').setLevel(logging.CRITICAL)
warnings.filterwarnings('ignore')

def main():
    print("Loading models and setting up bounds...")
    b_longum_path = os.path.join("models", "B_longum.xml")
    a_hallii_path = os.path.join("models", "A_hallii.xml")

    # open acetate and butyrate for secretion
    bounds = {
        'EX_h2o_e_': (-1000.0, 1000.0), 
        'EX_h_e_': (-1000.0, 1000.0),
        'EX_glc_D_e_': (0.0, 0.0),  
        'EX_ac_e_': (0.0, 1000.0),  
        'EX_but_e_': (0.0, 1000.0)  
    }

    m_bl = load_agora_model(b_longum_path, bounds)
    m_ah = load_agora_model(a_hallii_path, bounds)

    # Vmax and Km for Michaelis-Menten kinetics
    kinetics = {
        'b_longum': {
            'EX_glc_D_e_': KineticParams(10.0, 0.5)
        },
        'a_hallii': {
            'EX_glc_D_e_': KineticParams(8.0, 0.6),
            'EX_ac_e_': KineticParams(12.0, 1.2)  
        }
    }

    # approximate diffusion rates (cm^2/h)
    diff_coeffs = {
        'b_longum': 1e-5,  
        'a_hallii': 1e-5,
        'glc': 0.03,       
        'ac': 0.04,
        'but': 0.035
    }

    grid = (50, 50)
    
    print("Initializing 2D PDE grid...")
    pde = SpatialPDEGrid(
        grid_shape=grid,
        dx=0.1,             
        dt=0.05, # keeps CFL stable
        diffusion_coeffs=diff_coeffs,
        model_b_longum=m_bl,
        model_a_hallii=m_ah,
        kinetics=kinetics
    )

    # seed bacteria in the exact center
    cx, cy = grid[0]//2, grid[1]//2
    
    x_init = {
        'b_longum': np.zeros(grid),
        'a_hallii': np.zeros(grid)
    }
    x_init['b_longum'][cx, cy] = 0.1
    x_init['a_hallii'][cx, cy] = 0.1

    # blanket the grid in glucose
    s_init = {
        'glc': np.ones(grid) * 50.0,  
        'ac': np.zeros(grid),
        'but': np.zeros(grid)
    }

    pde.initialize_state(x_init, s_init)

    print("Running 24h spatial simulation... (this may take a few minutes)")
    history = pde.simulate(t_max=24.0)
    
    print("Simulation complete. Saving artifacts...")
    
    # generate heatmap png
    plot_spatial_snapshot(history, step_index=-1)

    # dump raw data for later analysis
    with open("simulation_data.pkl", "wb") as f:
        pickle.dump(history, f)
        
    print("Done. Saved spatial_snapshot_t24.0.png and simulation_data.pkl.")

if __name__ == "__main__":
    main()