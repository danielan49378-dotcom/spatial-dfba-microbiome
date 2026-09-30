# test_growth.py
# Unit test to verify A. hallii can physically grow and produce butyrate in isolation

import os
import cobra
import logging

# mute solver warnings
logging.getLogger("cobra").setLevel(logging.ERROR)

def test_growth():
    print("Loading A. hallii...")
    model = cobra.io.read_sbml_model(os.path.join("models", "A_hallii.xml"))
    
    # force feed glucose
    try:
        model.reactions.get_by_id("EX_glc_D_e_").lower_bound = -10.0
    except KeyError:
        print("Error: No glucose exchange found.")

    sol = model.optimize()
    
    print("\n--- SANITY CHECK ---")
    print(f"Growth (mu): {sol.objective_value:.4f}")
    
    if sol.objective_value > 0:
        try:
            print(f"Butyrate Flux: {sol.fluxes['EX_but_e_']:.4f}")
        except KeyError:
            print("Error: Butyrate exchange missing.")
    else:
        print("Bacteria starved. Check if background amino acids/vitamins are blocked.")

if __name__ == "__main__":
    test_growth()