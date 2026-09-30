# quick utility to figure out AGORA2's weird metabolite naming conventions
import os
import cobra
import logging

# mute cobra warnings
logging.getLogger('cobra').setLevel(logging.CRITICAL)

def find_mets():
    path = os.path.join("models", "B_longum.xml")
    print(f"scanning {path} for exchange reactions...")
    
    model = cobra.io.read_sbml_model(path)
    
    print("\n--- AGORA2 Exchange Mappings ---")
    for rxn in model.exchanges:
        id_l = rxn.id.lower()
        name_l = rxn.name.lower()
        
        if 'glucose' in name_l or 'glc' in id_l:
            print(f"Glc: {rxn.id} | {rxn.name}")
        elif 'acetate' in name_l or ('ac' in id_l and 'mac' not in id_l):
            print(f"Ac:  {rxn.id} | {rxn.name}")
        elif 'butyrate' in name_l or 'but' in id_l:
            print(f"But: {rxn.id} | {rxn.name}")
        elif 'water' in name_l or 'h2o' in id_l:
            print(f"H2O: {rxn.id} | {rxn.name}")

if __name__ == "__main__":
    find_mets()