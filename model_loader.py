# model_loader.py
# Safely loads AGORA2 SBML files and bypasses generic mass imbalances
import logging
import cobra
from cobra.util.array import create_stoichiometric_matrix

logger = logging.getLogger('model_loader')

def validate_matrix(model, tol=1e-6):
    """Checks S-matrix and logs generic R-group reactions without crashing."""
    s_mat = create_stoichiometric_matrix(model)
    logger.info(f"Loaded {model.id} S-matrix: {s_mat.shape[0]} mets x {s_mat.shape[1]} rxns")
    
    bad_rxns = []
    for rxn in model.reactions:
        if rxn.boundary: 
            continue
            
        # catch pseudo-reactions and R-groups common in AGORA2
        bal = rxn.check_mass_balance()
        if any(abs(v) > tol for v in bal.values()):
            bad_rxns.append(rxn.id)
            
    if bad_rxns:
        logger.warning(f"{model.id}: Found {len(bad_rxns)} imbalanced rxns. Proceeding anyway.")

def set_bounds(model, bounds):
    """Applies environment constraints to exchange boundaries."""
    for rxn_id, (lb, ub) in bounds.items():
        try:
            rxn = model.reactions.get_by_id(rxn_id)
            rxn.bounds = (lb, ub)
        except KeyError:
            # quietly skip if model doesn't have this specific exchange
            pass

def load_agora_model(file_path, bounds):
    model = cobra.io.read_sbml_model(file_path)
    validate_matrix(model)
    set_bounds(model, bounds)
    model.solver = 'glpk' 
    return model