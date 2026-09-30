# prcc_sensitivity.py
# Uses Latin Hypercube Sampling and PRCC to find correlation between parameters and metabolite yield

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import rankdata, t, qmc
from concurrent.futures import ProcessPoolExecutor, as_completed

def generate_lhs(bounds, n_samples, seed=42):
    k = len(bounds)
    sampler = qmc.LatinHypercube(d=k, seed=seed)
    sample_unit = sampler.random(n=n_samples)
    
    l_bounds = np.array([b[0] for b in bounds.values()])
    u_bounds = np.array([b[1] for b in bounds.values()])
    
    scaled = qmc.scale(sample_unit, l_bounds, u_bounds)
    return pd.DataFrame(scaled, columns=list(bounds.keys()))

def calculate_prcc(X, y, param_names):
    N, k = X.shape
    
    if N <= k + 2:
        raise ValueError("Not enough samples for degrees of freedom.")

    # If biology flatlines (e.g. bacteria starve/max out), variance is zero and math explodes. 
    # This catches it before numpy crashes.
    if np.std(y) < 1e-8:
        print("Biological failure: Output variance is effectively zero. Check growth limits.")
        return pd.DataFrame()

    D = np.column_stack((X, y))
    ranked_D = np.apply_along_axis(rankdata, 0, D)
    C = np.corrcoef(ranked_D, rowvar=False)
    
    try:
        W = np.linalg.inv(C)
    except np.linalg.LinAlgError:
        W = np.linalg.pinv(C)
        
    prcc_vals = np.zeros(k)
    p_vals = np.zeros(k)
    dof = N - k - 1
    
    for i in range(k):
        r = -W[i, -1] / np.sqrt(W[i, i] * W[-1, -1])
        prcc_vals[i] = r
        
        if r**2 > 0.999999:
            t_stat = np.inf
            p_val = 0.0
        else:
            t_stat = r * np.sqrt(dof / (1.0 - r**2))
            p_val = 2 * t.sf(np.abs(t_stat), df=dof)
            
        p_vals[i] = p_val
        
    return pd.DataFrame({'Parameter': param_names, 'PRCC': prcc_vals, 'p_value': p_vals}).set_index('Parameter')

def plot_prcc(prcc_df, title="PRCC Sensitivity Analysis", p_thresh=0.05):
    if prcc_df.empty:
        print("No PRCC data to plot.")
        return
        
    df_sorted = prcc_df.reindex(prcc_df['PRCC'].abs().sort_values(ascending=True).index)
    colors = ['#1f77b4' if p < p_thresh else '#d62728' for p in df_sorted['p_value']]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(df_sorted.index, df_sorted['PRCC'], color=colors, edgecolor='black')
    ax.axvline(0, color='black', linewidth=1.2)
    ax.set_xlabel('Partial Rank Correlation Coefficient (PRCC)')
    ax.set_title(title)
    ax.set_xlim(-1.1, 1.1)
    plt.tight_layout()
    plt.show()

def run_global_sensitivity(simulator_fn, bounds, n_samples, max_workers=4):
    lhs_df = generate_lhs(bounds, n_samples)
    results = np.zeros(n_samples)
    
    print(f"Spawning {max_workers} workers for {n_samples} LHS samples...")
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_idx = {executor.submit(simulator_fn, row): idx for idx, row in lhs_df.iterrows()}
        for future in as_completed(future_to_idx):
            idx = future_to_idx[future]
            try:
                results[idx] = future.result()
            except Exception:
                results[idx] = np.nan
                
    valid_mask = ~np.isnan(results)
    valid_X = lhs_df.values[valid_mask]
    valid_y = results[valid_mask]
    
    return calculate_prcc(valid_X, valid_y, list(lhs_df.columns))