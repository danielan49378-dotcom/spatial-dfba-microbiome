# visualizer.py
# generates 2D heatmaps of the spatial dFBA grid

import matplotlib.pyplot as plt

def plot_spatial_snapshot(history, step_index=-1):
    # grab the requested snapshot (default to the final state)
    snap = history[step_index]
    t = snap['time']
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(f'Spatial Cross-Feeding (t = {t:.1f}h)', fontsize=18, fontweight='bold')
    
    # --- biomass (row 1) ---
    im0 = axes[0, 0].imshow(snap['x']['b_longum'], cmap='Greens', origin='lower')
    axes[0, 0].set_title('B. longum (gDW/L)')
    fig.colorbar(im0, ax=axes[0, 0], fraction=0.046, pad=0.04)
    
    im1 = axes[0, 1].imshow(snap['x']['a_hallii'], cmap='Purples', origin='lower')
    axes[0, 1].set_title('A. hallii (gDW/L)')
    fig.colorbar(im1, ax=axes[0, 1], fraction=0.046, pad=0.04)
    
    axes[0, 2].axis('off') # leave top right empty
    
    # --- metabolites (row 2) ---
    im2 = axes[1, 0].imshow(snap['s']['glc'], cmap='Blues', origin='lower')
    axes[1, 0].set_title('Glucose (mmol/L)')
    fig.colorbar(im2, ax=axes[1, 0], fraction=0.046, pad=0.04)
    
    im3 = axes[1, 1].imshow(snap['s']['ac'], cmap='Oranges', origin='lower')
    axes[1, 1].set_title('Acetate (mmol/L)')
    fig.colorbar(im3, ax=axes[1, 1], fraction=0.046, pad=0.04)
    
    im4 = axes[1, 2].imshow(snap['s']['but'], cmap='Reds', origin='lower')
    axes[1, 2].set_title('Butyrate (mmol/L)')
    fig.colorbar(im4, ax=axes[1, 2], fraction=0.046, pad=0.04)
    
    # strip axis ticks for clean heatmaps
    for ax in axes.flat:
        if ax.has_data():
            ax.set_xticks([])
            ax.set_yticks([])
            
    plt.tight_layout()
    
    # save to disk instead of popping open a window to save rendering resources
    filename = f"spatial_snapshot_t{t:.1f}.png"
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close(fig) # free RAM