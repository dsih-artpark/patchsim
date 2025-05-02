import matplotlib.pyplot as plt
import os
import math

def plot_patch_subplots(t_range, out_ode, patches, output_dir, model_name):
    """
    Plots all patches as subplots in a single figure and saves the figure.
    """
    n = len(patches)
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5*ncols, 4*nrows))
    axes = axes.flatten() if n > 1 else [axes]
    for i, patch in enumerate(patches):
        ax = axes[i]
        for c in ["S", "I", "R"]:
            ax.plot(t_range, out_ode[f"{c}_{i}"], label=c)
        ax.set_title(f"Patch {patch} (ODE)")
        ax.set_xlabel("Time")
        ax.set_ylabel("Count")
        ax.legend()
    # Hide unused subplots
    for j in range(i+1, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f"all_patches_{model_name}_ode.png"))
    plt.close()


