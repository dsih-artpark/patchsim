import math
import os

import matplotlib.pyplot as plt


def plot_patch_subplots(
    t_range,
    out_ode,
    patches,
    output_dir,
    model_name,
    patch_parameters=None,
    compartments=None,
    solver="ode",
):
    """
    Plots all patches as subplots in a single figure and saves the figure.
    """
    n = len(patches)
    if n == 0:
        raise ValueError("`patches` must contain at least one patch")
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4 * nrows))
    axes = axes.flatten() if n > 1 else [axes]
    for i, patch in enumerate(patches):
        ax = axes[i]
        # Plot the model's actual compartments; fall back to those present for this patch.
        comps = compartments or [k[: -len(f"_{i}")] for k in out_ode if k.endswith(f"_{i}")]
        for c in comps:
            ax.plot(t_range, out_ode[f"{c}_{i}"], label=c)
        solver_label = "ODE" if solver == "ode" else "Discrete"
        title = f"Patch {patch} ({solver_label})"
        if patch_parameters and patch in patch_parameters:
            params = patch_parameters[patch]
            param_str = ", ".join(f"{k}={v}" for k, v in params.items())
            title += f"\n({param_str})"
        ax.set_title(title)
        ax.set_xlabel("Time")
        ax.set_ylabel("Count")
        ax.legend()
    # Hide unused subplots (use n instead of loop index i)
    for j in range(n, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout()
    os.makedirs(output_dir, exist_ok=True)
    plt.savefig(os.path.join(output_dir, f"patch_timeseries_{model_name}_{solver}.png"))
    plt.close()
