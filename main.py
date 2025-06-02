import numpy as np
import qutip as qt
from noise import two_qubit_noise
from simulate import sweep_parameters
from visualize import (
    plot_fidelity_surface,
    plot_yield_surface,
    plot_fidelity_contour,
    plot_yield_contour,
    plot_fidelity_vs_gamma,
    plot_yield_vs_gamma,
    plot_fidelity_vs_p,
    plot_yield_vs_p,
    # Newly added plotting functions for differences:
    plot_diff_surface,
    plot_diff_contour,
    plot_diff_vs_gamma,
    plot_diff_vs_p
)

def compute_default_metrics(gammas, ps):
    """
    Compute default (unmitigated) fidelity and yield for each (gamma, p).
    default_fidelity = fidelity of noisy Bell pair
    default_yield = 1.0 (no purification)
    """
    results_default = {}
    bell = qt.bell_state('00')  # |Φ⁺⟩
    for gamma in gammas:
        for p in ps:
            noisy_rho = two_qubit_noise(bell.proj(), gamma, p)
            default_fidelity = float(qt.fidelity(noisy_rho, bell))
            default_yield = 1.0
            results_default[(gamma, p)] = {
                'fidelity': default_fidelity,
                'yield': default_yield
            }
    return results_default

def compute_difference(results_default, results_purify):
    """
    Given default and purification results, compute difference:
      diff_fidelity = F_purify - F_default
      diff_yield    = yield_purify - yield_default (where yield_default=1)
    """
    diff_results = {}
    for key in results_purify:
        default = results_default[key]
        purify = results_purify[key]
        diff_results[key] = {
            'diff_fidelity': purify['avg_fidelity'] - default['fidelity'],
            'diff_yield': purify['avg_yield'] - default['yield']
        }
    return diff_results

def main():
    # Define the noise‐parameter grid and simulation settings
    gammas = np.linspace(0.01, 0.2, 10)   # 10 values from 0.01 to 0.20
    ps = np.linspace(0.01, 0.2, 10)       # 10 values from 0.01 to 0.20
    depth = 2                             # Purification depth
    trials = 500                          # Monte Carlo trials per (gamma, p)

    # 1. Compute default (unmitigated) performance
    results_default = compute_default_metrics(gammas, ps)

    # 2. Run Monte Carlo sweep for adaptive purification
    results_purify = sweep_parameters(gammas, ps, depth, trials)

    # 3. Compute Δ (difference) between purified and default
    diff_results = compute_difference(results_default, results_purify)

    # --- Plotting Purified Results ---

    # 3D surface plots (Adaptive Purification)
    plot_fidelity_surface(results_purify, gammas, ps)
    plot_yield_surface(results_purify, gammas, ps)

    # 2D contour plots (Adaptive Purification)
    plot_fidelity_contour(results_purify, gammas, ps)
    plot_yield_contour(results_purify, gammas, ps)

    # 2D line plots vs. gamma for fixed p (Adaptive Purification)
    fixed_p_idx = 0  # index into ps for p
    plot_fidelity_vs_gamma(results_purify, gammas, ps, fixed_p_index=fixed_p_idx)
    plot_yield_vs_gamma(results_purify, gammas, ps, fixed_p_index=fixed_p_idx)

    # 2D line plots vs. p for fixed gamma (Adaptive Purification)
    fixed_gamma_idx = 0  # index into gammas for gamma
    plot_fidelity_vs_p(results_purify, gammas, ps, fixed_gamma_index=fixed_gamma_idx)
    plot_yield_vs_p(results_purify, gammas, ps, fixed_gamma_index=fixed_gamma_idx)

    # --- Plotting Differences (Adaptive minus Default) ---

    # 3D surface plots: Δ Fidelity and Δ Yield
    plot_diff_surface(diff_results, gammas, ps, metric='diff_fidelity', title='Surface: Δ Fidelity')
    plot_diff_surface(diff_results, gammas, ps, metric='diff_yield', title='Surface: Δ Yield')

    # 2D contour plots: Δ Fidelity and Δ Yield
    plot_diff_contour(diff_results, gammas, ps, metric='diff_fidelity', title='Contour: Δ Fidelity')
    plot_diff_contour(diff_results, gammas, ps, metric='diff_yield', title='Contour: Δ Yield')

    # 2D line plots vs. gamma for fixed p: Δ Fidelity and Δ Yield
    plot_diff_vs_gamma(
        diff_results, gammas, ps,
        metric='diff_fidelity',
        fixed_p_index=fixed_p_idx,
        title=f'Δ Fidelity vs Gamma (p={ps[fixed_p_idx]:.2f})'
    )
    plot_diff_vs_gamma(
        diff_results, gammas, ps,
        metric='diff_yield',
        fixed_p_index=fixed_p_idx,
        title=f'Δ Yield vs Gamma (p={ps[fixed_p_idx]:.2f})'
    )

    # 2D line plots vs. p for fixed gamma: Δ Fidelity and Δ Yield
    plot_diff_vs_p(
        diff_results, gammas, ps,
        metric='diff_fidelity',
        fixed_gamma_index=fixed_gamma_idx,
        title=f'Δ Fidelity vs p (γ={gammas[fixed_gamma_idx]:.2f})'
    )
    plot_diff_vs_p(
        diff_results, gammas, ps,
        metric='diff_yield',
        fixed_gamma_index=fixed_gamma_idx,
        title=f'Δ Yield vs p (γ={gammas[fixed_gamma_idx]:.2f})'
    )

if __name__ == '__main__':
    main()
