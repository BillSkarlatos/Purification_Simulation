# data_extraction.py

import csv
import numpy as np
import qutip as qt
from noise import two_qubit_noise
from throughput_improvement import run_static, run_adaptive

"""
This script computes static (fixed-depth DEJMPS) and adaptive (purified) fidelity
and yield over a grid of (gamma, p) noise parameters. It outputs a CSV file
with columns:
    gamma, p,
    F_static, Y_static,
    F_adapt,  Y_adapt,
    ΔF (F_adapt − F_static),
    ΔY (Y_adapt − Y_static)

Adjust GRID_SIZE and TRIALS as needed. By default, we use a 6×6 grid (36 points)
over γ,p ∈ [0.01, 0.20] and 200 Monte Carlo trials per point.
"""

# ===== USER PARAMETERS =====
GRID_SIZE = 6       # creates GRID_SIZE × GRID_SIZE combinations (must be ≥ 6 to get ≥ 30 points)
TRIALS = 200        # Monte Carlo trials per (gamma, p)
STATIC_DEPTH = 2    # fixed DEJMPS rounds for static protocol
LOOKUP_DEPTH = 2    # fixed depth for adaptive protocol (as in run_adaptive)

OUTPUT_CSV = "data_extraction_results.csv"
# =============================

def simulate_point(gamma: float, p: float, trials: int, static_depth: int, lookup_depth: int):
    """
    For given (gamma, p), perform 'trials' Monte Carlo runs to estimate:
      - static fidelity & yield (depth = static_depth)
      - adaptive fidelity & yield (lookup_depth for run_adaptive)
    Returns:
      avg_F_static, avg_Y_static, avg_F_adapt, avg_Y_adapt
    """
    fidelities_static = []
    yields_static = []
    fidelities_adapt = []
    yields_adapt = []

    # Pre-construct noisy Bell pairs each trial inside loops
    bell = qt.bell_state('00')

    for _ in range(trials):
        # Generate two fresh noisy Bell pairs for static run
        rho1_static = two_qubit_noise(bell.proj(), gamma, p)
        rho2_static = two_qubit_noise(bell.proj(), gamma, p)
        f_s, y_s = run_static(rho1_static, rho2_static, gamma, p, static_depth)
        fidelities_static.append(f_s)
        yields_static.append(y_s)

        # Generate two fresh noisy Bell pairs for adaptive run
        rho1_ad = two_qubit_noise(bell.proj(), gamma, p)
        rho2_ad = two_qubit_noise(bell.proj(), gamma, p)
        f_a, y_a = run_adaptive(rho1_ad, rho2_ad, gamma, p, lookup_depth)
        fidelities_adapt.append(f_a)
        yields_adapt.append(y_a)

    avg_F_static = np.mean(fidelities_static)
    avg_Y_static = np.mean(yields_static)
    avg_F_adapt  = np.mean(fidelities_adapt)
    avg_Y_adapt  = np.mean(yields_adapt)

    return avg_F_static, avg_Y_static, avg_F_adapt, avg_Y_adapt


def main():
    # Define grid of gamma and p values
    gammas = np.linspace(0.01, 0.20, GRID_SIZE)
    ps     = np.linspace(0.01, 0.20, GRID_SIZE)

    # Open CSV for writing
    with open(OUTPUT_CSV, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        # Header
        writer.writerow([
            "gamma", "p",
            "F_static", "Y_static",
            "F_adapt", "Y_adapt",
            "ΔF",       "ΔY"
        ])

        # Loop over grid
        for gamma in gammas:
            for p_val in ps:
                F_s, Y_s, F_a, Y_a = simulate_point(
                    gamma, p_val, TRIALS, STATIC_DEPTH, LOOKUP_DEPTH
                )
                delta_F = F_a - F_s
                delta_Y = Y_a - Y_s

                writer.writerow([
                    f"{gamma:.4f}",
                    f"{p_val:.4f}",
                    f"{F_s:.6f}",
                    f"{Y_s:.6f}",
                    f"{F_a:.6f}",
                    f"{Y_a:.6f}",
                    f"{delta_F:.6f}",
                    f"{delta_Y:.6f}"
                ])

    print(f"Data extraction complete. Results saved to '{OUTPUT_CSV}'.")


if __name__ == "__main__":
    main()
