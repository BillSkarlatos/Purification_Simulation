import numpy as np
import csv
from main import compute_default_metrics, compute_difference
from simulate import sweep_parameters

def main():
    # Define the noise‐parameter grid and simulation settings (matching main.py)
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

    # 4. Write all data to a CSV file
    # Columns: γ, p, F, Y, ΔF, ΔY
    with open('data.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['gamma', 'p', 'F', 'Y', 'delta_F', 'delta_Y'])
        for gamma in gammas:
            for p in ps:
                F = results_purify[(gamma, p)]['avg_fidelity']
                Y = results_purify[(gamma, p)]['avg_yield']
                delta_F = diff_results[(gamma, p)]['diff_fidelity']
                delta_Y = diff_results[(gamma, p)]['diff_yield']
                writer.writerow([gamma, p, F, Y, delta_F, delta_Y])

    print("Data extraction complete. Saved to data.csv.")

if __name__ == '__main__':
    main()
