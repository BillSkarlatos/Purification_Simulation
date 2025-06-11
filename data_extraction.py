#!/usr/bin/env python3
"""
export_data_updated.py: export the metrics used in adaptive.py graphs to CSV,
using projector-based fidelity to match the paper, and also record pre-purification fidelity.
"""

import csv
import time
import numpy as np
import qutip as qt
from noise import two_qubit_noise, amp_damp_kraus, phase_damp_kraus
from purification import dejmps_purify, apply_filter

# Load lookup table from adaptive.py preprocessing
lookup = np.load('lookup_table.npy')

# Grid settings (must match adaptive.py)
gammas = np.linspace(0.01, 0.20, 20)
ps     = np.linspace(0.01, 0.20, 20)
# Probe and timing settings
N_probe   = 50
M_probe   = 50
cycle_time = 0.01  # seconds per cycle

# Prepare Bell state and projection
bell_proj = qt.bell_state('00').proj()

# Single-qubit noise (as in adaptive.py)
def single_qubit_noise(rho, gamma, p):
    A = amp_damp_kraus(gamma)
    B = phase_damp_kraus(p)
    K_single = [b * a for a in A for b in B]
    out = qt.Qobj(np.zeros((2, 2)), dims=[[2], [2]])
    for K in K_single:
        out += K * rho * K.dag()
    return out

# Channel estimation (unchanged)
def estimate_channel(gamma_true, p_true, N, M):
    n1 = 0
    proj1 = qt.basis(2,1) * qt.basis(2,1).dag()
    for _ in range(N):
        rho_out = single_qubit_noise(proj1, gamma_true, 0)
        if np.random.rand() < float((proj1 * rho_out).tr().real):
            n1 += 1
    gamma_hat = 1 - n1/N

    m_minus = 0
    plus = (qt.basis(2,0)+qt.basis(2,1)).unit()
    minus = (qt.basis(2,0)-qt.basis(2,1)).unit()
    proj_minus = minus * minus.dag()
    for _ in range(M):
        rho_out = single_qubit_noise(plus*plus.dag(), 0, p_true)
        if np.random.rand() < float((proj_minus * rho_out).tr().real):
            m_minus += 1
    p_hat = 2*m_minus/M

    return gamma_hat, p_hat, N+M

# Static purification cascade (depth d, filter alpha)
def static_purification(gamma, p, d, alpha):
    # Compute raw fidelity before any purification
    rho0 = two_qubit_noise(bell_proj, gamma, p)
    F_raw = float((bell_proj * rho0).tr().real)

    # Initialize survivors to raw pairs
    survivors = [rho0 for _ in range(2**d)]
    for _ in range(d):
        survivors = [rho_f for rho in survivors
                     for rho_f, pf in [apply_filter(rho, alpha)]
                     if rho_f is not None and np.random.rand() < pf]
        np.random.shuffle(survivors)
        new = []
        for idx in range(0, len(survivors)-1, 2):
            rho_p, pp = dejmps_purify(survivors[idx], survivors[idx+1])
            if rho_p is not None and np.random.rand() < pp:
                new.append(rho_p)
        survivors = new
        if not survivors:
            break

    init, final = 2**d, len(survivors)
    # Yield: fraction of surviving pairs
    Y = final/init if init > 0 else 0
    # Purified fidelity: average overlap of survivors if any; otherwise zero
    F_static = np.mean([float((bell_proj * rho).tr().real) for rho in survivors]) if final > 0 else 0.0
    return F_raw, F_static, Y

# Helper: nearest index on grid
def find_index(arr, val): return np.abs(arr - val).argmin()

# Main export routine
def main():
    start = time.time()
    with open('adaptive_data_updated.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            'gamma', 'p',
            'd_static', 'alpha_static', 'F_raw', 'F_static', 'Y_static', 'R_static',
            'd_adaptive', 'alpha_adaptive', 'F_adaptive', 'Y_adaptive', 'R_adaptive',
            'delta_F', 'delta_Y', 'delta_T'
        ])
        for i, g in enumerate(gammas):
            for j, p in enumerate(ps):
                # Static
                d_st, a_st = int(lookup[i,j,0]), float(lookup[i,j,1])
                F_raw, F_st, Y_st = static_purification(g, p, d_st, a_st)
                R_st = (2**d_st)/cycle_time * Y_st

                # Adaptive
                gh, ph, probes = estimate_channel(g, p, N_probe, M_probe)
                ii = find_index(gammas, gh)
                jj = find_index(ps, ph)
                d_ad, a_ad = int(lookup[ii,jj,0]), float(lookup[ii,jj,1])
                _, F_ad, Y_ad = static_purification(g, p, d_ad, a_ad)
                R_ad = (2**d_ad)/cycle_time * Y_ad - probes/cycle_time

                # Differences
                dF = F_ad - F_st
                dY = Y_ad - Y_st
                dT = R_ad - R_st

                writer.writerow([
                    g, p,
                    d_st, a_st, F_raw, F_st, Y_st, R_st,
                    d_ad, a_ad, F_ad, Y_ad, R_ad,
                    dF, dY, dT
                ])
    print(f"Export complete in {time.time()-start:.2f}s to 'adaptive_data_updated.csv'.")

if __name__ == '__main__':
    np.random.seed(42)
    main()