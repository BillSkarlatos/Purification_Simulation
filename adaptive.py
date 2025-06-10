#!/usr/bin/env python3
"""
adaptive.py: adaptive DEJMPS purification using precomputed static-lookup table
with single-qubit channel estimation, correct throughput accounting, and contour plots.
"""

import time
import numpy as np
import qutip as qt
from qutip import fidelity
from noise import two_qubit_noise, amp_damp_kraus, phase_damp_kraus
from purification import dejmps_purify, apply_filter
from visualize import plot_diff_surface, plot_diff_contour

# Simulation configuration
gammas    = np.linspace(0.01, 0.20, 20)
ps        = np.linspace(0.01, 0.20, 20)
Ftarget   = 0.90
# Number of probe qubits
N_probe   = 50
M_probe   = 50
# Cycle time (seconds)
cycle_time = 0.01

# Pre-constructed ideal Bell state
bell_proj  = qt.bell_state('00').proj()
bell_state = qt.bell_state('00')

# Single-qubit noise: compose amplitude- and phase-damping
def single_qubit_noise(rho, gamma, p):
    A = amp_damp_kraus(gamma)
    B = phase_damp_kraus(p)
    K_single = [b * a for a in A for b in B]
    out = qt.Qobj(np.zeros((2, 2)), dims=[[2], [2]])
    for K in K_single:
        out += K * rho * K.dag()
    return out

# Channel estimation via single-qubit probes
def estimate_channel(gamma_true, p_true, N, M):
    n1 = 0
    proj1 = qt.basis(2,1) * qt.basis(2,1).dag()
    for _ in range(N):
        rho_out = single_qubit_noise(proj1, gamma_true, 0)
        if np.random.rand() < float((proj1 * rho_out).tr()):
            n1 += 1
    gamma_hat = 1 - n1/N

    m_minus = 0
    plus = (qt.basis(2,0)+qt.basis(2,1)).unit()
    minus = (qt.basis(2,0)-qt.basis(2,1)).unit()
    proj_minus = minus*minus.dag()
    for _ in range(M):
        rho_out = single_qubit_noise(plus*plus.dag(), 0, p_true)
        if np.random.rand() < float((proj_minus * rho_out).tr()):
            m_minus += 1
    p_hat = 2*m_minus/M

    return gamma_hat, p_hat, N+M

# Static purification cascade given depth and alpha
def static_purification(gamma, p, d, alpha):
    survivors = [two_qubit_noise(bell_proj, gamma, p) for _ in range(2**d)]
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
    Y = final/init if init>0 else 0
    F = np.mean([fidelity(rho, bell_state) for rho in survivors]) if final>0 else 0.0
    return F, Y

# Load precomputed static lookup table
lookup = np.load('lookup_table.npy')

# Find nearest grid index
def find_index(arr, val):
    return np.abs(arr - val).argmin()

# Main simulation
def simulate():
    start = time.time()
    diff_results = {}

    for g in gammas:
        for p in ps:
            i = find_index(gammas, g)
            j = find_index(ps, p)
            d_st, a_st = int(lookup[i,j,0]), float(lookup[i,j,1])
            F_st, Y_st = static_purification(g, p, d_st, a_st)
            R_st = (2**d_st)/cycle_time * Y_st

            gh, ph, probes = estimate_channel(g, p, N_probe, M_probe)
            ii = find_index(gammas, gh)
            jj = find_index(ps, ph)
            d_ad, a_ad = int(lookup[ii,jj,0]), float(lookup[ii,jj,1])
            F_ad, Y_ad = static_purification(g, p, d_ad, a_ad)
            R_ad = (2**d_ad)/cycle_time * Y_ad - probes/cycle_time

            diff_results[(g, p)] = {
                'ΔF':    F_ad - F_st,
                'ΔY':    Y_ad - Y_st,
                'ΔThroughput': R_ad - R_st
            }

    # Plot surfaces
    plot_diff_surface(diff_results, gammas, ps, metric='ΔF', title='Fidelity Difference')
    plot_diff_surface(diff_results, gammas, ps, metric='ΔY', title='Yield Difference')
    plot_diff_surface(diff_results, gammas, ps, metric='ΔThroughput', title='Throughput Difference')
    # Plot contour maps
    plot_diff_contour(diff_results, gammas, ps, metric='ΔF', title='Fidelity Difference Contour')
    plot_diff_contour(diff_results, gammas, ps, metric='ΔY', title='Yield Difference Contour')
    plot_diff_contour(diff_results, gammas, ps, metric='ΔThroughput', title='Throughput Difference Contour')

    print(f"Simulation completed in {time.time()-start:.2f} seconds.")

if __name__=='__main__':
    np.random.seed(42)
    simulate()
