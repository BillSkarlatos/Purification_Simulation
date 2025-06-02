import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
from tqdm import trange

# =============================================================================
# 1. Basic noise + DEJMPS definitions (same as before)
# =============================================================================

def amp_damp_kraus(gamma):
    """Single-qubit Kraus operators for amplitude damping (probability gamma)."""
    return [
        qt.Qobj([[1, 0], [0, np.sqrt(1 - gamma)]]),
        qt.Qobj([[0, np.sqrt(gamma)], [0, 0]])
    ]

def phase_damp_kraus(p):
    """Single-qubit Kraus operators for phase damping (probability p)."""
    return [
        np.sqrt(1 - p) * qt.qeye(2),
        np.sqrt(p) * qt.sigmaz() / 2
    ]

def two_qubit_noise(rho, gamma, p):
    """
    Apply amplitude-damping (gamma) and phase-damping (p) to a 2-qubit state rho.
    Ensures dims = [[2,2],[2,2]].
    """
    if rho.dims == [[4], [4]]:
        rho = qt.Qobj(rho.full(), dims=[[2, 2], [2, 2]])
    K1 = amp_damp_kraus(gamma)
    K2 = phase_damp_kraus(p)
    out = qt.Qobj(np.zeros((4, 4)), dims=[[2, 2], [2, 2]])
    for k1 in K1:
        for k2 in K2:
            K = qt.tensor(k1, k2)
            out += K * rho * K.dag()
    return out

def cnot_4qubit(control, target, num_qubits=4):
    """
    Build a 4-qubit CNOT with `control` and `target` indices.
    """
    P0 = qt.basis(2, 0) * qt.basis(2, 0).dag()
    P1 = qt.basis(2, 1) * qt.basis(2, 1).dag()
    X  = qt.sigmax()
    ops = []
    for proj, op_on_control in [(P0, qt.qeye(2)), (P1, X)]:
        op_list = []
        for i in range(num_qubits):
            if i == control:
                op_list.append(proj)
            elif i == target:
                op_list.append(op_on_control if proj == P1 else qt.qeye(2))
            else:
                op_list.append(qt.qeye(2))
        ops.append(qt.tensor(op_list))
    return ops[0] + ops[1]

def dejmps_purify(rho1, rho2):
    """
    One round of DEJMPS on two 2-qubit pairs (rho1, rho2).
    Returns (purified_2qubit_rho, success_probability).
    """
    rho_combined = qt.tensor(rho1, rho2)  # 4-qubit system
    CNOT_A = cnot_4qubit(0, 2, 4)
    CNOT_B = cnot_4qubit(1, 3, 4)
    # Apply Alice's CNOT then Bob's CNOT
    rho_after = CNOT_B * (CNOT_A * rho_combined * CNOT_A.dag()) * CNOT_B.dag()
    # Project onto |00> or |11> on qubits 2 & 3
    P0 = qt.basis(2, 0) * qt.basis(2, 0).dag()
    P1 = qt.basis(2, 1) * qt.basis(2, 1).dag()
    proj_00 = qt.tensor(qt.qeye(2), qt.qeye(2), P0, P0)
    proj_11 = qt.tensor(qt.qeye(2), qt.qeye(2), P1, P1)
    rho_00 = proj_00 * rho_after * proj_00
    p_00   = rho_00.tr()
    rho_11 = proj_11 * rho_after * proj_11
    p_11   = rho_11.tr()
    success_prob = float(p_00 + p_11)
    if success_prob == 0:
        return None, 0.0
    rho_success = (rho_00 + rho_11) / success_prob
    # Trace out measured qubits 2 & 3
    purified_rho = rho_success.ptrace([0, 1])
    return purified_rho, success_prob

# =============================================================================
# 2. “Adaptive” function (for simplicity, we fix depth=2 here; replace with your lookup logic)
# =============================================================================

def run_adaptive(rho1, rho2, gamma, p, lookup_depth=2):
    """
    Perform `lookup_depth` rounds of purification (this stands in for your real adaptive lookup).
    Returns final fidelity and total yield.
    """
    bell = qt.bell_state('00')
    current_rho = rho1
    total_yield = 1.0
    for _ in range(lookup_depth):
        current_rho, prob = dejmps_purify(current_rho, rho2)
        if prob == 0 or current_rho is None:
            return 0.0, 0.0
        total_yield *= prob
        rho2 = two_qubit_noise(bell.proj(), gamma, p)  # refresh second pair each round
    fidelity = float(qt.fidelity(current_rho, bell))
    return fidelity, total_yield

# =============================================================================
# 3. “Static” function: run exactly `d` rounds of DEJMPS
# =============================================================================

def run_static(rho1, rho2, gamma, p, d):
    """
    Perform exactly d rounds of DEJMPS (no adaptivity).
    Returns (fidelity, yield).
    """
    bell = qt.bell_state('00')
    current_rho = rho1
    total_yield = 1.0
    for _ in range(d):
        current_rho, prob = dejmps_purify(current_rho, rho2)
        if prob == 0 or current_rho is None:
            return 0.0, 0.0
        total_yield *= prob
        rho2 = two_qubit_noise(bell.proj(), gamma, p)
    fidelity = float(qt.fidelity(current_rho, bell))
    return fidelity, total_yield

# =============================================================================
# 4. Main Monte Carlo loop: compute (F_adapt,Y_adapt) and best static (Y_static)
# =============================================================================

def compute_relative_improvement(gammas, ps, depths_static, trials=200):
    """
    For each (gamma,p) pair, compute:
      - (F_adapt,Y_adapt) using a fixed 'lookup_depth' (e.g. 2)
      - For each d in depths_static, (F_d, Y_d)
      - Pick the static depth d* whose F_d* >= F_adapt with maximal Y_d*
      - Compute R = (Y_adapt - Y_d*) / Y_d* (if Y_d*>0), else large positive number
    
    Returns three 2D arrays (over grid):  F_adapt_grid, Y_adapt_grid, R_grid.
    """
    bell = qt.bell_state('00')
    A = len(gammas); B = len(ps)

    F_adapt_grid = np.zeros((B, A))
    Y_adapt_grid = np.zeros((B, A))
    R_grid       = np.zeros((B, A))  # relative improvement in percentage

    # Pre-allocate storage for static results:
    # static_F[d_idx, j, i], static_Y[d_idx, j, i]
    static_F = np.zeros((len(depths_static), B, A))
    static_Y = np.zeros((len(depths_static), B, A))
    
    # 1) First: compute static results for every d in depths_static (one Monte Carlo per depth)
    for di, d in enumerate(depths_static):
        print(f"Running static depth = {d} ...")
        for j, p_val in enumerate(ps):
            for i, gamma in enumerate(gammas):
                fidelities = []
                yields = []
                for _ in range(trials):
                    rho1 = two_qubit_noise(bell.proj(), gamma, p_val)
                    rho2 = two_qubit_noise(bell.proj(), gamma, p_val)
                    f_d, y_d = run_static(rho1, rho2, gamma, p_val, d)
                    fidelities.append(f_d)
                    yields.append(y_d)
                static_F[di, j, i] = np.mean(fidelities)
                static_Y[di, j, i] = np.mean(yields)

    # 2) Next: compute adaptive results and choose best static for each (i,j)
    for j, p_val in enumerate(ps):
        for i, gamma in enumerate(gammas):
            # (a) Adaptive run
            fidelities_ad = []
            yields_ad     = []
            for _ in range(trials):
                rho1 = two_qubit_noise(bell.proj(), gamma, p_val)
                rho2 = two_qubit_noise(bell.proj(), gamma, p_val)
                f_ad, y_ad = run_adaptive(rho1, rho2, gamma, p_val, lookup_depth=2)
                fidelities_ad.append(f_ad)
                yields_ad.append(y_ad)
            F_ad = np.mean(fidelities_ad)
            Y_ad = np.mean(yields_ad)
            F_adapt_grid[j, i] = F_ad
            Y_adapt_grid[j, i] = Y_ad

            # (b) Find among static depths the one(s) meeting F_d >= F_ad
            candidate_yields = []
            for di, d in enumerate(depths_static):
                if static_F[di, j, i] >= F_ad:
                    candidate_yields.append(static_Y[di, j, i])
                else:
                    # If it cannot reach adaptive fidelity, disqualify
                    candidate_yields.append(0.0)

            Y_static_best = max(candidate_yields)
            if Y_static_best > 0:
                R_grid[j, i] = 100.0 * (Y_ad - Y_static_best) / Y_static_best
            else:
                # If no static depth can match the adaptive fidelity,
                # we assign a large positive improvement (e.g. 100%).
                R_grid[j, i] = 100.0

    return F_adapt_grid, Y_adapt_grid, R_grid

# =============================================================================
# 5. Run everything and plot
# =============================================================================

if __name__ == "__main__":
    # Grid of noise parameters
    gammas = np.linspace(0.01, 0.2, 10)
    ps     = np.linspace(0.01, 0.2, 10)
    depths_static = [1, 2, 3]   # try up to 3 rounds statically
    trials = 200               # Monte Carlo trials per depth
    
    # Compute grids
    F_adapt, Y_adapt, R_rel = compute_relative_improvement(gammas, ps, depths_static, trials)

    # 2D contour plot of relative improvement
    G, P = np.meshgrid(gammas, ps)
    plt.figure(figsize=(6, 5))
    contour = plt.contourf(G, P, R_rel, levels=np.linspace(0, 100, 21), cmap='RdYlBu_r')
    plt.colorbar(contour, label='Relative Throughput Improvement (%)')
    plt.xlabel(r'$\gamma$ (amplitude damping)')
    plt.ylabel(r'$p$ (phase damping)')
    plt.title('Relative Throughput Improvement: Adaptive vs Best Static')
    plt.tight_layout()
    plt.show()
