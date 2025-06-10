import qutip as qt
import numpy as np

def cnot_4qubit(control, target, num_qubits=4):
    """
    Construct a CNOT acting on a num_qubits-qubit system,
    where 'control' and 'target' are the indices (0-based) of the qubits.

    Returns a qutip.Qobj representing that 2^num_qubits × 2^num_qubits unitary.
    """
    # |0><0| and |1><1| on a single qubit:
    P0 = qt.basis(2, 0) * qt.basis(2, 0).dag()
    P1 = qt.basis(2, 1) * qt.basis(2, 1).dag()
    X  = qt.sigmax()

    # We'll build:  (P0_control ⊗ I_rest) + (P1_control ⊗ X_on_target ⊗ I_others)
    ops = []

    for proj, op_on_control in [(P0, qt.qeye(2)), (P1, X)]:
        op_list = []
        for i in range(num_qubits):
            if i == control:
                op_list.append(proj)
            elif i == target:
                # if proj == P1, we want X on the target; if proj == P0, we want I on target
                op_list.append(op_on_control if proj == P1 else qt.qeye(2))
            else:
                op_list.append(qt.qeye(2))
        ops.append(qt.tensor(op_list))

    # Sum the two pieces
    return ops[0] + ops[1]


def dejmps_purify(rho1, rho2):
    """
    Perform one round of DEJMPS purification on two 2-qubit density matrices rho1 and rho2.
    Returns (purified_rho, success_prob).

    Inputs:
        rho1, rho2: qutip.Qobj (4×4), each a noisy Bell pair for qubits [A1,B1] and [A2,B2]
    Outputs:
        purified_rho: qutip.Qobj (2-qubit) if success_prob > 0; else None
        success_prob: float, probability of successful post-selection
    """
    # Combine the two pairs into a 4-qubit system: ordering [A1, B1, A2, B2]
    rho_combined = qt.tensor(rho1, rho2)

    # Build Alice's CNOT (control=0 → target=2) and Bob's CNOT (control=1 → target=3)
    CNOT_A = cnot_4qubit(control=0, target=2, num_qubits=4)
    CNOT_B = cnot_4qubit(control=1, target=3, num_qubits=4)

    # Apply Alice's then Bob's
    rho_after = CNOT_B * (CNOT_A * rho_combined * CNOT_A.dag()) * CNOT_B.dag()

    # Projectors for measuring qubits 2 (A2) and 3 (B2) in the Z-basis
    P0 = qt.basis(2, 0) * qt.basis(2, 0).dag()
    P1 = qt.basis(2, 1) * qt.basis(2, 1).dag()

    # We only keep cases where both qubits 2&3 are 00 or both are 11
    proj_00 = qt.tensor(qt.qeye(2), qt.qeye(2), P0, P0)
    proj_11 = qt.tensor(qt.qeye(2), qt.qeye(2), P1, P1)

    rho_00 = proj_00 * rho_after * proj_00
    p_00   = rho_00.tr()
    rho_11 = proj_11 * rho_after * proj_11
    p_11   = rho_11.tr()

    success_prob = float(p_00 + p_11)
    if success_prob == 0:
        return None, 0.0

    # Normalize the post-selected state
    rho_success = (rho_00 + rho_11) / success_prob

    # Trace out the measured qubits (indices 2 and 3) to keep only qubits [0, 1]
    purified_rho = rho_success.ptrace([0, 1])
    return purified_rho, success_prob

def apply_filter(rho, alpha):
    """
    Apply the non-unitary local filter F(α)⊗F(α) to a 2-qubit state.
    Returns (rho_filtered, success_prob).
    """
    # single-qubit filter
    F = np.diag([np.sqrt(alpha), np.sqrt(1-alpha)])
    F2 = qt.tensor(qt.Qobj(F), qt.Qobj(F))
    rho_f = F2 * rho * F2.dag()
    p_filt = float(rho_f.tr())
    if p_filt == 0:
        return None, 0.0
    return (rho_f / p_filt), p_filt
