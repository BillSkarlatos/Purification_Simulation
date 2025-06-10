import numpy as np
import qutip as qt

def amp_damp_kraus(gamma):
    """
    Returns the single-qubit Kraus operators for amplitude damping with probability gamma.
    """
    return [
        qt.Qobj([[1, 0], [0, np.sqrt(1 - gamma)]]),
        qt.Qobj([[0, np.sqrt(gamma)], [0, 0]])
    ]

def phase_damp_kraus(p):
    """
    Returns the single-qubit Kraus operators for phase damping with probability p.
    """
    return [
        np.sqrt(1 - p) * qt.qeye(2),
        np.sqrt(p) * qt.sigmaz()
    ]

def two_qubit_noise(rho, gamma, p):
    """
    Applies amplitude-damping (gamma) and phase-damping (p) noise to a 2-qubit density matrix rho.
    Returns the resulting noisy density matrix.
    """
    # Ensure rho has composite dims [[2,2],[2,2]]
    if rho.dims == [[4], [4]]:
        rho = qt.Qobj(rho.full(), dims=[[2,2], [2,2]])

    K1 = amp_damp_kraus(gamma)
    K2 = phase_damp_kraus(p)
    out = qt.Qobj(np.zeros((4, 4)), dims=[[2,2], [2,2]])
    for k1 in K1:
        for k2 in K2:
            K = qt.tensor(k1, k2)
            out += K * rho * K.dag()
    return out
