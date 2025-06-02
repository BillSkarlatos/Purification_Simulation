import numpy as np
import qutip as qt
from noise import two_qubit_noise
from purification import dejmps_purify
from tqdm import trange

# Pre-define the Bell state |Φ⁺⟩
bell = qt.bell_state('00')

def simulate_run(gamma, p, depth):
    """
    Simulate one run of entanglement purification under noise parameters (gamma, p).
    Returns (final_fidelity, total_yield).

    Inputs:
        gamma: float, amplitude-damping probability
        p: float, phase-damping probability
        depth: int, number of DEJMPS rounds
    Outputs:
        final_fidelity: float, fidelity of the final 2-qubit state w.r.t. |Φ⁺⟩
        total_yield: float, success probability after 'depth' rounds
    """
    # Generate two noisy Bell pairs
    rho1 = two_qubit_noise(bell.proj(), gamma, p)
    rho2 = two_qubit_noise(bell.proj(), gamma, p)

    total_yield = 1.0
    current_rho = rho1
    for _ in range(depth):
        current_rho, prob = dejmps_purify(current_rho, rho2)
        if prob == 0 or current_rho is None:
            return 0.0, 0.0
        total_yield *= prob

    # Compute fidelity with respect to |Φ⁺⟩
    fidelity = qt.fidelity(current_rho, bell)
    return float(fidelity), float(total_yield)

def sweep_parameters(gammas, ps, depth, trials=100):
    """
    Perform a Monte Carlo sweep over lists of gamma and p values for a fixed purification depth.
    """
    results = {}
    for gamma in gammas:
        for p_val in ps:
            fidelities = []
            yields = []
            for _ in range(trials):
                f, y = simulate_run(gamma, p_val, depth)
                fidelities.append(f)
                yields.append(y)
            results[(gamma, p_val)] = {
                'avg_fidelity': np.mean(fidelities),
                'avg_yield': np.mean(yields)
            }
    return results

if __name__ == '__main__':
    # Example usage:
    gammas = np.linspace(0.01, 0.2, 10)
    ps = np.linspace(0.01, 0.2, 10)
    depth = 2
    trials = 500
    results = sweep_parameters(gammas, ps, depth, trials)
    # You can save 'results' or pass to visualization routines
