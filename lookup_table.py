import numpy as np
import qutip as qt
from noise import two_qubit_noise
from purification import dejmps_purify, apply_filter
from multiprocessing import Pool, cpu_count

# Simulation parameters
gammas   = np.linspace(0.01, 0.20, 20)
ps       = np.linspace(0.01, 0.20, 20)
depths   = [1, 2, 3]
alphas   = np.linspace(0.1, 0.9, 20)
T        = 500
Ftarget  = 0.90

# Bell state projector and reference state
bell_proj  = qt.bell_state('00').proj()
bell_state = qt.bell_state('00')

# Initialize best trackers
best_fid   = np.full((20, 20), -np.inf)
best_yield = np.full((20, 20), -np.inf)
best_d     = np.zeros((20, 20), dtype=np.uint8)
best_alpha = np.zeros((20, 20), dtype=np.float32)

# Single-trial run

def run_trial(gamma, p, d, alpha):
    # initial noisy pair
    rho = two_qubit_noise(bell_proj, gamma, p)
    total_y = 1.0
    for _ in range(d):
        fr, pf = apply_filter(rho, alpha)
        total_y *= pf
        if pf == 0:
            return 0.0, 0.0   # no survivors
        rho2 = two_qubit_noise(bell_proj, gamma, p)
        fr2, pf2 = apply_filter(rho2, alpha)
        total_y *= pf2
        if pf2 == 0:
            return 0.0, 0.0
        purho, psucc = dejmps_purify(fr, fr2)
        total_y *= psucc
        if purho is None or psucc == 0:
            return 0.0, 0.0
        rho = purho

    # use overlap, not qutip.fidelity
    F = float((bell_proj * rho).tr().real)
    return F, total_y

# Task processor for a single (i,j,d,alpha)
def process_task(args):
    i, j, gamma, p_val, d, alpha = args
    f_list, y_list = [], []
    for _ in range(T):
        f, y = run_trial(gamma, p_val, d, alpha)
        f_list.append(f)
        y_list.append(y)
    return (i, j, d, alpha, np.mean(f_list), np.mean(y_list))

if __name__ == '__main__':
    # Prepare tasks
    tasks = [(i, j, gamma, p_val, d, alpha)
             for i, gamma in enumerate(gammas)
             for j, p_val in enumerate(ps)
             for d in depths
             for alpha in alphas]

    # Parallel execution
    with Pool(cpu_count()) as pool:
        for i, j, d, alpha, avg_f, avg_y in pool.imap_unordered(process_task, tasks, chunksize=50):
            # If fidelity target met, choose by yield
            if avg_f >= Ftarget:
                if avg_y > best_yield[i, j]:
                    best_yield[i, j] = avg_y
                    best_fid[i, j]   = avg_f
                    best_d[i, j]     = d
                    best_alpha[i, j] = alpha
            # Otherwise choose by fidelity
            elif avg_f > best_fid[i, j]:
                best_fid[i, j]   = avg_f
                best_yield[i, j] = avg_y
                best_d[i, j]     = d
                best_alpha[i, j] = alpha

    # Stack into lookup array and save
    lookup = np.stack((best_d, best_alpha), axis=2)
    np.save('lookup_table.npy', lookup)
    print("Lookup table generated and saved (parallelized across %d cores)." % cpu_count())
