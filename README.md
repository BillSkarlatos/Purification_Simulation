# Adaptive DEJMPS Purification Simulation

This repository contains Python code for simulating and analyzing an **adaptive** DEJMPS entanglement purification protocol under combined amplitude- and phase-damping noise. The adaptive scheme uses a precomputed lookup table to select the optimal purification depth and filter strength based on real-time channel estimates.

---

## Repository Structure

* `noise.py`             : Defines single- and two-qubit noise channels (amplitude & phase damping).
* `purification.py`      : Implements DEJMPS pair purification and local filtering operations.
* `adaptive.py`          : Main simulation driver comparing static vs adaptive purification. Generates 3D surface & contour plots.
* `export_data_updated.py`: Exports simulation metrics (fidelity, yield, throughput, deltas) to CSV (`adaptive_data_updated.csv`).
* `lookup_table.py`      : Generates the static lookup table (`lookup_table.npy`) via Monte Carlo trials over depth and filter parameters.
* `visualize.py`         : Plotting utilities for 3D surfaces and contour maps of differences.

## Dependencies

* Python 3.8+
* [NumPy](https://numpy.org/)
* [QuTiP](https://qutip.org/) (for quantum objects & fidelity)
* [Matplotlib](https://matplotlib.org/) (for plotting)
* [multiprocessing](https://docs.python.org/3/library/multiprocessing.html) (standard library)

Install via pip:

```bash
pip install -r deps.txt
```

---

## Usage

### 1. Generate Lookup Table

This precomputes optimal purification parameters `(d, α)` for each noise grid point.

```bash
python lookup_table.py
# Produces lookup_table.npy
```

### 2. Run Adaptive vs Static Simulation

Compare static vs adaptive performance and generate figures.

```bash
python adaptive.py
# Outputs 3D surface and contour PNGs under `figures/`
```

### 3. Export Metrics to CSV

Extract fidelity, yield, throughput, and deltas for further analysis or table generation.

```bash
python export_data_updated.py
# Creates adaptive_data_updated.csv
```

---

## File Descriptions

* **lookup\_table.py**: Monte Carlo loops over depths \[1,2,3] and filter strengths α∈\[0.1,0.9] to maximize either fidelity target (≥0.9) or fidelity. Saves `lookup_table.npy` of shape (20,20,2).

* **adaptive.py**:

  1. Loads `lookup_table.npy`.
  2. Estimates noise parameters (γ,p) via single-qubit probes.
  3. Looks up optimal `(d,α)` and applies static purification.
  4. Computes adaptive purification and throughput cost.
  5. Stores `ΔF`, `ΔY`, `ΔThroughput` and plots surfaces & contours.

* **export\_data\_updated.py**:

  * Uses the same static/adaptive routines.
  * Writes `adaptive_data_updated.csv` with columns:
    `gamma, p, d_static, alpha_static, F_static, Y_static, R_static, d_adaptive, alpha_adaptive, F_adaptive, Y_adaptive, R_adaptive, delta_F, delta_Y, delta_T`.

* **noise.py**:

  * `amp_damp_kraus(γ)`, `phase_damp_kraus(p)` return Kraus operators.
  * `two_qubit_noise(ρ, γ, p)` applies both.

* **purification.py**:

  * `apply_filter(ρ, α) -> (ρ_filtered, p_success)` implements local filtering.
  * `dejmps_purify(ρ1,ρ2) -> (ρ_purified, p_success)` implements one DEJMPS round.

* **visualize.py**:

  * `plot_diff_surface(...)` and `plot_diff_contour(...)` for generating PNGs.

---

## Future Work

* Extend lookup to deeper purification depths.
* Automate trigger logic in QKD software based on measured \$(\hatγ,\hat p)\$.
* Integrate real-time lookup in hardware testbed.

---
