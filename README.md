# Purification_Simulation

This repository contains a Monte Carlo simulation of adaptive vs. static entanglement purification over noisy fiber‐optic channels (amplitude‐ and phase‐damping).  All of the generated plots live under `images/`.  Below, you will find:

- A list of all images (file names) in `images/` and what each one illustrates.
- Step‐by‐step instructions on how to install dependencies and run the code to reproduce the data and plots.

---

## Images

Below are all of the plots located in the `images/` folder. Each image is embedded so you can see it directly in the README. The relative path to each file is `images/<filename>.png`.

---

#### `fidelity_surface.png`  
*3D surface plot of average purified‐pair fidelity (adaptive protocol) as a function of amplitude‐damping (γ) and phase‐damping (p).*

![3D surface plot of average purified‐pair fidelity](images/3D_Fidelity.png)

---

#### `yield_surface.png`  
*3D surface plot of average purification yield (adaptive protocol) as a function of γ and p.*

![3D surface plot of average purification yield](images/3D_Yield.png)

---

#### `fidelity_contour.png`  
*2D contour map showing average purified‐pair fidelity (adaptive) over the (γ, p) grid.*

![2D contour map of average purified‐pair fidelity](images/Avg_Fidelity.png)

---

#### `yield_contour.png`  
*2D contour map showing average purification yield (adaptive) over the (γ, p) grid.*

![2D contour map of average purification yield](images/Ang_Yield.png)

---

#### `fidelity_vs_gamma.png`  
*Line plot of average purified‐pair fidelity vs. γ (for a fixed p) under the adaptive protocol.*

![Line plot of fidelity vs. gamma](images/Fideliy_vs_gamma_p001.png)

---

#### `yield_vs_gamma.png`  
*Line plot of average purification yield vs. γ (for a fixed p) under the adaptive protocol.*

![Line plot of yield vs. gamma](images/Yield_vs_gamma_p001.png)

---

#### `fidelity_vs_p.png`  
*Line plot of average purified‐pair fidelity vs. p (for a fixed γ) under the adaptive protocol.*

![Line plot of fidelity vs. p](images/Fidelity_vs_p_gamma001.png)

---

#### `yield_vs_p.png`  
*Line plot of average purification yield vs. p (for a fixed γ) under the adaptive protocol.*

![Line plot of yield vs. p](images/Yield_vs_p_gamma001.png)

---

#### `diff_fidelity_surface.png`  
*3D surface plot of the fidelity difference (adaptive minus default/static) over γ and p.*

![3D surface plot of fidelity difference](images/3D_delta_fidelity.png)

---

#### `diff_yield_surface.png`  
*3D surface plot of the yield difference (adaptive minus default/static) over γ and p.*

![3D surface plot of yield difference](images/3D_delta_yield.png)

---

#### `diff_fidelity_contour.png`  
*2D contour map of fidelity difference (adaptive – default/static) over the (γ, p) grid.*

![2D contour map of fidelity difference](images/delts_fidelity.png)

---

#### `diff_yield_contour.png`  
*2D contour map of yield difference (adaptive – default/static) over the (γ, p) grid.*

![2D contour map of yield difference](images/delta_yield.png)

---

#### `diff_fidelity_vs_gamma.png`  
*Line plot of fidelity difference vs. γ (for a fixed p).*

![Line plot of fidelity difference vs. gamma](images/delta_fidelity_vs_gamma.png)

---

#### `diff_yield_vs_gamma.png`  
*Line plot of yield difference vs. γ (for a fixed p).*

![Line plot of yield difference vs. gamma](images/delta_yield_vs_gamma.png)

---

#### `diff_fidelity_vs_p.png`  
*Line plot of fidelity difference vs. p (for a fixed γ).*

![Line plot of fidelity difference vs. p](images/delta_fidelity_vs_p.png)

---

#### `diff_yield_vs_p.png`  
*Line plot of yield difference vs. p (for a fixed γ).*

![Line plot of yield difference vs. p](images/delta_yield_vs_p.png)

---

## How to Run the Code

Below are step‐by‐step instructions to install all dependencies, run the simulation, extract data, and regenerate all plots.

> **Note:** All commands assume you are in the repository root.

### 1. Clone the Repository

```bash
git clone https://github.com/BillSkarlatos/Purification_Simulation.git
cd Purification_Simulation
```

### 2. Install Dependencies

```bash
pip install --upgrade pip
pip install -r deps.txt
```

### 3. Run the code

```bash
python main.py
python throughput_improvement.py
```
