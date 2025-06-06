# Purification_Simulation

This repository contains a Monte Carlo simulation of adaptive vs. static entanglement purification over noisy fiber‐optic channels (amplitude‐ and phase‐damping).  All of the generated plots live under `images/`.  Below, you will find:

- A list of all images (file names) in `images/` and what each one illustrates.
- Step‐by‐step instructions on how to install dependencies and run the code to reproduce the data and plots.

---

## Images

All of the plots below are located in the `images/` folder.  The relative path to any given image is `images/<filename>.png`.  

1. **`images/fidelity_surface.png`**  
   *3D surface plot of average purified‐pair fidelity (adaptive protocol) as a function of amplitude‐damping (γ) and phase‐damping (p).*

2. **`images/yield_surface.png`**  
   *3D surface plot of average purification yield (adaptive protocol) as a function of γ and p.*

3. **`images/fidelity_contour.png`**  
   *2D contour map showing average purified‐pair fidelity (adaptive) over the (γ, p) grid.*

4. **`images/yield_contour.png`**  
   *2D contour map showing average purification yield (adaptive) over the (γ, p) grid.*

5. **`images/fidelity_vs_gamma.png`**  
   *Line plot of average purified‐pair fidelity vs. γ (for a fixed p) under the adaptive protocol.*

6. **`images/yield_vs_gamma.png`**  
   *Line plot of average purification yield vs. γ (for a fixed p) under the adaptive protocol.*

7. **`images/fidelity_vs_p.png`**  
   *Line plot of average purified‐pair fidelity vs. p (for a fixed γ) under the adaptive protocol.*

8. **`images/yield_vs_p.png`**  
   *Line plot of average purification yield vs. p (for a fixed γ) under the adaptive protocol.*

9. **`images/diff_fidelity_surface.png`**  
   *3D surface plot of the fidelity difference (adaptive minus default/static) over γ and p.*

10. **`images/diff_yield_surface.png`**  
    *3D surface plot of the yield difference (adaptive minus default/static) over γ and p.*

11. **`images/diff_fidelity_contour.png`**  
    *2D contour map of fidelity difference (adaptive – default/static) over the (γ, p) grid.*

12. **`images/diff_yield_contour.png`**  
    *2D contour map of yield difference (adaptive – default/static) over the (γ, p) grid.*

13. **`images/diff_fidelity_vs_gamma.png`**  
    *Line plot of fidelity difference vs. γ (for a fixed p).*

14. **`images/diff_yield_vs_gamma.png`**  
    *Line plot of yield difference vs. γ (for a fixed p).*

15. **`images/diff_fidelity_vs_p.png`**  
    *Line plot of fidelity difference vs. p (for a fixed γ).*

16. **`images/diff_yield_vs_p.png`**  
    *Line plot of yield difference vs. p (for a fixed γ).*

---

## How to Run the Code

Below are step‐by‐step instructions to install all dependencies, run the simulation, extract data, and regenerate all plots.

> **Note:** All commands assume you are in the repository root.

### 1. Clone the Repository

```bash
git clone https://github.com/BillSkarlatos/Purification_Simulation.git
cd Purification_Simulation

pip install --upgrade pip
pip install -r deps.txt

python main.py
```

