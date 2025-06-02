import numpy as np
import matplotlib.pyplot as plt

def plot_fidelity_surface(results, gammas, ps):
    """
    Plot a 3D surface of average fidelity over (gamma, p) grid.
    
    Inputs:
        results: dict mapping (gamma, p) tuples to {'avg_fidelity': ..., 'avg_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
    """
    G, P = np.meshgrid(gammas, ps)
    Fdata = np.array([[results[(g, p)]['avg_fidelity'] for g in gammas] for p in ps])

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(G, P, Fdata)
    ax.set_xlabel('Gamma (damping)')
    ax.set_ylabel('p (dephasing)')
    ax.set_zlabel('Avg. Fidelity')
    ax.set_title('3D Surface: Average Fidelity')
    plt.tight_layout()
    plt.show()

def plot_yield_surface(results, gammas, ps):
    """
    Plot a 3D surface of average yield (success probability) over (gamma, p) grid.
    
    Inputs:
        results: dict mapping (gamma, p) tuples to {'avg_fidelity': ..., 'avg_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
    """
    G, P = np.meshgrid(gammas, ps)
    Ydata = np.array([[results[(g, p)]['avg_yield'] for g in gammas] for p in ps])

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(G, P, Ydata)
    ax.set_xlabel('Gamma (damping)')
    ax.set_ylabel('p (dephasing)')
    ax.set_zlabel('Avg. Yield')
    ax.set_title('3D Surface: Average Yield')
    plt.tight_layout()
    plt.show()

def plot_fidelity_contour(results, gammas, ps, levels=20):
    """
    Plot a 2D contour map of average fidelity over (gamma, p) grid.
    
    Inputs:
        results: dict mapping (gamma, p) tuples to {'avg_fidelity': ..., 'avg_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        levels: number of contour levels
    """
    G, P = np.meshgrid(gammas, ps)
    Fdata = np.array([[results[(g, p)]['avg_fidelity'] for g in gammas] for p in ps])

    plt.figure(figsize=(8, 6))
    contour = plt.contourf(G, P, Fdata, levels=levels)
    plt.colorbar(contour, label='Avg. Fidelity')
    plt.xlabel('Gamma (damping)')
    plt.ylabel('p (dephasing)')
    plt.title('Contour: Average Fidelity')
    plt.tight_layout()
    plt.show()

def plot_yield_contour(results, gammas, ps, levels=20):
    """
    Plot a 2D contour map of average yield (success probability) over (gamma, p) grid.
    
    Inputs:
        results: dict mapping (gamma, p) tuples to {'avg_fidelity': ..., 'avg_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        levels: number of contour levels
    """
    G, P = np.meshgrid(gammas, ps)
    Ydata = np.array([[results[(g, p)]['avg_yield'] for g in gammas] for p in ps])

    plt.figure(figsize=(8, 6))
    contour = plt.contourf(G, P, Ydata, levels=levels)
    plt.colorbar(contour, label='Avg. Yield')
    plt.xlabel('Gamma (damping)')
    plt.ylabel('p (dephasing)')
    plt.title('Contour: Average Yield')
    plt.tight_layout()
    plt.show()

def plot_fidelity_vs_gamma(results, gammas, ps, fixed_p_index=0):
    """
    Plot a 2D line graph of fidelity vs gamma for a fixed p value.
    
    Inputs:
        results: dict mapping (gamma, p) tuples to {'avg_fidelity': ..., 'avg_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        fixed_p_index: index in ps to fix p for the plot
    """
    p_val = ps[fixed_p_index]
    Fdata = [results[(g, p_val)]['avg_fidelity'] for g in gammas]

    plt.figure(figsize=(8, 6))
    plt.plot(gammas, Fdata, marker='o')
    plt.xlabel('Gamma (damping)')
    plt.ylabel('Avg. Fidelity')
    plt.title(f'Fidelity vs Gamma (p={p_val:.2f})')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_yield_vs_gamma(results, gammas, ps, fixed_p_index=0):
    """
    Plot a 2D line graph of yield vs gamma for a fixed p value.
    
    Inputs:
        results: dict mapping (gamma, p) tuples to {'avg_fidelity': ..., 'avg_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        fixed_p_index: index in ps to fix p for the plot
    """
    p_val = ps[fixed_p_index]
    Ydata = [results[(g, p_val)]['avg_yield'] for g in gammas]

    plt.figure(figsize=(8, 6))
    plt.plot(gammas, Ydata, marker='o')
    plt.xlabel('Gamma (damping)')
    plt.ylabel('Avg. Yield')
    plt.title(f'Yield vs Gamma (p={p_val:.2f})')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_fidelity_vs_p(results, gammas, ps, fixed_gamma_index=0):
    """
    Plot a 2D line graph of fidelity vs p for a fixed gamma value.
    
    Inputs:
        results: dict mapping (gamma, p) tuples to {'avg_fidelity': ..., 'avg_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        fixed_gamma_index: index in gammas to fix gamma for the plot
    """
    gamma_val = gammas[fixed_gamma_index]
    Fdata = [results[(gamma_val, p)]['avg_fidelity'] for p in ps]

    plt.figure(figsize=(8, 6))
    plt.plot(ps, Fdata, marker='o')
    plt.xlabel('p (dephasing)')
    plt.ylabel('Avg. Fidelity')
    plt.title(f'Fidelity vs p (gamma={gamma_val:.2f})')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_yield_vs_p(results, gammas, ps, fixed_gamma_index=0):
    """
    Plot a 2D line graph of yield vs p for a fixed gamma value.
    
    Inputs:
        results: dict mapping (gamma, p) tuples to {'avg_fidelity': ..., 'avg_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        fixed_gamma_index: index in gammas to fix gamma for the plot
    """
    gamma_val = gammas[fixed_gamma_index]
    Ydata = [results[(gamma_val, p)]['avg_yield'] for p in ps]

    plt.figure(figsize=(8, 6))
    plt.plot(ps, Ydata, marker='o')
    plt.xlabel('p (dephasing)')
    plt.ylabel('Avg. Yield')
    plt.title(f'Yield vs p (gamma={gamma_val:.2f})')
    plt.grid(True)
    plt.tight_layout()
    plt.show()


# Difference plotting functions

def plot_diff_surface(diff_results, gammas, ps, metric, title):
    """
    Plot a 3D surface of difference metric (diff_fidelity or diff_yield) over (gamma, p).
    
    Inputs:
        diff_results: dict mapping (gamma, p) to {'diff_fidelity': ..., 'diff_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        metric: 'diff_fidelity' or 'diff_yield'
        title: title for the plot
    """
    G, P = np.meshgrid(gammas, ps)
    Ddata = np.array([[diff_results[(g, p)][metric] for g in gammas] for p in ps])

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(G, P, Ddata)
    ax.set_xlabel('Gamma (damping)')
    ax.set_ylabel('p (dephasing)')
    ax.set_zlabel(metric)
    ax.set_title(title)
    plt.tight_layout()
    plt.show()

def plot_diff_contour(diff_results, gammas, ps, metric, title, levels=20):
    """
    Plot a 2D contour map of difference metric (diff_fidelity or diff_yield).
    
    Inputs:
        diff_results: dict mapping (gamma, p) to {'diff_fidelity': ..., 'diff_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        metric: 'diff_fidelity' or 'diff_yield'
        title: title for the plot
        levels: number of contour levels
    """
    G, P = np.meshgrid(gammas, ps)
    Ddata = np.array([[diff_results[(g, p)][metric] for g in gammas] for p in ps])

    plt.figure(figsize=(8, 6))
    contour = plt.contourf(G, P, Ddata, levels=levels)
    plt.colorbar(contour, label=metric)
    plt.xlabel('Gamma (damping)')
    plt.ylabel('p (dephasing)')
    plt.title(title)
    plt.tight_layout()
    plt.show()

def plot_diff_vs_gamma(diff_results, gammas, ps, metric, fixed_p_index=0, title=''):
    """
    Plot a 2D line graph of difference metric vs gamma for a fixed p.
    
    Inputs:
        diff_results: dict mapping (gamma, p) to {'diff_fidelity': ..., 'diff_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        metric: 'diff_fidelity' or 'diff_yield'
        fixed_p_index: index in ps to fix p
        title: title for the plot
    """
    p_val = ps[fixed_p_index]
    Ddata = [diff_results[(g, p_val)][metric] for g in gammas]

    plt.figure(figsize=(8, 6))
    plt.plot(gammas, Ddata, marker='o')
    plt.xlabel('Gamma (damping)')
    plt.ylabel(metric)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def plot_diff_vs_p(diff_results, gammas, ps, metric, fixed_gamma_index=0, title=''):
    """
    Plot a 2D line graph of difference metric vs p for a fixed gamma.
    
    Inputs:
        diff_results: dict mapping (gamma, p) to {'diff_fidelity': ..., 'diff_yield': ...}
        gammas: 1D array of gamma values
        ps: 1D array of p values
        metric: 'diff_fidelity' or 'diff_yield'
        fixed_gamma_index: index in gammas to fix gamma
        title: title for the plot
    """
    gamma_val = gammas[fixed_gamma_index]
    Ddata = [diff_results[(gamma_val, p)][metric] for p in ps]

    plt.figure(figsize=(8, 6))
    plt.plot(ps, Ddata, marker='o')
    plt.xlabel('p (dephasing)')
    plt.ylabel(metric)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()
    plt.show()
