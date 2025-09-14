# viz.py
# Visualisations: repliement 2D, courbes d'énergie, comparatifs.

from typing import List, Tuple, Optional
import json
import matplotlib.pyplot as plt
import numpy as np
from hp import Coord

def _compute_bounds(coords: List[Coord], pad: int = 2):
    xs = [x for x, _ in coords]
    ys = [y for _, y in coords]
    xmin, xmax = min(xs) - pad, max(xs) + pad
    ymin, ymax = min(ys) - pad, max(ys) + pad
    return xmin, xmax, ymin, ymax

def plot_fold_2d(coords: List[Coord], hp_seq: str, energy: int, outpath: str):
    plt.figure(figsize=(6, 6))
    xmin, xmax, ymin, ymax = _compute_bounds(coords, pad=2)

    # Grille carrée
    ax = plt.gca()
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlim([xmin - 0.5, xmax + 0.5])
    ax.set_ylim([ymin - 0.5, ymax + 0.5])
    ax.set_xticks(range(xmin, xmax + 1))
    ax.set_yticks(range(ymin, ymax + 1))
    ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.5)

    # Segments gris
    xs = [x for x, _ in coords]
    ys = [y for _, y in coords]
    plt.plot(xs, ys, '-', color='0.6', linewidth=2, zorder=1, label='Chaîne')

    # Points H/P
    Hx, Hy, Px, Py = [], [], [], []
    for (x, y), t in zip(coords, hp_seq):
        if t == "H":
            Hx.append(x); Hy.append(y)
        else:
            Px.append(x); Py.append(y)

    # H = disque noir plein, P = cercle blanc cerclé noir
    plt.scatter(Hx, Hy, s=120, c='black', edgecolors='black', zorder=3, label='H hydrophobe')
    plt.scatter(Px, Py, s=120, facecolors='white', edgecolors='black', zorder=3, label='P polaire')

    plt.title(f"Repliement 2D — Énergie = {energy}")
    # Légende explicite
    handles, labels = ax.get_legend_handles_labels()
    # dédoublonne
    seen = {}
    h2, l2 = [], []
    for h, l in zip(handles, labels):
        if l not in seen:
            seen[l] = True
            h2.append(h); l2.append(l)
    plt.legend(h2, l2, loc='best', frameon=True)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_energy(energies: List[int], outpath: str, title: str = "Énergie (courante) vs steps", label: Optional[str] = None):
    steps = np.arange(1, len(energies) + 1)
    plt.figure(figsize=(7.5, 4.5))
    plt.plot(steps, energies, label=label if label else "énergie")
    plt.xlabel("Step")
    plt.ylabel("Énergie")
    plt.title(title)
    if label:
        plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_best_energy(best_energies: List[int], outpath: str, title: str = "Meilleure énergie (best-so-far) vs steps", label: Optional[str] = None):
    steps = np.arange(1, len(best_energies) + 1)
    plt.figure(figsize=(7.5, 4.5))
    plt.plot(steps, best_energies, label=label if label else "best-so-far")
    plt.xlabel("Step")
    plt.ylabel("Énergie (meilleure rencontrée)")
    plt.title(title)
    if label:
        plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_energy_compare(
    energies_mc: List[int],
    energies_remc: List[int],
    outpath: str,
    title: str = "Comparatif énergie courante (MC vs REMC)",
):
    steps_mc = np.arange(1, len(energies_mc) + 1)
    steps_remc = np.arange(1, len(energies_remc) + 1)
    m = max(len(steps_mc), len(steps_remc))
    plt.figure(figsize=(8, 5))
    plt.plot(steps_mc, energies_mc, label="MC")
    plt.plot(steps_remc, energies_remc, label="REMC")
    plt.xlabel("Step")
    plt.ylabel("Énergie")
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def plot_best_energy_compare(
    best_mc: List[int],
    best_remc: List[int],
    outpath: str,
    title: str = "Comparatif meilleure énergie (best-so-far) — MC vs REMC",
):
    steps_mc = np.arange(1, len(best_mc) + 1)
    steps_remc = np.arange(1, len(best_remc) + 1)
    plt.figure(figsize=(8, 5))
    plt.plot(steps_mc, best_mc, label="MC")
    plt.plot(steps_remc, best_remc, label="REMC")
    plt.xlabel("Step")
    plt.ylabel("Énergie (best-so-far)")
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(outpath, dpi=200)
    plt.close()

def save_best_fold_json(coords: List[Coord], energy: int, outpath: str):
    data = {"energy": int(energy), "coords": [[int(x), int(y)] for (x, y) in coords]}
    with open(outpath, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
