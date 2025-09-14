# monte_carlo.py
# Algorithme Monte-Carlo simple (Metropolis) pour le modèle HP 2D.

from typing import Dict, List, Tuple, Optional
import numpy as np
from hp import Conformation, energy_hp
from moves import MovesEngine, MovesConfig

def metropolis_accept(delta_E: int, T: float, rng: np.random.Generator) -> bool:
    if delta_E <= 0:
        return True
    # exp(-ΔE / T)
    return rng.random() < np.exp(-float(delta_E) / max(T, 1e-12))


def run_mc(
    hp_seq: str,
    steps: int = 5000,
    temperature: float = 1.0,
    move_mode: str = "hybrid",
    rho: float = 0.5,
    seed: Optional[int] = None,
) -> Dict[str, object]:
    """
    Exécute MC sur une conformation initiale (ligne).
    steps ≤ 10_000 (sera tronqué si plus).
    Retour:
      {
        'final_coords': List[(x,y)],
        'best_coords': List[(x,y)],
        'energies': List[int],            # énergie courante par step
        'best_energies': List[int],       # best-so-far par step
        'best_energy': int
      }
    """
    steps = int(min(steps, 10_000))
    rng = np.random.default_rng(seed)
    moves = MovesEngine(MovesConfig(rho=rho, rng=rng))

    n = len(hp_seq)
    # Conformation initiale: ligne
    coords = [(i, 0) for i in range(n)]
    E = energy_hp(coords, hp_seq)

    best_coords = coords.copy()
    best_E = E

    energies = []
    best_energies = []

    for t in range(steps):
        proposal = moves.attempt_move(coords, mode=move_mode)
        if proposal is None:
            # pas de move possible -> trace et continue
            energies.append(E)
            best_energies.append(best_E)
            continue

        E_new = energy_hp(proposal, hp_seq)
        if metropolis_accept(E_new - E, temperature, rng):
            coords = proposal
            E = E_new
            if E < best_E:
                best_E = E
                best_coords = coords.copy()

        energies.append(E)
        best_energies.append(best_E)

    return {
        "final_coords": coords,
        "best_coords": best_coords,
        "energies": energies,
        "best_energies": best_energies,
        "best_energy": best_E,
    }
