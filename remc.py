# remc.py
# Replica Exchange Monte Carlo (REMC) pour HP 2D, températures linéaires.

from typing import Dict, List, Tuple, Optional
import numpy as np
from hp import Conformation, energy_hp
from moves import MovesEngine, MovesConfig

def metropolis_accept(delta_E: int, T: float, rng: np.random.Generator) -> bool:
    if delta_E <= 0:
        return True
    return rng.random() < np.exp(-float(delta_E) / max(T, 1e-12))


def run_remc(
    hp_seq: str,
    steps: int = 5000,
    n_replicas: int = 8,
    tmin: float = 0.5,
    tmax: float = 3.0,
    exchange_every: int = 10,
    move_mode: str = "hybrid",
    rho: float = 0.5,
    seed: Optional[int] = None,
) -> Dict[str, object]:
    """
    REMC standard:
      - Réplicas avec T linéairement espacées entre [tmin..tmax].
      - 1 tentative de mouvement par réplica et par step.
      - Echanges pair/impair toutes les 'exchange_every' itérations:
           accepte avec p = min(1, exp((1/Tj - 1/Ti) * (Ei - Ej))).
    Trace:
      - energies: énergie courante du réplica à Tmin (référence) par step.
      - best_energies: meilleure énergie rencontrée (tous réplicas) par step.
    """
    steps = int(min(steps, 10_000))
    n_replicas = max(2, int(n_replicas))
    rng = np.random.default_rng(seed)
    moves = MovesEngine(MovesConfig(rho=rho, rng=rng))

    # Températures linéaires
    temps = np.linspace(tmin, tmax, n_replicas, dtype=float)

    # Initialiser réplicas: même conformation de départ (ligne)
    n = len(hp_seq)
    init_coords = [(i, 0) for i in range(n)]
    init_E = energy_hp(init_coords, hp_seq)

    replicas_coords: List[List[Tuple[int, int]]] = [init_coords.copy() for _ in range(n_replicas)]
    replicas_E = [init_E for _ in range(n_replicas)]

    best_E = init_E
    best_coords = init_coords.copy()

    energies = []
    best_energies = []

    even = True  # alternance des paires d'échange

    for step in range(1, steps + 1):
        # Un move par réplica
        for r in range(n_replicas):
            prop = moves.attempt_move(replicas_coords[r], mode=move_mode)
            if prop is None:
                continue
            E_new = energy_hp(prop, hp_seq)
            if metropolis_accept(E_new - replicas_E[r], temps[r], rng):
                replicas_coords[r] = prop
                replicas_E[r] = E_new
                if E_new < best_E:
                    best_E = E_new
                    best_coords = prop.copy()

        # Trace: énergie du réplica "froid" (indice 0)
        energies.append(replicas_E[0])
        best_energies.append(best_E)

        # Echanges
        if step % exchange_every == 0:
            start = 0 if even else 1
            for i in range(start, n_replicas - 1, 2):
                j = i + 1
                Ei, Ej = replicas_E[i], replicas_E[j]
                Ti, Tj = temps[i], temps[j]
                delta = (1.0 / Tj - 1.0 / Ti) * (Ei - Ej)
                acc = 1.0 if delta >= 0 else np.exp(delta)
                if rng.random() < acc:
                    # swap états des réplicas i et j
                    replicas_coords[i], replicas_coords[j] = replicas_coords[j], replicas_coords[i]
                    replicas_E[i], replicas_E[j] = replicas_E[j], replicas_E[i]
            even = not even

    # Choisir meilleur sur tous réplicas pour le "final"
    # (ici on prend le réplica froid comme final pour cohérence de la trace)
    final_coords = replicas_coords[0]
    final_E = replicas_E[0]
    if best_E < final_E:
        final_coords = best_coords.copy()

    return {
        "final_coords": final_coords,
        "best_coords": best_coords,
        "energies": energies,
        "best_energies": best_energies,
        "best_energy": best_E,
    }
