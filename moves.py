# moves.py
# Implémentation des mouvements VSHD (End/Corner/Crankshaft) et Pull 2D complets,
# + moteur de mouvements avec mode hybride contrôlé par rho.

from typing import List, Tuple, Optional
from dataclasses import dataclass
import numpy as np
from hp import Coord, neighbours4, is_self_avoiding, manhattan

@dataclass
class MovesConfig:
    rho: float = 0.5  # prob d'utiliser Pull en mode "hybrid"
    rng: Optional[np.random.Generator] = None

class MovesEngine:
    def __init__(self, config: MovesConfig):
        self.config = config
        self.rng = config.rng if config.rng is not None else np.random.default_rng()

    # ---------- API publics -------------

    def attempt_move(self, coords: List[Coord], mode: str = "hybrid") -> Optional[List[Coord]]:
        """
        mode ∈ {'vshd','pull','hybrid'}.
        Retourne une nouvelle liste de coords si succès, sinon None.
        """
        if mode not in ("vshd", "pull", "hybrid"):
            raise ValueError("mode de mouvements invalide (attendu: 'vshd', 'pull', 'hybrid').")

        if mode == "vshd":
            return self._attempt_vshd(coords)
        elif mode == "pull":
            return self._attempt_pull(coords)
        else:  # hybrid
            if self.rng.random() < self.config.rho:
                res = self._attempt_pull(coords)
                if res is not None:
                    return res
                # si pull échoue, tente VSHD
                return self._attempt_vshd(coords)
            else:
                res = self._attempt_vshd(coords)
                if res is not None:
                    return res
                return self._attempt_pull(coords)

    # ---------- VSHD --------------------

    def _attempt_vshd(self, coords: List[Coord]) -> Optional[List[Coord]]:
        """
        Essaie un mouvement parmi {end, corner, crankshaft} choisi aléatoirement.
        """
        n = len(coords)
        occ = set(coords)
        kinds = ["end", "corner", "crankshaft"]
        self.rng.shuffle(kinds)

        for kind in kinds:
            if kind == "end":
                # choisir extrémité 0 ou n-1
                for end_idx in self.rng.permutation([0, n-1]):
                    newc = coords.copy()
                    nb_of = 1 if end_idx == 0 else n-2
                    anchor = coords[nb_of]
                    # Positions libres adjacentes à l'ancre
                    candidates = [p for p in neighbours4(anchor) if p not in occ or p == coords[end_idx]]
                    # Retire la position actuelle de l'extrémité (sinon ça reviendrait au même)
                    candidates = [p for p in candidates if p != coords[end_idx]]
                    self.rng.shuffle(candidates)
                    for cand in candidates:
                        newc[end_idx] = cand
                        if is_self_avoiding(newc):
                            return newc
                # sinon, continue
            elif kind == "corner":
                # pour chaque i ∈ [1..n-2], si corner possible -> déplace i
                order = list(range(1, n-1))
                self.rng.shuffle(order)
                for i in order:
                    im1 = coords[i-1]; ip1 = coords[i+1]
                    # l'intersection des voisins4(im1) ∩ voisins4(ip1) donne le coin potentiel
                    v1 = set(neighbours4(im1))
                    v2 = set(neighbours4(ip1))
                    candidates = list(v1.intersection(v2))
                    # enlever positions occupées sauf la position actuelle de i
                    candidates = [c for c in candidates if c not in occ or c == coords[i]]
                    # enlever aussi ip1 et im1 pour éviter cassure (ils sont tjs adjacents entre eux)
                    candidates = [c for c in candidates if c != im1 and c != ip1 and c != coords[i]]
                    self.rng.shuffle(candidates)
                    for cand in candidates:
                        newc = coords.copy()
                        newc[i] = cand
                        if is_self_avoiding(newc):
                            return newc
            else:  # crankshaft
                # Implémentation 2D classique: si quatre résidus consécutifs k..k+3 forment un motif en "U",
                # on tente une rotation à 180° du "U".
                # Stratégie robuste: pour i=1..n-3 (i = centre gauche), on vérifie les positions autour.
                order = list(range(0, n-3))
                self.rng.shuffle(order)
                for k in order:
                    p0, p1, p2, p3 = coords[k], coords[k+1], coords[k+2], coords[k+3]
                    # Un "U" si p0 adjacent à p1, p1 adjacent à p2, p2 adjacent à p3, et p0 et p3 sont à distance 2
                    if not (manhattan(p0, p1) == manhattan(p1, p2) == manhattan(p2, p3) == 1):
                        continue
                    if manhattan(p0, p3) != 2:
                        continue
                    # Les positions cibles (flip 180°) pour p1 et p2 sont celles opposées par rapport à la ligne (p0 -> p3)
                    # On déduit le vecteur entre p0 et p3 (deux pas orthogonaux): composantes (±2,0) ou (0,±2) ou (±1,±1)? Non: à distance 2 -> soit (2,0),(0,2) ou (1,1).
                    # Pour un "U", p0->p3 doit être (±2,0) ou (0,±2) OU diagonal (±1,±1). Les cas diag sont très rares pour un vrai U, on les ignore.
                    dx = p3[0] - p0[0]
                    dy = p3[1] - p0[1]
                    if (abs(dx), abs(dy)) == (2, 0):
                        # U horizontal ; p1 et p2 devraient occuper les deux cases adjacentes au milieu (entre p0 et p3) en haut ou en bas.
                        midx = (p0[0] + p3[0]) // 2
                        y = p1[1]  # ligne courante
                        # flip vers l'autre côté (y +/- 1)
                        dy_flip = 1 if p1[1] == p0[1] else -1
                        cand1 = (midx, y + dy_flip)
                        cand2 = (midx + 1 if dx > 0 else midx - 1, y + dy_flip)
                    elif (abs(dx), abs(dy)) == (0, 2):
                        # U vertical ; flip gauche/droite
                        midy = (p0[1] + p3[1]) // 2
                        x = p1[0]
                        dx_flip = 1 if p1[0] == p0[0] else -1
                        cand1 = (x + dx_flip, midy)
                        cand2 = (x + dx_flip, midy + 1 if dy > 0 else midy - 1)
                    else:
                        continue  # on ignore les cas diagonaux pour stabilité

                    # Les nouvelles positions pour p1 et p2 (ordre sans importance)
                    targets = [cand1, cand2]
                    if (targets[0] in coords) or (targets[1] in coords):
                        # autorise si les cibles correspondent aux anciennes positions de p1/p2 ; sinon collision
                        pass

                    # Build tentative conformation
                    newc = coords.copy()
                    # Vérifie que les cibles sont libres ou occupées par p1/p2
                    occ = set(coords)
                    if all((t not in occ or t in (p1, p2)) for t in targets):
                        newc[k+1] = targets[0]
                        newc[k+2] = targets[1]
                        if is_self_avoiding(newc):
                            return newc
                # fin crankshaft
        return None

    # ---------- Pull moves (2D complets) --------------------

    def _attempt_pull(self, coords: List[Coord]) -> Optional[List[Coord]]:
        """
        Implémentation 2D :
        - Choisit i ∈ [1..n-2].
        - Détermine les deux carrés possibles (rotations ±90°) autour du segment (i,i+1).
        - Si C occupé par i-1 -> corner move équivalent.
        - Sinon, déplace i->L, i-1->C et propage j = i-2,..,0 vers la position de j+2 (ancienne).
        - Abandon si collision durant la propagation.
        """
        n = len(coords)
        if n < 3:
            return None
        occ0 = set(coords)

        i_candidates = list(range(1, n - 1))
        self.rng.shuffle(i_candidates)

        for i in i_candidates:
            pi = coords[i]
            pip1 = coords[i + 1]

            dx = pip1[0] - pi[0]
            dy = pip1[1] - pi[1]
            if abs(dx) + abs(dy) != 1:
                continue  # chaîne invalide, mais on sécurise

            # Deux carrés potentiels: s1 = (dy,-dx), s2 = (-dy,dx)
            shifts = [(dy, -dx), (-dy, dx)]
            self.rng.shuffle(shifts)

            for sx, sy in shifts:
                L = (pip1[0] + sx, pip1[1] + sy)
                C = (pi[0] + sx, pi[1] + sy)

                # L doit être libre ou être la position de i (rare)
                if L in occ0 and L != pi:
                    continue

                # Cas corner: C occupé par i-1 -> simple move i->L
                if i - 1 >= 0 and C == coords[i - 1]:
                    newc = coords.copy()
                    newc[i] = L
                    if is_self_avoiding(newc):
                        return newc
                    else:
                        continue

                # Sinon, C doit être libre
                if C in occ0:
                    continue

                # Essai de propagation
                newc = coords.copy()
                old = coords.copy()
                occ = set(coords)

                # Applique i->L, i-1->C
                newc[i] = L
                occ.add(L)
                occ.discard(pi)

                if i - 1 >= 0:
                    newc[i - 1] = C
                    occ.add(C)
                    occ.discard(coords[i - 1])

                # Si i-2 est déjà adjacent au nouveau (i-1), on s'arrête
                if i - 2 >= 0 and manhattan(newc[i - 2], newc[i - 1]) == 1:
                    if is_self_avoiding(newc):
                        return newc
                    else:
                        continue

                # Propagation j = i-2 .. 0 : (xj,yj) <- ancienne position de j+2
                # On vérifie les collisions au fur et à mesure.
                valid = True
                j = i - 2
                while j >= 0:
                    target = old[j + 2]
                    if target in occ and target != newc[j]:
                        valid = False
                        break
                    # Libère ancienne position de j et occupe target
                    occ.discard(newc[j])
                    newc[j] = target
                    occ.add(target)
                    # Si ce déplacement rétablit l'adjacence avec j+1, on peut arrêter tôt
                    if j - 1 >= 0 and manhattan(newc[j - 1], newc[j]) == 1:
                        # On a rétabli la connectivité en remontant suffisamment
                        # mais il faut poursuivre pour s'assurer qu'on n'a pas cassé plus loin.
                        pass
                    j -= 1

                if valid and is_self_avoiding(newc):
                    return newc
                # sinon on tente une autre configuration
        return None
