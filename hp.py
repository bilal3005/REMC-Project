# hp.py
# Conversion AA -> HP 
# parsing d'entrée (FASTA ou séquence brute), utilitaires.

from dataclasses import dataclass
from typing import List, Tuple, Optional, Iterable
import os
import numpy as np

# Import facultatif : si BioPython n'est pas présent, on gère la séquence brute
try:
    from Bio import SeqIO  # type: ignore
except Exception:  # pragma: no cover
    SeqIO = None  # type: ignore

Coord = Tuple[int, int]

HYDROPHOBIC = set(list("VIFLMCW"))
POLAR = set(list("DEKRHYSTNQGAP"))
VALID_AA = HYDROPHOBIC.union(POLAR)

AA_TO_HP = {aa: "H" for aa in HYDROPHOBIC}
AA_TO_HP.update({aa: "P" for aa in POLAR})


def is_hp_string(seq: str) -> bool:
    s = seq.strip().upper()
    return all(ch in ("H", "P") for ch in s) and len(s) > 0


def convert_aa_to_hp(seq: str) -> str:
    """
    Convertit une séquence d'acides aminés (AA) en HP selon la table stricte.
    Refuse toute lettre en dehors de V,I,F,L,M,C,W (H) et D,E,K,R,H,Y,S,T,N,Q,G,A,P (P).
    """
    s = seq.strip().upper().replace(" ", "")
    if is_hp_string(s):
        raise ValueError("Séquences HP (H/P) interdites en entrée : fournissez une séquence AA (acides aminés).")
    hp = []
    for i, ch in enumerate(s):
        if ch not in VALID_AA:
            raise ValueError(
                f"Lettre AA invalide '{ch}' à la position {i+1}. "
                "Table stricte: H={V,I,F,L,M,C,W} ; P={D,E,K,R,H,Y,S,T,N,Q,G,A,P}."
            )
        hp.append(AA_TO_HP[ch])
    return "".join(hp)


def parse_input(input_arg: str) -> str:
    """
    Accepte soit un chemin vers un fichier FASTA, soit une séquence brute AA.
    """
    # Si le chemin existe, on tente FASTA
    if os.path.exists(input_arg):
        if SeqIO is None:
            raise RuntimeError("BioPython n'est pas installé, impossible de lire un FASTA. "
                               "Installez biopython ou fournissez une séquence brute AA.")
        records = list(SeqIO.parse(input_arg, "fasta"))
        if not records:
            raise ValueError(f"Aucune séquence trouvée dans le FASTA : {input_arg}")
        # Concatène si plusieurs entrées (classique : 1 seule)
        aa_seq = "".join(str(rec.seq) for rec in records)
        return convert_aa_to_hp(aa_seq)
    else:
        # Séquence brute
        return convert_aa_to_hp(input_arg)


def manhattan(p: Coord, q: Coord) -> int:
    return abs(p[0] - q[0]) + abs(p[1] - q[1])


def neighbours4(p: Coord) -> List[Coord]:
    x, y = p
    return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]


def is_self_avoiding(coords: List[Coord]) -> bool:
    """Vérifie qu'il n'y a pas d'auto-intersection et que la chaîne est connectée."""
    if len(set(coords)) != len(coords):
        return False
    for i in range(len(coords) - 1):
        if manhattan(coords[i], coords[i + 1]) != 1:
            return False
    return True


def initial_line(n: int) -> List[Coord]:
    """Conformation initiale simple: ligne horizontale sur x = 0..n-1, y = 0."""
    return [(i, 0) for i in range(n)]


def energy_hp(coords: List[Coord], hp_seq: str) -> int:
    """
    Énergie HP 2D: -1 par contact H-H entre voisins orthogonaux non consécutifs.
    Impl. O(n) via table d'occupation.
    """
    occ = {coords[i]: i for i in range(len(coords))}
    E = 0
    for i, p in enumerate(coords):
        if hp_seq[i] != "H":
            continue
        for nb in neighbours4(p):
            j = occ.get(nb, None)
            if j is None:
                continue
            if hp_seq[j] != "H":
                continue
            if abs(i - j) == 1:
                continue  # voisins le long de la chaîne : pas comptés
            if j > i:
                E -= 1
    return E


@dataclass
class Conformation:
    hp_seq: str
    coords: List[Coord]
    energy: Optional[int] = None

    def compute_energy(self) -> int:
        self.energy = energy_hp(self.coords, self.hp_seq)
        return self.energy

    def copy(self) -> "Conformation":
        return Conformation(self.hp_seq, list(self.coords), self.energy if self.energy is not None else None)
