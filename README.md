# HP 2D Protein Folding — MC & REMC

**Student:** Bilal D.

**Class:** M2BI

**Course:** Gestion de projet informatique et programmation avancée

**Teachers:**  Dr. Costas Bouyioukos and Pr. Jean-Christophe Gelly 

**Date:** 09/09/2025

---

This project implements the Hydrophobic–Polar lattice protein folding model in Python.
It features both Monte Carlo (MC) and Replica Exchange Monte Carlo (REMC) algorithms, full VSHD and Pull move sets, and clear visualization outputs.

Note: we restricted ourselves to the 2D lattice version of the HP model. Extending the implementation to 3D would be possible, but was not completed here due to time constraints.

Designed to reproduce the spirit of:  
**Thachuk C, Shmygelska A, Hoos HH (2007). _A replica exchange Monte Carlo algorithm for protein folding in the HP model._ BMC Bioinformatics.** https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-342 

---

## Features

- **Input:** FASTA file _or_ raw amino acid (AA) string.
- **Strict AA→HP mapping:**  
  - Hydrophobic (**H**): `V, I, F, L, M, C, W`  
  - Polar (**P**): `D, E, K, R, H, Y, S, T, N, Q, G, A, P`  
  Any other letter → **error**.
- **Move sets:**  
  - **VSHD**: End / Corner / Crankshaft (2D)  
  - **Pull moves (2D, complete)** with propagation/backtracking  
  - **Hybrid** mode controlled by `--rho` (probability of using Pull vs VSHD)
- **Algorithms:**  
  - **MC:** ≤ 10,000 steps. Each step = attempt a move + Metropolis test.  
  - **REMC:** Replicas at temperatures linearly spaced in `[Tmin..Tmax]`. One move/replica/step. Exchange attempts every `exchange_every` steps between adjacent pairs, accepted with:  
- **Outputs:** images + CSV traces, per-run timestamped directory.
- **Clear plots:** final fold in a square grid, energy vs steps, best-so-far vs steps, and MC/REMC comparisons.

---

## Installation

```bash
git clone <this-repo>
cd <this-repo>
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Dependencies (see `requirements.txt`):
- `numpy`
- `matplotlib`
- `biopython` (only needed if you read FASTA files; raw AA strings work without it)

---

## Usage (CLI)

The main entry point is `main.py`.

### Inputs
- **FASTA file**: `--input examples/sample.fasta`
- **Raw AA sequence**: `--input "MVLSTNQDEKRIHYGAFWLCIVMFSTNQ"`  

### Algorithm selection
- `--algo mc` : Monte Carlo only  
- `--algo remc` : Replica Exchange only  
- `--algo both` : Runs both and produces comparison plots

### Move set
- `--moves vshd` | `--moves pull` | `--moves hybrid` (default)  
- `--rho 0.5` controls the hybrid selection probability of Pull (otherwise VSHD).

### MC parameters
- `--mc-steps 5000` (capped at 10000)
- `--mc-T 1.0`

### REMC parameters
- `--remc-steps 5000` (capped at 10000)
- `--replicas 8`
- `--tmin 0.5` `--tmax 3.0` (linear temperature ladder)
- `--exchange-every 10`



### Examples

```bash
# MC only (hybrid moves, rho=0.6)
python main.py --input examples/sample.fasta --algo mc --mc-steps 5000 --mc-T 1.0 --moves hybrid --rho 0.6

# REMC only (8 replicas, exchanges every 10 steps)
python main.py --input "MVLSTNQDEKRIHYGAFWLCIVMFSTNQ" --algo remc --remc-steps 5000 --replicas 8 --tmin 0.5 --tmax 3.0 --exchange-every 10

# Run both and compare
python main.py --input examples/sample.fasta --algo both --mc-steps 5000 --remc-steps 5000
```

---

## Provided Example Sequences

The `examples/` folder contains four real protein sequences in FASTA format.  
These cover different lengths and structural types, useful for testing MC and REMC performance.


| `examples/sample.fasta` | 27 aa | **Toy protein example** — a small synthetic sequence
| `examples/1YRF.fasta` | 35 aa | **HP35 (Villin headpiece subdomain)** : a small all-α domain
| `examples/P01542.fasta` | 46 aa | **Crambin** — a plant protein, very compact and highly hydrophobic
| `examples/POCG47.fasta` | 76 aa | **Ubiquitin (1–76)** — the entire ubiquitin protein, highly stable
| `examples/P00648.fasta` | 110 aa | **Barnase (1–110)** — a bacterial nuclease fragment, with complex architecture


---


## AA → HP Conversion (Strict)

Input must be **amino acid letters** (uppercase/lowercase accepted).

- **H (Hydrophobic):** `V, I, F, L, M, C, W`  
- **P (Polar):** `D, E, K, R, H, Y, S, T, N, Q, G, A, P`  

Any other letter triggers an error with a clear message pointing to the invalid position.

---

## Model & Energy

- 2D square lattice (self-avoiding walk).  
- **Energy:**  
  \[ E = -\#\{\text{non-consecutive orthogonal H–H contacts}\} \]  
- Immediate neighbors along the chain are **not** counted.  
- Implementation uses an occupancy map for \(O(n)\) energy evaluation per conformation.

---

## 🔀 Move Sets

### VSHD (2D)
- **End:** relocate an end monomer to another free neighbor of the second/penultimate monomer.
- **Corner:** pivot a corner monomer into the other lattice corner bridging its neighbors.
- **Crankshaft:** 180° flip of a 4-monomer “U” motif when geometrically feasible.

### Pull moves (2D, complete)
Implements two-square rotations around the segment \((i, i+1)\) with **propagation/backtracking** along the chain to maintain connectivity and avoid collisions, following standard 2D Pull-move mechanics.

### Hybrid
`--moves hybrid` randomly selects between Pull and VSHD at each attempt, with `--rho` being the probability to try Pull first (fallbacks to VSHD if the chosen move fails).

---

## Algorithms

### Monte Carlo (MC)
- Each step:
  1. Attempt a move (per move set).
  2. Compute \(\Delta E\), accept with Metropolis:  
     \[ p = \min(1, e^{-\Delta E / T}) \]

### Replica Exchange Monte Carlo (REMC)
- `R = --replicas` independent replicas at temperatures \(T_1, \dots, T_R\) linearly spaced from `Tmin` to `Tmax`.
- Each step: **1 move per replica** using Metropolis at its \(T_r\).
- Every `--exchange-every` steps: attempt pairwise exchanges among adjacent replicas, alternating even/odd pairs.  
  Acceptance between replicas \(i\) and \(j\):  
  \[ p = \min\left(1,\ \exp\left[\left(\frac{1}{T_j}-\frac{1}{T_i}\right)\left(E_i - E_j\right)\right]\right) \]

---

## Outputs

All outputs go to a timestamped directory, e.g.:
```
examples/output/run_YYYYmmdd_HHMMSS/
```

### For single algorithm runs (`--algo mc` or `--algo remc`)
- `fold_final.png` — 2D fold.  
- `energy_plot.png` — current energy vs steps.
- `best_energy_plot.png` — best-so-far energy vs steps.
- `energy_trace.csv` — columns: `step,energy,best_energy`.
- `best_fold.json` — best conformation coordinates and energy.

### For `--algo both`
- `fold_final_mc.png`, `energy_plot_mc.png`, `best_energy_plot_mc.png`, `energy_trace_mc.csv`, `best_fold_mc.json`
- `fold_final_remc.png`, `energy_plot_remc.png`, `best_energy_plot_remc.png`, `energy_trace_remc.csv`, `best_fold_remc.json`
- `energy_compare.png` — MC vs REMC current energy.
- `best_energy_compare.png` — MC vs REMC best-so-far.

---

## Code Organization

```
.
├── main.py              # CLI + orchestration
├── hp.py                # AA→HP conversion, input parsing, energy, utilities
├── moves.py             # VSHD + full 2D Pull moves (+ hybrid by rho)
├── monte_carlo.py       # MC algorithm
├── remc.py              # REMC algorithm 
├── viz.py               # plots
├── requirements.txt     
├── examples/
│   ├── sample.fasta     # toy example
│   └── output/          # results go here
```

---

## Tips & Notes

- **FASTA parsing** needs `biopython`. If it's missing, either install it or pass a raw AA sequence.
- Set a **random seed** with `--seed` for reproducible trajectories on the same platform/Python/Numpy versions.
- For REMC, ensure a reasonable **temperature ladder**: too narrow → rare exchanges; too wide → poor low-T exploration. Start with defaults and tune `--replicas`, `--tmin`, `--tmax`.
- Pull moves can be more **disruptive** and help escape traps; VSHD can fine-tune packings. The **hybrid** often works better than either alone.
- Energy is **non-positive** (more negative is better).

---


## Reference

- Thachuk C, Shmygelska A, Hoos HH (2007). _A replica exchange Monte Carlo algorithm for protein folding in the HP model._ **BMC Bioinformatics**.  doi:10.1186/1471-2105-8-342 https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-8-342





