# main.py
# CLI + orchestration:
# - Lecture séquence (FASTA/AA) 
# - Exécution MC, REMC, ou BOTH
# - Visualisations et sorties (images + CSV + JSON)

import argparse
import csv
import os
from datetime import datetime
from typing import Tuple, Optional

from hp import parse_input, energy_hp
from monte_carlo import run_mc
from remc import run_remc
from viz import (
    plot_fold_2d,
    plot_energy,
    plot_best_energy,
    plot_energy_compare,
    plot_best_energy_compare,
    save_best_fold_json,
)

DEFAULT_OUTPUT_DIR = os.path.join("examples", "output")


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def write_trace_csv(out_csv: str, energies, best_energies):
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["step", "energy", "best_energy"])
        for i, (e, b) in enumerate(zip(energies, best_energies), start=1):
            w.writerow([i, e, b])


def make_run_dir() -> str:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    rd = os.path.join(DEFAULT_OUTPUT_DIR, f"run_{ts}")
    ensure_dir(rd)
    return rd


def main():
    parser = argparse.ArgumentParser(
        description="HP 2D folding (MC & REMC)."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Chemin FASTA OU séquence brute AA .",
    )
    parser.add_argument(
        "--algo",
        choices=["mc", "remc", "both"],
        default="mc",
        help="Algorithme à exécuter: MC, REMC, ou les deux avec comparatif.",
    )
    parser.add_argument("--seed", type=int, default=None, help="Graine RNG.")

    # Mouvements
    parser.add_argument(
        "--moves",
        choices=["vshd", "pull", "hybrid"],
        default="hybrid",
        help="Jeu de mouvements: VSHD, Pull ou Hybride.",
    )
    parser.add_argument(
        "--rho",
        type=float,
        default=0.5,
        help="Pour 'hybrid': probabilité d'utiliser Pull (0..1).",
    )

    # MC
    parser.add_argument("--mc-steps", type=int, default=5000, help="Nombre de steps MC (<= 10000).")
    parser.add_argument("--mc-T", type=float, default=1.0, help="Température (MC).")

    # REMC
    parser.add_argument("--remc-steps", type=int, default=5000, help="Nombre de steps REMC (<= 10000).")
    parser.add_argument("--replicas", type=int, default=8, help="Nombre de réplicas.")
    parser.add_argument("--tmin", type=float, default=0.5, help="Température min (REMC).")
    parser.add_argument("--tmax", type=float, default=3.0, help="Température max (REMC).")
    parser.add_argument(
        "--exchange-every", type=int, default=10, help="Périodicité des tentatives d'échange entre réplicas."
    )

    args = parser.parse_args()

    # Conversion AA->HP (strict)
    hp_seq = parse_input(args.input)
    n = len(hp_seq)
    if n < 2:
        raise SystemExit("La séquence doit contenir au moins 2 résidus.")

    run_dir = make_run_dir()

    # --- MC ---
    res_mc = None
    if args.algo in ("mc", "both"):
        res_mc = run_mc(
            hp_seq=hp_seq,
            steps=args.mc_steps,
            temperature=args.mc_T,
            move_mode=args.moves,
            rho=args.rho,
            seed=args.seed,
        )
        # sorties MC
        e_mc = res_mc["energies"]
        b_mc = res_mc["best_energies"]
        bestE_mc = int(res_mc["best_energy"])
        best_coords_mc = res_mc["best_coords"]

        # images
        plot_fold_2d(best_coords_mc, hp_seq, bestE_mc, os.path.join(run_dir, "fold_final_mc.png" if args.algo == "both" else "fold_final.png"))
        plot_energy(e_mc, os.path.join(run_dir, "energy_plot_mc.png" if args.algo == "both" else "energy_plot.png"), title="Énergie courante (MC)", label="MC")
        plot_best_energy(b_mc, os.path.join(run_dir, "best_energy_plot_mc.png" if args.algo == "both" else "best_energy_plot.png"), title="Meilleure énergie (MC)", label="MC")

        # CSV + JSON
        write_trace_csv(os.path.join(run_dir, "energy_trace_mc.csv" if args.algo == "both" else "energy_trace.csv"), e_mc, b_mc)
        save_best_fold_json(best_coords_mc, bestE_mc, os.path.join(run_dir, "best_fold_mc.json" if args.algo == "both" else "best_fold.json"))

    # --- REMC ---
    res_remc = None
    if args.algo in ("remc", "both"):
        res_remc = run_remc(
            hp_seq=hp_seq,
            steps=args.remc_steps,
            n_replicas=args.replicas,
            tmin=args.tmin,
            tmax=args.tmax,
            exchange_every=args.exchange_every,
            move_mode=args.moves,
            rho=args.rho,
            seed=args.seed,
        )
        # sorties REMC
        e_remc = res_remc["energies"]
        b_remc = res_remc["best_energies"]
        bestE_remc = int(res_remc["best_energy"])
        best_coords_remc = res_remc["best_coords"]

        # images
        plot_fold_2d(best_coords_remc, hp_seq, bestE_remc, os.path.join(run_dir, "fold_final_remc.png" if args.algo == "both" else "fold_final.png"))
        plot_energy(e_remc, os.path.join(run_dir, "energy_plot_remc.png" if args.algo == "both" else "energy_plot.png"), title="Énergie courante (REMC)", label="REMC")
        plot_best_energy(b_remc, os.path.join(run_dir, "best_energy_plot_remc.png" if args.algo == "both" else "best_energy_plot.png"), title="Meilleure énergie (REMC)", label="REMC")

        # CSV + JSON
        write_trace_csv(os.path.join(run_dir, "energy_trace_remc.csv" if args.algo == "both" else "energy_trace.csv"), e_remc, b_remc)
        save_best_fold_json(best_coords_remc, bestE_remc, os.path.join(run_dir, "best_fold_remc.json" if args.algo == "both" else "best_fold.json"))

    # --- Comparatifs si BOTH ---
    if args.algo == "both" and res_mc is not None and res_remc is not None:
        plot_energy_compare(
            res_mc["energies"], res_remc["energies"], os.path.join(run_dir, "energy_compare.png")
        )
        plot_best_energy_compare(
            res_mc["best_energies"], res_remc["best_energies"], os.path.join(run_dir, "best_energy_compare.png")
        )

    print(f"Sorties disponibles dans: {run_dir}")


if __name__ == "__main__":
    main()
