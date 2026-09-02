#!/usr/bin/env python3
"""Compute the family-table inputs for the SATLIB Bejing archive.

This diagnostic is isolated from the MATLAB runner. It reruns the affine-FR auxiliary
LP, records the selected inequality matrix N=[-d B], and saves the six
nonzero cases for an independent sparse-rank calculation.
"""

from __future__ import annotations

import csv
import importlib.util
import sys
import tarfile
import time
from pathlib import Path

import gurobipy as gp
import numpy as np
import scipy.io
from scipy.sparse.csgraph import structural_rank


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "archives" / "Bejing.tar.gz"
OUTPUT_DIR = HERE / "results" / "bejing_rank_matrices"
OUTPUT_CSV = HERE / "results" / "bejing_metrics.csv"


def load_screen_module():
    path = HERE / "screen_satlib.py"
    spec = importlib.util.spec_from_file_location("satlib_screen", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main() -> None:
    screen = load_screen_module()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    rows = []
    with tarfile.open(ARCHIVE, "r:gz") as archive:
        for member in screen.instance_members(archive):
            source = f"{ARCHIVE.name}::{member.name}"
            data = archive.extractfile(member).read()
            n, declared, clauses = screen.parse_dimacs(data, source)
            c_matrix = screen.build_dual_matrix(n, clauses)

            start = time.perf_counter()
            model = gp.Model("bejing_affine_fr_metrics")
            model.Params.OutputFlag = 0
            model.Params.TimeLimit = 300
            model.Params.Threads = 1
            inequalities = c_matrix.shape[1]
            u = model.addMVar(inequalities, lb=0.0, ub=1.0, name="u")
            v = model.addMVar(inequalities, lb=0.0, name="v")
            model.setObjective(u.sum(), gp.GRB.MAXIMIZE)
            model.addConstr(c_matrix @ u + c_matrix @ v == np.zeros(n + 1))
            model.optimize()
            if model.Status != gp.GRB.OPTIMAL:
                raise RuntimeError(f"{source}: Gurobi status {model.Status}")

            selected = np.flatnonzero(np.asarray(u.X) + np.asarray(v.X) > 0.5)
            # Columns of C are rows of [-d B].
            n_matrix = c_matrix[:, selected].T.tocsr()
            nonzero = bool(n_matrix.nnz)
            srank = int(structural_rank(n_matrix)) if nonzero else 0
            mat_name = ""
            if nonzero:
                mat_name = Path(member.name).stem + ".mat"
                scipy.io.savemat(
                    OUTPUT_DIR / mat_name,
                    {"N": n_matrix, "n": n, "declared_clauses": declared},
                    do_compression=True,
                )
            rows.append(
                {
                    "member": member.name,
                    "num_vars": n,
                    "num_clauses": declared,
                    "selected_rows": n_matrix.shape[0],
                    "nonzeros_in_N": n_matrix.nnz,
                    "structural_rank": srank,
                    "has_nonzero_reduction": nonzero,
                    "lp_runtime_seconds": model.Runtime,
                    "wall_runtime_seconds": time.perf_counter() - start,
                    "matrix_file": mat_name,
                }
            )
            print(
                member.name,
                f"selected={n_matrix.shape[0]}",
                f"sprank={srank}",
                f"nonzero={nonzero}",
                flush=True,
            )

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
