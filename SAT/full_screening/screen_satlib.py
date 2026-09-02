#!/usr/bin/env python3
"""Resumable affine-FR-only screening of the SATLIB archive collection.

This is an isolated screening harness.  It reads the original .tar.gz files
directly, does not extract into or modify SAT/data, and does not solve an SDP.
The auxiliary LP matches the split formulation in ../../src/affineFR.m.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import tarfile
import time
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import gurobipy as gp
import numpy as np
import scipy
import scipy.sparse as sp


DEFAULT_ARCHIVE_DIR = Path(__file__).resolve().parent / "archives"
DEFAULT_RESULTS_DIR = Path(__file__).resolve().parent / "results"


@dataclass
class ScreenRecord:
    archive: str
    family: str
    member: str
    instance_id: str
    num_vars: int | str
    declared_clauses: int | str
    parsed_clauses: int | str
    eligible: bool | str
    skip_reason: str
    solver_status: str
    solver_objective: float | str
    support_count: int | str
    has_nonzero_reduction: bool | str
    primal_residual: float | str
    solver_runtime_seconds: float | str
    wall_runtime_seconds: float
    message: str


def parse_dimacs(data: bytes, source: str) -> tuple[int, int, list[list[int]]]:
    """Parse DIMACS CNF as a zero-terminated token stream."""
    text = data.decode("ascii", errors="strict")
    num_vars: int | None = None
    declared_clauses: int | None = None
    clauses: list[list[int]] = []
    current_clause: list[int] = []

    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if not line or line.startswith("c"):
            continue
        # Several SATLIB archives use '%' as an end-of-instance marker.
        if line.startswith("%"):
            if current_clause:
                raise ValueError(f"unterminated clause before % in {source}")
            break
        parts = line.split()
        if parts[0] == "p":
            if len(parts) < 4 or parts[1].lower() != "cnf":
                raise ValueError(f"invalid problem line at {source}:{line_number}")
            if num_vars is not None:
                raise ValueError(f"multiple problem lines in {source}")
            num_vars = int(parts[2])
            declared_clauses = int(parts[3])
            if num_vars < 0 or declared_clauses < 0:
                raise ValueError(f"negative header count in {source}")
            continue
        if num_vars is None:
            raise ValueError(
                f"clause data precedes the problem line at {source}:{line_number}"
            )
        for token in parts:
            literal = int(token)
            if literal == 0:
                # Some legacy SATLIB files append a standalone zero after all
                # header-declared clauses.  It is an end marker, not an extra
                # empty clause.  A zero encountered before the declared count
                # is still preserved as a genuine empty clause.
                if current_clause or len(clauses) < declared_clauses:
                    clauses.append(current_clause)
                    current_clause = []
            else:
                if len(clauses) >= declared_clauses and not current_clause:
                    raise ValueError(
                        f"nonzero data follows all declared clauses at "
                        f"{source}:{line_number}"
                    )
                if abs(literal) > num_vars:
                    raise ValueError(
                        f"literal {literal} exceeds {num_vars} variables "
                        f"at {source}:{line_number}"
                    )
                current_clause.append(literal)

    if num_vars is None or declared_clauses is None:
        raise ValueError(f"missing DIMACS problem line in {source}")
    if current_clause:
        raise ValueError(f"unterminated clause at end of {source}")
    if len(clauses) != declared_clauses:
        raise ValueError(
            f"clause-count mismatch in {source}: header={declared_clauses}, "
            f"parsed={len(clauses)}"
        )
    return num_vars, declared_clauses, clauses


def build_dual_matrix(num_vars: int, clauses: list[list[int]]) -> sp.csr_matrix:
    """Build C=[A_ineq'; b_ineq'] for the split auxiliary LP.

    The SAT model has clause inequalities A*x >= b and bounds 0 <= x <= 1.
    FRAformat.m converts these to A_ineq*x <= b_ineq before affineFR.m
    constructs C.
    """
    num_clauses = len(clauses)
    num_inequalities = num_clauses + 2 * num_vars
    rows: list[int] = []
    cols: list[int] = []
    values: list[float] = []

    for clause_index, clause in enumerate(clauses):
        negative_count = 0
        for literal in clause:
            variable_index = abs(literal) - 1
            # Original A coefficient is sign(literal); FRAformat negates A.
            coefficient = -1.0 if literal > 0 else 1.0
            rows.append(variable_index)
            cols.append(clause_index)
            values.append(coefficient)
            negative_count += literal < 0
        # Original rhs is 1-negative_count; FRAformat negates it.
        b_ineq = float(negative_count - 1)
        if b_ineq:
            rows.append(num_vars)
            cols.append(clause_index)
            values.append(b_ineq)

    for variable_index in range(num_vars):
        lower_col = num_clauses + variable_index
        upper_col = num_clauses + num_vars + variable_index
        # -x_i <= 0
        rows.append(variable_index)
        cols.append(lower_col)
        values.append(-1.0)
        # x_i <= 1
        rows.extend((variable_index, num_vars))
        cols.extend((upper_col, upper_col))
        values.extend((1.0, 1.0))

    matrix = sp.coo_matrix(
        (values, (rows, cols)),
        shape=(num_vars + 1, num_inequalities),
        dtype=float,
    ).tocsr()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    return matrix


def status_name(status: int) -> str:
    names = {
        gp.GRB.LOADED: "LOADED",
        gp.GRB.OPTIMAL: "OPTIMAL",
        gp.GRB.INFEASIBLE: "INFEASIBLE",
        gp.GRB.INF_OR_UNBD: "INF_OR_UNBD",
        gp.GRB.UNBOUNDED: "UNBOUNDED",
        gp.GRB.CUTOFF: "CUTOFF",
        gp.GRB.ITERATION_LIMIT: "ITERATION_LIMIT",
        gp.GRB.NODE_LIMIT: "NODE_LIMIT",
        gp.GRB.TIME_LIMIT: "TIME_LIMIT",
        gp.GRB.SOLUTION_LIMIT: "SOLUTION_LIMIT",
        gp.GRB.INTERRUPTED: "INTERRUPTED",
        gp.GRB.NUMERIC: "NUMERIC",
        gp.GRB.SUBOPTIMAL: "SUBOPTIMAL",
        gp.GRB.INPROGRESS: "INPROGRESS",
        gp.GRB.USER_OBJ_LIMIT: "USER_OBJ_LIMIT",
        gp.GRB.WORK_LIMIT: "WORK_LIMIT",
        gp.GRB.MEM_LIMIT: "MEM_LIMIT",
    }
    return names.get(status, f"STATUS_{status}")


def solve_affine_fr_screen(
    num_vars: int,
    clauses: list[list[int]],
    time_limit: float,
    threads: int,
) -> dict[str, object]:
    """Solve the split LP and decide whether its exposing matrix is nonzero."""
    c_matrix = build_dual_matrix(num_vars, clauses)
    num_inequalities = c_matrix.shape[1]

    model = gp.Model("satlib_affine_fr_screen")
    model.Params.OutputFlag = 0
    model.Params.TimeLimit = time_limit
    model.Params.Threads = threads
    u = model.addMVar(num_inequalities, lb=0.0, ub=1.0, name="u")
    v = model.addMVar(num_inequalities, lb=0.0, name="v")
    model.setObjective(u.sum(), gp.GRB.MAXIMIZE)
    model.addConstr(c_matrix @ u + c_matrix @ v == np.zeros(num_vars + 1))
    model.optimize()

    if model.Status == gp.GRB.INFEASIBLE:
        model.Params.NumericFocus = 3
        model.optimize()

    result: dict[str, object] = {
        "status": status_name(model.Status),
        "objective": "",
        "support_count": "",
        "has_nonzero_reduction": "",
        "primal_residual": "",
        "runtime": model.Runtime,
    }
    if model.Status not in (gp.GRB.OPTIMAL, gp.GRB.SUBOPTIMAL):
        return result

    u_value = np.asarray(u.X)
    v_value = np.asarray(v.X)
    y_value = u_value + v_value
    equality_residual = np.linalg.norm(c_matrix @ y_value)
    bound_residual = math.sqrt(
        np.linalg.norm(np.minimum(u_value, 0.0)) ** 2
        + np.linalg.norm(np.minimum(v_value, 0.0)) ** 2
        + np.linalg.norm(np.maximum(u_value - 1.0, 0.0)) ** 2
    )
    primal_residual = math.hypot(equality_residual, bound_residual)
    selected = np.flatnonzero(y_value > 0.5)
    if selected.size:
        nonzero_selected_columns = c_matrix[:, selected].getnnz(axis=0) > 0
        has_nonzero_reduction = bool(np.any(nonzero_selected_columns))
    else:
        has_nonzero_reduction = False

    result.update(
        objective=float(model.ObjVal),
        support_count=int(selected.size),
        has_nonzero_reduction=has_nonzero_reduction,
        primal_residual=float(primal_residual),
    )
    if primal_residual >= 1e-5:
        result["status"] = "BIG_RESIDUAL"
        result["has_nonzero_reduction"] = ""
    return result


def archive_family(path: Path) -> str:
    name = path.name
    return name[:-7] if name.endswith(".tar.gz") else path.stem


def instance_members(archive: tarfile.TarFile) -> list[tarfile.TarInfo]:
    return sorted(
        (
            member
            for member in archive.getmembers()
            if member.isfile() and member.name.lower().endswith((".cnf", ".sat"))
        ),
        key=lambda member: member.name,
    )


def load_processed(csv_path: Path, retry_errors: bool) -> set[str]:
    if not csv_path.exists():
        return set()
    with csv_path.open(newline="", encoding="utf-8") as stream:
        return {
            row["instance_id"]
            for row in csv.DictReader(stream)
            if not retry_errors or row["solver_status"] != "ERROR"
        }


def append_record(csv_path: Path, record: ScreenRecord) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    exists = csv_path.exists()
    with csv_path.open("a", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(asdict(record)))
        if not exists:
            writer.writeheader()
        writer.writerow(asdict(record))
        stream.flush()


def write_run_metadata(
    results_dir: Path, archive_dir: Path, args: argparse.Namespace
) -> None:
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "archive_dir": str(archive_dir),
        "results_dir": str(results_dir),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "gurobi": ".".join(map(str, gp.gurobi.version())),
        "max_vars": args.max_vars,
        "max_clauses": args.max_clauses,
        "time_limit": args.time_limit,
        "threads": args.threads,
        "method": (
            "Correct DIMACS token-stream parsing followed by the split "
            "affine-FR auxiliary LP; no SDP solve"
        ),
    }
    (results_dir / "run_metadata.json").write_text(
        json.dumps(metadata, indent=2) + "\n", encoding="utf-8"
    )


def summarize(csv_path: Path, output_path: Path) -> dict[str, object]:
    rows: list[dict[str, str]] = []
    if csv_path.exists():
        with csv_path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
    # Retried instances are appended to preserve the audit history.  For the
    # current summary, retain only the most recent record for each instance.
    latest_rows: dict[str, dict[str, str]] = {}
    for row in rows:
        latest_rows[row["instance_id"]] = row
    rows = list(latest_rows.values())
    by_family: dict[str, dict[str, int]] = {}
    for row in rows:
        family = row["family"]
        counts = by_family.setdefault(
            family,
            {
                "processed": 0,
                "eligible": 0,
                "reduced": 0,
                "skipped_size": 0,
                "errors_or_incomplete": 0,
            },
        )
        counts["processed"] += 1
        if row["eligible"].lower() == "true":
            counts["eligible"] += 1
        if row["has_nonzero_reduction"].lower() == "true":
            counts["reduced"] += 1
        if row["skip_reason"] == "SIZE_LIMIT":
            counts["skipped_size"] += 1
        if row["solver_status"] not in ("OPTIMAL", "SUBOPTIMAL", "SKIPPED_SIZE"):
            counts["errors_or_incomplete"] += 1
    summary = {
        "updated_utc": datetime.now(timezone.utc).isoformat(),
        "total_records": len(rows),
        "families_with_reduction": sorted(
            family for family, counts in by_family.items() if counts["reduced"] > 0
        ),
        "by_family": dict(sorted(by_family.items())),
    }
    output_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def choose_archives(archive_dir: Path, names: Iterable[str]) -> list[Path]:
    all_archives = sorted(archive_dir.glob("*.tar.gz"), key=lambda path: path.name)
    requested = list(names)
    if not requested:
        return all_archives
    lookup = {archive_family(path): path for path in all_archives}
    lookup.update({path.name: path for path in all_archives})
    missing = [name for name in requested if name not in lookup]
    if missing:
        raise ValueError(f"unknown archives: {', '.join(missing)}")
    return [lookup[name] for name in requested]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive-dir", type=Path, default=DEFAULT_ARCHIVE_DIR)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--archives", nargs="*", default=[])
    parser.add_argument("--max-vars", type=int, default=100_000)
    parser.add_argument("--max-clauses", type=int, default=150_000)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--max-instances", type=int, default=0)
    parser.add_argument(
        "--retry-errors",
        action="store_true",
        help="retry instances whose most recent checkpoint has status ERROR",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    archive_dir = args.archive_dir.resolve()
    results_dir = args.results_dir.resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    csv_path = results_dir / "instance_screening.csv"
    summary_path = results_dir / "screening_summary.json"
    write_run_metadata(results_dir, archive_dir, args)
    processed = load_processed(csv_path, retry_errors=args.retry_errors)
    archives = choose_archives(archive_dir, args.archives)
    if not archives:
        raise FileNotFoundError(
            f"no .tar.gz SATLIB archives found in {archive_dir}; "
            "supply --archive-dir or follow the setup in README.md"
        )
    print(
        f"Screening {len(archives)} archive(s); {len(processed)} existing "
        "checkpoint record(s).",
        flush=True,
    )

    completed_this_run = 0
    for archive_index, archive_path in enumerate(archives, start=1):
        family = archive_family(archive_path)
        with tarfile.open(archive_path, mode="r:gz") as archive:
            members = instance_members(archive)
            print(
                f"[{archive_index}/{len(archives)}] {family}: "
                f"{len(members)} instance(s)",
                flush=True,
            )
            for member_index, member in enumerate(members, start=1):
                instance_id = f"{archive_path.name}::{member.name}"
                if instance_id in processed:
                    continue
                wall_start = time.perf_counter()
                record = ScreenRecord(
                    archive=archive_path.name,
                    family=family,
                    member=member.name,
                    instance_id=instance_id,
                    num_vars="",
                    declared_clauses="",
                    parsed_clauses="",
                    eligible="",
                    skip_reason="",
                    solver_status="",
                    solver_objective="",
                    support_count="",
                    has_nonzero_reduction="",
                    primal_residual="",
                    solver_runtime_seconds="",
                    wall_runtime_seconds=0.0,
                    message="",
                )
                try:
                    extracted = archive.extractfile(member)
                    if extracted is None:
                        raise ValueError("archive member could not be read")
                    data = extracted.read()
                    num_vars, declared_clauses, clauses = parse_dimacs(
                        data, instance_id
                    )
                    record.num_vars = num_vars
                    record.declared_clauses = declared_clauses
                    record.parsed_clauses = len(clauses)
                    if num_vars > args.max_vars or len(clauses) > args.max_clauses:
                        record.eligible = False
                        record.skip_reason = "SIZE_LIMIT"
                        record.solver_status = "SKIPPED_SIZE"
                    else:
                        record.eligible = True
                        result = solve_affine_fr_screen(
                            num_vars,
                            clauses,
                            time_limit=args.time_limit,
                            threads=args.threads,
                        )
                        record.solver_status = str(result["status"])
                        record.solver_objective = result["objective"]
                        record.support_count = result["support_count"]
                        record.has_nonzero_reduction = result[
                            "has_nonzero_reduction"
                        ]
                        record.primal_residual = result["primal_residual"]
                        record.solver_runtime_seconds = result["runtime"]
                except Exception as error:  # Preserve the failure in the audit trail.
                    record.solver_status = "ERROR"
                    record.message = f"{type(error).__name__}: {error}"
                record.wall_runtime_seconds = time.perf_counter() - wall_start
                append_record(csv_path, record)
                processed.add(instance_id)
                completed_this_run += 1
                if member_index % 100 == 0 or member_index == len(members):
                    print(
                        f"  {member_index}/{len(members)}; "
                        f"latest={record.solver_status}; "
                        f"reduction={record.has_nonzero_reduction}",
                        flush=True,
                    )
                if args.max_instances and completed_this_run >= args.max_instances:
                    summary = summarize(csv_path, summary_path)
                    print(
                        f"Pilot limit reached. Total records: "
                        f"{summary['total_records']}",
                        flush=True,
                    )
                    return 0
        summarize(csv_path, summary_path)

    summary = summarize(csv_path, summary_path)
    print(
        f"Completed. Total records: {summary['total_records']}; "
        f"families with reduction: "
        f"{len(summary['families_with_reduction'])}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
