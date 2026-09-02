# Full SATLIB affine-FR screening

This directory contains the screening used to select the ten SAT families in
Table 5. It reads the 97 original SATLIB archives directly and does not modify
or populate the curated `SAT/data` directory.

The screening has two stages:

1. Parse each DIMACS file as a zero-terminated token stream and verify its
   declared variable and clause counts.
2. For every instance within the limits of 100,000 variables and 150,000
   clauses, solve only the split affine-FR auxiliary LP. No SDP is formed or
   solved.

`results/instance_screening.csv` is appended after every instance, making the
run resumable. `results/screening_summary.json` is refreshed after every
archive.

The isolated Python implementation uses Gurobi's matrix API and matches the
split LP in `../../src/affineFR.m` together with the SAT formulation and bound
handling in `../src/SAT_ILP_FEA.m` and `../../src/FRAformat.m`.

## Archive setup

The SATLIB archives are third-party benchmark data and are not included in
this repository. Download the 97 original `.tar.gz` archives from SATLIB and
either place them in `archives/` or pass their directory explicitly with
`--archive-dir`. The script stops with a clear error if it finds no archives.

The screening requires Python 3.9 or later, NumPy, SciPy, `gurobipy`, and a
working Gurobi license. Run the following commands from this directory.

Example pilot:

```sh
python3 screen_satlib.py --archive-dir /path/to/SATLIB/archives \
    --archives aim inductive-inference
```

Full resumable run:

```sh
python3 screen_satlib.py --archive-dir /path/to/SATLIB/archives
```

The completed-run findings and exceptions are recorded in
`SCREENING_REPORT.md`.
