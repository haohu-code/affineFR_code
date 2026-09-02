# SAT experiments

This directory contains the MATLAB code and SATLIB data used for the SAT
experiments in the paper.

## Directory layout

- `src/`: reusable parsers, model builders, facial-reduction routines, and
  solver wrappers.
- `scripts/`: experiment entry points.
- `data/`: the ten SAT families used in Table 5.
- `results/`: generated MATLAB results, summaries, LaTeX rows, and logs.
- `full_screening/`: Python code and saved results for the 50,306-instance
  SATLIB screening used to select the ten families in Table 5.

The shared affine-FR routines are in `../src/`, and the bundled YALMIP
dependency is in `../YALMIP-master/`.

## Reproduce Table 5

Start MATLAB and run:

```matlab
cd('affineFR_code/SAT')  % from the repository root
run('scripts/run_table5.m')
```

The script processes the ten families in Table 5 and writes the per-instance
results (`.mat`) and family summary (`.csv`) to `results`. The family order
and aggregation rules are encoded in `src/summarize_table5.m`.

Table 5 needs only the reduced matrix order. Therefore, `run_table5` sets
`options.computeV = false`: it solves the auxiliary LP and computes the rank
using sparse QR, but does not construct a facial range vector. The reported
time is the sum of the auxiliary-LP time and sparse-rank computation time.

The timing results depend on the machine, MATLAB release, and Gurobi version.
The matrix orders and reduction counts should be reproducible up to numerical
tolerances.

The complete family-selection screening is documented in
`full_screening/README.md`. The original 97 SATLIB archives are not bundled;
the screening script reads them directly from a user-supplied archive
directory.

## Fixed objectives for Table 6

The exact linear objective vectors used for the 14 weighted-SAT instances in
Table 6 are stored in `data/table6_objectives.mat`. They are indexed by
instance name and can be loaded independently of the historical experiment
results:

```matlab
addpath('src')
weights = load_table6_objective('anomaly');
```

The MAT file also records the objective-generation rule and its provenance.
Table 6 runners load these fixed vectors instead of generating new random
objectives.

## Reproduce Table 6

From the `affineFR_code/SAT` directory, run:

```matlab
run('scripts/run_table6.m')
```

The runner uses the fixed objectives above and checkpoints a timestamped raw
MAT file after every completed instance. It also writes a CSV summary, a MAT
summary, and a LaTeX table to `results/`. The unreduced SDPs for the five
`par8` instances and `ssa0432-003` are deliberately not constructed or
solved because they previously exhausted memory. Their raw status is recorded
as `SKIPPED_OOM_SAFETY`, while the LaTeX table retains the historical
`OOM` label.

The reference outputs reported in the manuscript were generated on
2026-09-02 with MATLAB R2024b Update 6, Gurobi 11.0.0, Mosek 11.0, and the
bundled YALMIP 20230622. Solver times, and bounds returned for numerically
sensitive problems, can vary with the software environment.

To regenerate the CSV and LaTeX outputs from a completed raw file without
solving any optimization problem, use:

```matlab
addpath('src')
load('results/table6_raw_TIMESTAMP.mat','table6Results')
generate_table6_outputs(table6Results,'results','TIMESTAMP')
```

The Table 6 SDP experiments require Mosek and YALMIP; all preprocessing
experiments require Gurobi.
