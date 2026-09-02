# SATLIB affine-FR screening report

## Scope and method

- Source: all 97 original SATLIB `.tar.gz` archives.
- Inventory: 50,306 `.cnf` or `.sat` instances.
- DIMACS parsing: clauses are token sequences terminated by `0`; the parsed
  clause count is checked against the problem-line declaration.
- Size limits: at most 100,000 variables and at most 150,000 clauses.
- Computation: the split auxiliary LP used by `../../src/affineFR.m`, with one
  Gurobi thread and a 300-second time limit per instance.
- Screening test: an instance is marked as reduced when the inequalities
  selected by the auxiliary LP produce a nonzero exposing matrix in
  `affineFR.m`.
- No SDP was formed or solved.

The original archives were read directly. Nothing was extracted into or
changed in the package's curated `SAT/data` folder.

## Overall outcome

| Item | Count |
|---|---:|
| Archive families | 97 |
| Instances inventoried | 50,306 |
| Instances solved to LP optimality | 50,296 |
| Instances with a nonzero affine-FR reduction | 101 |
| Instances solved with no reduction | 50,195 |
| Instances skipped by the size limits | 8 |
| Instances reaching the 300-second limit | 1 |
| Nonconforming DIMACS instances | 1 |

Ten families contain at least one detected reduction:

| Family | Inventoried | LP eligible | Reduction detected | Notes |
|---|---:|---:|---:|---|
| Bejing | 16 | 16 | 6 | Newly identified by the full screening |
| QG | 22 | 22 | 22 | |
| bf | 4 | 4 | 4 | |
| blocksworld | 7 | 7 | 7 | |
| bmc | 13 | 7 | 6 | Six size skips; one LP time limit |
| hanoi | 2 | 2 | 2 | |
| jnh | 50 | 50 | 27 | |
| logistics | 4 | 4 | 4 | |
| parity | 30 | 30 | 15 | |
| ssa | 8 | 8 | 8 | |

The `inductive-inference` family has 41 instances and none has a detected
reduction under the corrected DIMACS parser. Its earlier apparent reductions
were caused by parsing clauses line by line instead of reading zero-terminated
clauses.

## Exceptions

- `bmc.tar.gz::bmc-ibm-4.cnf` reached the 300-second auxiliary-LP time limit.
  Its reduction status is unresolved.
- `dubois.tar.gz::dubois100.cnf` declares 800 clauses, but strict DIMACS
  token-stream parsing finds only 598 terminators. The file has exactly 800
  noncomment data lines, 202 of which omit the required terminal `0`. A
  separate diagnostic treating each data line as one three-literal clause
  solved optimally and found no reduction, but that recovery is not included
  in the primary strict-parser result.
- Six `bmc` instances and two `gcp-large` instances exceed a size limit. They
  were inventoried but not sent to the auxiliary LP.

## Files

- `results/instance_screening.csv`: instance-level records. Retry records are
  retained for auditability; for duplicate `instance_id` values, the last row
  is the final result.
- `results/screening_summary.json`: deduplicated family-level summary.
- `results/run_metadata.json`: environment, limits, and solver versions.
- `screen_satlib.py`: isolated, resumable screening harness.

This is a family-selection screening. It determines whether the auxiliary LP
selects inequalities that make the exposing matrix nonzero. It does not yet
compute the exact reduced matrix order or rerun the SDP experiments required
for a revised numerical table.
