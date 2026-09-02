# Affine facial reduction experiments

This package contains the MATLAB implementation and numerical experiments for
the paper. The MIPLIB and SAT experiments use the same affine facial-reduction
implementation.

## Directory layout

- `src/`: shared affine and partial facial-reduction routines.
- `MIPLIB/`: MIPLIB mixed-binary preprocessing experiments.
- `SAT/`: SAT preprocessing and SDP experiments.
- `SAT/full_screening/`: full SATLIB family-selection screening and its saved
  instance-level summaries.
- `YALMIP-master/`: bundled YALMIP dependency for the SAT SDP experiments.

The shared `affineFR.m` and `FRAcompute.m` are the SAT versions. They support
the optional rank-only mode used for the large SAT screening, while their
default arguments retain the behavior used by the MIPLIB experiments.

## Requirements

- MATLAB (tested with R2024b Update 6)
- Gurobi with its MATLAB interface on the MATLAB path (tested with 11.0.0)
- Mosek for the SAT SDP experiments (tested with 11.0)

The bundled YALMIP version is 20230622. The SAT benchmark instances used by
the supplied scripts are included. The MIPLIB 2017 instance files are not
included because of their size; users must download them separately and place
the required `.mps` files in `MIPLIB/data/` before running a MIPLIB script. See
`MIPLIB/README.md` for the exact setup and verification instructions.

## Quick tests

From the package directory, run the short MIPLIB preprocessing test with:

```matlab
run('MIPLIB/scripts/demo.m')
```

Reproduce the SAT Table 5 screening with:

```matlab
run('SAT/scripts/run_table5.m')
```

The ten families in Table 5 were selected by the full SATLIB screening in
`SAT/full_screening/`. Its saved summaries are included; reproducing that
screening requires downloading the original SATLIB archives separately. See
`SAT/full_screening/README.md`.

Run the complete MIPLIB experiment and generate its CSV and LaTeX tables with:

```matlab
run('MIPLIB/scripts/run_miplib.m')
```

Reproduce the 14 weighted-SAT experiments in Table 6 with:

```matlab
run('SAT/scripts/run_table6.m')
```

See `MIPLIB/README.md` and `SAT/README.md` for experiment-specific details.

## License and third-party material

The original Affine FR code is released under the MIT License; see `LICENSE`.
Bundled third-party software and benchmark data retain their own terms; see
`THIRD_PARTY.md` for attribution and citation information.
