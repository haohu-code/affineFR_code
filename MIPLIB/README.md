# MIPLIB preprocessing experiments

The scripts compare affine FR, partial FR with the nonnegative diagonal cone,
and partial FR with the diagonally dominant cone on MIPLIB mixed-binary
formulations.

## Required MIPLIB data (not included)

The MIPLIB 2017 instance files are not included in this repository because of
their size. Download `collection.zip` from the
[official MIPLIB download page](https://miplib.zib.de/download.html). The
included `prob_list.mat` contains the canonical list of the 332 mixed-binary
instances used in the paper after the stated filtering and exclusions.

Extract and, if necessary, decompress the required instances so that the
expected uncompressed `.mps` files appear directly under `MIPLIB/data/`. Do
not leave them inside an additional nested directory. For example, the file
for `10teams` must have the package-relative path
`MIPLIB/data/10teams.mps`.

From the package directory, verify the installation in MATLAB with:

```matlab
load('MIPLIB/prob_list.mat','prob_list')
dataDir = fullfile('MIPLIB','data');
installed = cellfun(@(p) isfile(fullfile(dataDir,p.name)),prob_list);
fprintf('%d of %d required instances are installed.\n',nnz(installed),numel(installed))
```

The complete dataset is ready when this prints `332 of 332 required instances
are installed.` The supplied runners also perform this check before beginning
an experiment and report the first missing filenames if the setup is
incomplete. MIPLIB remains third-party data and is subject to its own terms;
see `../THIRD_PARTY.md`.

From the package directory, run:

```matlab
run('MIPLIB/scripts/demo.m')
```

This short demonstration processes instances with at most 100 variables.

Run the complete experiment with:

```matlab
run('MIPLIB/scripts/run_miplib.m')
```

The complete experiment processes all listed instances with at most 10,000
variables and can take many hours. It checkpoints the raw `BIPtable` after each
completed instance and writes all outputs to `results/`. The generated outputs
are:

- timestamped raw and summary `.mat` files;
- a reduction-summary CSV;
- a timing-summary CSV;
- a complete per-instance CSV containing all methods and statuses;
- a LaTeX file containing the reduction and timing tables.

To regenerate the CSV and LaTeX tables from a saved raw checkpoint without
rerunning the optimization experiments, use:

```matlab
addpath('MIPLIB/scripts')
load('MIPLIB/results/miplib_raw_TIMESTAMP.mat','BIPtable')
generate_miplib_tables(BIPtable,'MIPLIB/results','TIMESTAMP')
```

`scripts/BIP_test.m` accepts the maximum number of variables and an optional
checkpoint filename. The benchmark list is stored in `prob_list.mat`, and the
separately downloaded MIPLIB files must be stored in `data/`.
