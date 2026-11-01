# Prepared mean and harmonic gauge component

This component supplies the explicit mean preparation, all five original
initial harmonic amplitudes, and the conditioned pure-second gauge. It
preserves the original physical second jet. The two arithmetic API files
retain exactly the reviewed scientific bytes. The conditioned representative
has a scalar displacement `-p xi`; the original one has zero scalar second
data. The first fields and their quadratic sources are unchanged.

The public destination is
`research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/`, with `preparation/` and
`gauge/` below it. The portable component may be located outside the source
repository. No Git executable, repository HEAD, private file, or network
access is needed. `--repo-root` selects the repository/extracted tree that
contains the retained `research/` sources. The explicit paths and exact
bytes are authenticated before any scientific module or SymPy import and
again after calculation. `ESSENTIAL_PINS.json` authenticates this component's
own essential proof, source-pin, infrastructure, scientific, and API files.
The integrating capsule supplies the outer manifest anchor for that pin file.

Use Python 3.10+ and `sympy==1.14.0` (only the preparation check needs SymPy).
From this gauge directory, choose the retained source root and fresh output
directory outside both the source tree and component:

```sh
repo=/path/to/extracted/HDblast
out=$(mktemp -d)
python -I -B checks/check_preparation.py --repo-root "$repo" --output "$out/preparation.json"
python -I -B checks/check_initial_mean.py --repo-root "$repo" --output "$out/initial.json"
python -I -B checks/check_conditioned.py --repo-root "$repo" --output "$out/conditioned.json"
python -I -B checks/run_controls.py --repo-root "$repo" --output "$out/controls.json"
python -I -B -O checks/check_preparation.py --repo-root "$repo" --output "$out/preparation-o.json"
python -I -B -O checks/check_initial_mean.py --repo-root "$repo" --output "$out/initial-o.json"
python -I -B -O checks/check_conditioned.py --repo-root "$repo" --output "$out/conditioned-o.json"
python -I -B -O checks/run_controls.py --repo-root "$repo" --output "$out/controls-o.json"
```

`--output -` writes the canonical JSON to stdout. File destinations must be
fresh and outside both source/component trees. The shipped `references/`
reports are immutable comparison records, never reproduction output paths.
Normal and optimized results use identical canonical data; optimization
flags, private paths, elapsed times, and temporary directory names are not
included in reports.

The control runner makes 14 real subprocess calls: three positive baselines,
nine scientific mutants, and two authentication mutants. Every negative
case requires its exact final witness: exit code 2, empty stderr, no output
file, and one complete JSON failure line with the expected stable witness
code and no success report or marker. An unrelated failure cannot
satisfy a scientific-mutation control. The science reports separately count
their own in-process symbolic/rational controls; these are not subprocess
replays. The initial gauge check covers every spacetime gauge component;
the conditioned check includes the fresh two-channel monopole contraction
bound. These finite checks do not prove the retained analytic IFT or
nonlinear evolution theorems.

For consumers, authenticate both API files before importing them. The
original callable is `mean_initial.initial_means`; the conditioned callable
is `conditioned_mean_initial.conditioned_initial_means`. Both return `U0`,
`V0=partial_x U`, and `U0_r` in order `(K0,L0,N0,J0,w0)` and use only
arithmetic on supplied exact/interval values. Coordinates satisfy
`dy=-rho dz`, `x=t/ell`, `r=z/ell`; the old gate's inverse-rho prose typo
is explicitly resolved in `MEAN_INITIAL_GAUGE.md`.

This component contains no evaluated evolved mean, sign result, or physical
origin mechanism. The response component separately supplies any certified
evolution/observable conclusion.
