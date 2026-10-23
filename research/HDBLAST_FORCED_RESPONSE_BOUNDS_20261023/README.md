# A certified local auxiliary forced response and global norm bounds

For the retained Einstein–scalar shell model, actual-source interval integration
and a bound for every remaining bulk coupling certify a strictly positive yy
probability sphere mean of the auxiliary forced response W at scaled
`(x,r)=(1/4,5/4)`. The certificate also covers the full rectangle
`|x-1/4|,|r-5/4|<=1/1000000`, where

```text
1.7e7 < mean(W_yy) < 6.5e10.
```

This is a coordinate metric mean before shell contact, not a physical force
or an invariant observable. All 23 bulk couplings are retained. The shell
and reflection contributions vanish on this region because its entire past
cone is separated from the shell, rather than by changing the boundary law.
The point enclosure and `|mean(W_yy)|<=3||W||_X` give the global bound

```text
||W||_X > 8.8e6,
||W||_X < 5.228e13.
```

The global norm is a scaled comparison norm containing first derivatives.
An independent weak-test source argument gives the weaker bound
`||W||_X>976562500000/1002109>974000`. These are not force or energy bounds. The
four-term true-operator Neumann sum has truncation tail below1035 in X, but
its coefficients and a numerical waveform are not evaluated by this package.

Read [the local response proof](compute/EXTERIOR_RESPONSE_PROOF.md),
[independent numerical review](numeric-review/NUMERICAL_REVIEW.md), and
[global-bound proof](CONDITIONAL_RESPONSE_BOUNDS.md) before interpreting the
numbers. The exact-area cone partition has twenty cells. Each actual source
enclosure includes the inherited first-field uncertainty. Full-domain waveform
construction and its complete error certificate remain unfinished. W alone is
not asserted to satisfy the separate full physical quadratic constraint system.

With Python3.12, run from any directory:

```sh
python -I -S -B reproduce.py \
  --manifest-sha256 <expected MANIFEST.json SHA256 from the publication record> \
  --output-dir /absolute/path/to/a/new/output-directory
```

The expected manifest hash is supplied outside the archive to avoid a circular
self-pin. The output directory must be new and outside the source tree.

The [technical note](manuscript/FORCED_RESPONSE_BOUNDS.pdf) presents the
calculation; its TeX source and displayed arithmetic are included.

The reproducer authenticates the exact regular-file inventory, then authenticates
and safely extracts all169 files of the pinned previously published geometric
source archive before any numeric import. It runs twelve top-level isolated
children: the actual-source producer, independent numerical reviewer, numerical
controls, displayed manuscript arithmetic, global-bound arithmetic and bound controls, each in normal and
optimized Python. It compares both typed JSON and raw bytes with saved golden
receipts, regenerates the displayed TeX identically, and reauthenticates both
source trees afterward. PDF typesetting is not part of replay. The old archive is an
immutable input dependency; its complete original scientific replay is not run.

Five bound controls reject a changed source-report byte, a dropped retarded
factor1/2, an incorrect derivative weight3, a missing value factor3, and use
of1-q in the weak lower-bound transfer. Nine separate numerical controls guard
source sectors, clipped areas, inherited authentication, normalization,
coupling feedback and the rectangle error ledger. These are fourteen targeted
in-process cases per execution mode, not fourteen separate child processes.
They do not prove every analytic assumption or form an exhaustive corruption
campaign.

No third-party Python dependencies, network calls, journal submission or
Zenodo operation are part of reproduction.
