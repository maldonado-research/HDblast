# Independent continuation reconstruction

The stdlib-only producer reconstructs four exact source ledgers from the pinned
old interval reports and checks the shared forced-energy, scalar-history and
partial null-constraint consequences of two separately supplied candidate
reports. It makes no scientific field queries and imports no project module.
Its 45 named gates are finite algebra/authentication/comparison gates, with
the analytic limitations explained in [PROOF.md](PROOF.md).

Run once from a repository checkout, replacing the absolute paths below:

```sh
python -I -S -B RESEARCH/independent/run_controls.py \
  --repo-root /ABSOLUTE/CHECKOUT \
  --source-candidate /ABSOLUTE/SOURCE_REPORT.json \
  --theory-candidate /ABSOLUTE/THEORY_REPORT.json \
  --output /ABSOLUTE/FRESH/CONTROLS.json
```

The wrapper runs the canonical normal/optimized pair once and eight actual
rejection children. Their reports are retained in CONTROLS.artifacts. All
inputs and the output argument must be absolute; output must be fresh and
outside the repository/component. Every child uses -I -S -B, adding -O for
optimized execution. The four rejection cases are changed old source bytes,
an altered candidate source budget, a malformed 1/0 candidate fraction and
an occupied output. Runtime failures produce one JSON line, empty stderr and
exit 2, with fresh output absent or existing bytes preserved.

The retained reference pair has SHA256
02a22609ca2f02e2c2445ae23e6c0c6ad1f03f411f4c609dce648bca0f8ee9da.
It binds source candidate 29043299ad58ac8233493b5f6da7a8f88c03894d60ddfb19cc5f914e2de53e09
and theory candidate 544932f53fc858e4fcafc4b278281643a9088fe9035aeaa441664d7b803d7ab5.
Internal static reviews and reproducible checks are distinct from external
specialist review, which has not been obtained here.
