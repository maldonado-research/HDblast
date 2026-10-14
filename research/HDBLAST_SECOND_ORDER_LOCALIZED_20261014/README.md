# Localized second-order constraint solvability and a monopole identity

Checkpoint `HDBLAST_SECOND_ORDER_LOCALIZED_20261014` continues the
13 October shell-Noether checkpoint. That checkpoint showed the ten Lorentz
(Taub-type) tests are *necessary* for a linear shell perturbation to extend
to second order.

**What is new**

* On an annulus of the time-symmetric slice that avoids the shell and the
  cone vertex, the local symmetry data are exactly four boosts and six
  rotations (Lemma 1). A localization theorem of Delay then makes the ten
  tests *sufficient* for compactly supported second-order initial data. The
  correction vanishes near the shell, so all second-order junction
  conditions hold with zero shell displacement (Theorem 2).
* For commuting quadrupoles, `[A,B] = 0`, such data always exist. For
  `A = B = diag(1,−1,0,0)` they are written in closed form (Section 4).
* An exact identity with the areal-radius weight `Ω r` fixes a monopole
  moment `m` of every compactly supported correction, for the packet written
  in the standard reconstruction gauge. For velocity packets in that gauge
  `m < 0`, so the spherically symmetric part of the correction is neither a
  compactly supported gauge mode nor a scalar correction (Theorem 3). `m`
  changes under first-order gauge changes of the packet, so this is not a
  gauge-invariant physical statement; an independent referee found this, and
  the note was corrected.

**Further continuation in this reviewed publication:**

* [Reflection symmetry](continuation/REFLECTION_FORMAL_INITIAL_DATA.md) removes the ten local compatibility obstructions at every finite initial-constraint order for commuting quadrupoles, yielding a consistent formal sequence.
* [Pre-contact dynamics](continuation/PRECONTACT_SECOND_ORDER_DYNAMICS.md) constructs second-order spacetime coefficients on a strict pre-contact slab, gives exact initial accelerations for a velocity quadrupole, and evaluates a tensor component of the quadratic scalar stress.
* [Exact angular Hamiltonian proof](manuscript/ANGULAR_EXACT_VERIFICATION.md) replaces reliance on finite sample points for the universal tangential correction identity.

**What is not established:** second-order evolution through shell contact,
convergence of the formal series, exact nonlinear existence, stability, energy positivity, and any
cosmological or Big Bang inference.

Read the [technical note](manuscript/SECOND_ORDER_LOCALIZED.md), then
[RESULTS.md](RESULTS.md) and [NEXT_CALCULATIONS.md](NEXT_CALCULATIONS.md).

## Layout

| Path | Contents |
|---|---|
| `basis/` | The original construction (code, results, note) from the companion session, text and JSON copied unchanged, except derived pickle caches are omitted; original upload hashes are retained in `reviews/UPLOADED_MANIFEST.json` |
| `new/` | Lemma 1 kernel analysis, full `l = 0` sector, displacement canonical form, fast exact verification |
| `independent/` | Numerical checker sharing no code with the original engine: 5D background re-derivation, direct numerical differentiation (40-digit evaluation, 70-digit differentiation), exact `ε`-series |
| `manuscript/` | Technical note (Markdown) |
| `continuation/` | Formal reflection proof, pre-contact evolution, independent algebra and receipts |
| `literature/` | Focused primary-source review, including 2025–2026 work, with applicability limits |
| `reviews/` | Independent referee report and the record of how each point was handled |

## Reproduce

The four `basis/results/*.pkl` files are excluded from the git repository by
its `*.pkl` policy. They are in the original standalone ZIP. `reproduce.py`
always regenerates `sources_cc.pkl` from source code, and recreates the
optional full verifier's other pickle inputs from text JSON. Existing pickle
caches are not loaded. The original ZIP manifest is a provenance record;
the repository manifest covers its distributed text and JSON files.


Python 3.9 or later, `sympy==1.14.0`, `mpmath==1.3.0`:

```
python reproduce.py          # complete default set, including independent l = 0 quadrature
python reproduce.py --core   # CI subset: omits only independent/indep_l0.py
python reproduce.py --full   # complete default set plus the original slow exact verifier
```

Each executed step must produce a fresh receipt and pass its value checks.
The replay includes the independent `f₀` comparison, rejects non-finite
numerical residuals, and records failures in `REPRODUCTION_RECEIPT.json`.
Receipts identify the selected mode and every excluded check;
`all_checks_passed` applies only to that selected set. The `--core` subset
does not independently verify the weighted radial integral.

Ricardo Maldonado, Independent Researcher, ORCID 0009-0009-3937-6527.
Prepared with substantial AI assistance (Anthropic Claude for the recovered checkpoint, OpenAI Codex and internal agents for the review and continuation); internal checks
are not external peer review.
