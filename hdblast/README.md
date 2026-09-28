# HDBLAST — start here

**Higher-dimensional blast research program · Ricardo Maldonado**

| Read | For |
|---|---|
| [HDBLAST in plain language (Sept 2026)](HDBLAST_IN_PLAIN_LANGUAGE.md) | Anyone: the idea, what has been found, what has failed, and what comes next |
| [Program map](PROGRAM_MAP.md) | Researchers: model definition, dated chronology, full claims ledger, falsification tests, file index |
| [Zenodo records](ZENODO_RECORDS.md) | The dated public record of every HDBLAST release |
| [Checkpoints](checkpoints/) | Full research packages (reports, code, data and internal AI referee reviews, not external peer review), 16–22 Sept 2026 |
| [Newest research round](../research/HDBLAST_CHECKPOINT_20260927/) | 27 Sept 2026 onward: stability, exact analysis, particle production, constraint control, new mechanisms (in progress; not yet reviewed) |
| [Public guide v20 (July 2026)](guides/HDBLAST_GENERAL_PUBLIC_GUIDE_V20.md) | The earlier guide, written after the v20 correction (archived; partly out of date) |

## The model in one paragraph

Five-dimensional Einstein gravity with one scalar field φ. The bulk is a mirror-symmetric (Z₂)
double copy, with a single shell (brane) that plays the role of our universe. The superpotential is
W(φ) = 1 − φ + φ³/3, the bulk potential is U = ½W′² − ⅔W², and the shell tension is
σ(φ) = 2W + δ(1 + cφ). The registered values δ = 0.001 and c = 0.5975949350280132 were fixed
**before** the dynamics were computed. Junction conditions:
K^μ_ν = (σ/6) δ^μ_ν and n·∂φ = −σ′(φ)/2.

## Status labels

| Label | Meaning |
|---|---|
| **Exact** | Symbolic or exact-arithmetic identity, verified by a saved script |
| **Certified** | Computer-assisted interval proof, conditional on the premises it lists |
| **Numerical** | Floating-point result with convergence or independent-integration evidence; no error certificate |
| **Conditional** | Holds only under stated extra assumptions |
| **Screening** | Early comparison with approximate likelihoods or inputs |
| **Negative** | A tested route that fails |
| **Withdrawn** | Publicly corrected |
| **Open** | Not established |

## Key results at a glance

- **Certified** (conditional on the M462 root enclosure): the registered shell solution exists and has an
  unstable scalar mode with −7.71788 < m²/H² < −7.71786 (growth rate 1.65719 H).
  **Numerical:** it is the only scalar bound state.
- **Exact:** all linearized Einstein components and the junction algebra were verified symbolically
  (Chat 11). The absence of unstable modes in the tensor, vector and special-harmonic sectors was
  shown with exact identities plus certified inputs (Chat 12), with negative controls. This is mode
  stability, not full linear stability.
- **Numerical:** in full nonlinear 5D evolution the shell has two fates, relaxation toward an empty
  de Sitter brane (the endpoint is approached but not reached in the simulations) or reversal and collapse. The 4D effective theory predicted the turnaround
  point (φ_b = −1.942) before the 5D code confirmed it.
- **Negative:** no radiation-dominated (hot) branch in the homogeneous, classical, matter-free roll-off.
  The shell's own motion is not, by itself, the source of a hot Big Bang.
- **Exact / numerical (22 Sept):** the constant φ = 1 endpoint fails the scalar boundary condition.
  A corrected static solution was found (numerical), and a matter extension was proposed with an
  exact energy-exchange law (no particle production computed yet).
- **Withdrawn (14 July 2026):** a claimed gravitational-radiation certification (PHYS-M309) was
  withdrawn after an audit found it tested the wrong wave direction.

The full ledger, with sources for every line, is in the [program map](PROGRAM_MAP.md#4-claims-ledger).
