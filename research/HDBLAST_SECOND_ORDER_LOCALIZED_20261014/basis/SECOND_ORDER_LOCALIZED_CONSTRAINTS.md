# Localized second-order constraint solvability for commuting shell quadrupoles

Ricardo Maldonado (Independent Researcher, ORCID 0009-0009-3937-6527) — technical
continuation, 13–14 October 2026 (America/Los_Angeles). Prepared with substantial
AI assistance (Anthropic Claude); the author retains responsibility.

## 0. Summary

The 13 October shell-Noether checkpoint proved that ten Lorentz quadratic
(Taub-type) conditions are **necessary** for a first-order shell perturbation to
be tangent to an exact solution family, and that, inside the fixed-profile
quadrupole family, all ten vanish exactly when `[A,B]=0`. It left open whether
those conditions are also **sufficient** for the second-order initial
constraint-and-interface problem, flagging the unrestricted interface adjoint
domain and range as the hard part.

This continuation settles the **initial-data** question for every first-order
state whose data vanish near the shell and the cone vertex — in particular for
every commuting quadrupole of the frozen family — by a localization argument:

1. On any annulus of the time-symmetric `T=0` slice, the local Killing initial
   data (KIDs) are exactly the ten Lorentz generators: four boost lapses
   `N=Ω x_a` and six rotations (Lemma 1, exact symbolic verification, including
   an explicit proof that no spherically symmetric static potential exists).
2. Hence (Theorem 2) compactly supported second-order data exist **iff** the ten
   per-copy local Taub integrals vanish. The solution vanishes near the shell,
   so every second-order interface condition holds trivially with zero
   second-order shell displacement. The shell adjoint-domain problem never
   arises. For fixed reflection parity the per-copy integrals are half the
   global ones, so the ten global tests of the previous checkpoint suffice.
3. For `A=B=diag(1,−1,0,0)` (the focused state of the previous checkpoint) we
   construct the second-order correction **explicitly**: closed-form polynomial
   tensors sphere by sphere, plus one radial quadrature with an explicit
   integrating factor. Exact SymPy verification gives `P_H[w]=S_H` and
   `P_M[w]=S_M` identically, for arbitrary radial profiles and every background
   solving the field equations.
4. Physical corollary: for a moving (velocity) quadrupole packet, the
   boost-weighted monopole of its second-order energy source is a sum of
   squares, strictly positive. No compactly supported correction can then
   leave the spherically symmetric part of the geometry unchanged: the packet
   must deform the monopole metric (a localized backreaction). The scalar
   field alone cannot screen it.

What is **not** established: second-order evolution *through* shell contact
(the second-order junction initial-boundary problem with quadratic interface
sources), third and higher orders, nonlinear existence, stability, physical
energy positivity of the full system, or any cosmological/Big Bang inference.

(Proof sections follow in the checkpoint manuscript.)
