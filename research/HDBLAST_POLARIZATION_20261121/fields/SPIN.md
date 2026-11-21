# Arbitrary holomorphic patterns at high angular wavenumber

A strongly oscillating angular field can change direction while a charged particle moves. Its magnitude over the brief flight interval can still be controlled in the population average. On the reference histories at the 1% work marker this gives a test for arbitrary linear combinations of the spatially spread modes at angular wavenumber 515, including a different spatial pattern in each of the eight preparation pulses.

## Local angular motion

Write d = κ−2. A polynomial P of degree d in the two complex coordinates on the unit three-sphere defines a horizontal curl eigenform through its real and imaginary parts in the fixed Hopf frame. Its pointwise vector magnitude is |P|. Use Haar probability measure for uniform births and the orthonormal basis

$$
Y_j(z)=\sqrt{(d+1){d\choose j}}\,z_1^{d-j}z_2^j,\qquad 0\leq j\leq d .
$$

An arbitrary P is allowed at each of the eight launch times. Relative phases and both transverse components are included. Patterns outside this holomorphic subspace of the three-sphere curl eigenspace constitute additional degrees of freedom.

For a birth point x choose the unique SU(2) matrix g whose first column is x. A rotation of the sphere transforms P into another degree-d polynomial in the same irreducible representation. Write bⱼ(g) for √(d+1) times its jth normalized coefficient. Haar averaging and Schur orthogonality give ‖bⱼ‖₂ = ‖P‖₂ for every j. For any point within angular geodesic distance ℓ of x, the modulus of its second coordinate in the rotated frame is at most sin ℓ and the first modulus is at most one. The triangle inequality pointwise and then in Haar L² gives

$$
\left\|\sup_{d(x,y)\leq\ell}|P(y)|\right\|_{L^2(dx)}
\leq A_d(\ell)\|P\|_2,\qquad
A_d(\ell)=\sum_{j=0}^d\sqrt{{d\choose j}}\sin^j\ell .
$$

The flight interval obeys ℓ < 0.00308 radians here. Rapid phase change along a Hopf fiber can rotate the electric and magnetic vectors. The local estimate allows every phase at each point on the trajectory.

## Bound on momenta

For wavenumber κ let W be the normalized eight-pulse source Gram matrix. Its quadratic norm summed over the angular coefficients is the modal squared launch amplitude. For each field coefficient vector Eκ(t) across pulse times, its dual norm is √(EκᵀW⁻¹Eκ). Cauchy–Schwarz bounds the angular L² norm of the field by the dual norm times the modal amplitude. Define a unit modal impulse envelope for a cohort born at s by
$$
H_\kappa(s)=\frac{Q}{a(T)}\int_s^T a(t)\left(\|E_\kappa(t)\|_{W^{-1}}+
\|B_\kappa(t)\|_{W^{-1}}\right)dt .
$$

Take the largest endpoint dual norm in each piecewise-linear field cell, the latest expansion factor and the earliest allowed birth time in that cell. The magnetic bound permits any direction and speed up to light speed. The existing axial force preserves the magnitude of transverse comoving momentum. Minkowski's inequality over the birth locations and Cauchy–Schwarz over mode amplitudes give
$$
\big\|(|p_\perp|/(mc)-p_*)_+\big\|_{L^2(dx)}\leq\epsilon\,R(s),\qquad
R(s)^2=\sum_\kappa A_{\kappa-2}(\ell_s)^2 H_\kappa(s)^2 .
$$

Here p★ bounds the birth kick and ε² is the total additional emitted energy in units of the original pulse. Markov's inequality gives a fraction at least max(0, 1−[εR(s)/(p−p★)]²) below each threshold p > p★. The bound on the population mean is p★ + ε ⟨R⟩, with the cohort weights taken from the archived production history. For the moving thermal target a boost along the axial direction preserves the transverse distribution. The thermal mean transverse momentum in units of mc follows by integrating the radial Jüttner distribution and averaging sin θ over birth directions,
$$
\langle p_\perp/(mc)\rangle_{\rm th}
=\frac{\pi}{4}\,\frac{e^{-1/\Theta}(2\Theta+6\Theta^2+6\Theta^3)}{K_2(1/\Theta)} .
$$
Θ is the dimensionless rest-frame temperature. The thermal CDF used in the comparison is given in the field methods note.

## Numerical results

Each test allows independent angular coefficients in every pulse. The first row has 8,224 real launch coefficients and the second row has 13,968.

| Curl blocks | ε required by mean, dedicated / shared | Minimum CDF gaps at ε = 0.02, dedicated / shared |
| --- | ---: | ---: |
| 515 | 0.211 / 0.111 | 76.6 / 58.5 percentage points |
| 360 and 515 | 0.0512 / 0.0296 | 21.6 / 6.8 percentage points |

Doubling the radial mesh gives a 74.8-point dedicated gap for wavenumber 515 and 31.1 points for the two-wavenumber family. These changes show the radial sensitivity of short angular wavelengths. The histories support the qualitative exclusion at ε = 0.02 across these refinements. The comparison across all ten sampled wavenumbers 40 to 515 produces a zero guaranteed CDF gap at ε = 0.02 under this bound. Adding lower wavenumbers can supply a much larger kick per unit launch energy.

The result applies to collisionless test particles, the prescribed background, and the eight sampled preparation waveforms. Heating through collisions, a change in the birth history, or collective fields generated by the particles can change the distribution.

```sh
OPENBLAS_NUM_THREADS=1 python fields/spin.py --history history --scan output/scan --output-dir output/fields
```

The command writes SPIN.json for the native and refined histories using the same input traces as the pulse study.
