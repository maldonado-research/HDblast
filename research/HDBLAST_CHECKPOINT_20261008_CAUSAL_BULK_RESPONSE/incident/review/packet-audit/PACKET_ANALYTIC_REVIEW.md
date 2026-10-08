# Independent packet and no-capture audit

Date: 2026-10-08. Scope: sections 4–6, their preparation assumptions, and the
bound-mode overlap in section 5 of `physical-model-incident/INCIDENT_SCATTERING_THEOREM.md`.
The version initially reviewed has SHA256
`d3c03b502503b22535c9487d5f9cc0e61fa02d4bd60b970f611af20bc4f42380`.
Sections 4–6 were reread after the source update with SHA256
`2914c9f287f60398bbe36846d30f04a6c7743aab535750534d13b3697b3ae6e5`.
The update explicitly mentions both position and velocity projections and the
insufficiency of setting `q=q_dot=0`; it preserves and clarifies the accepted
conclusion.
The main independent reviewer owns the final source freeze and control replay.

**Finding: the packet decay and no-capture statements are analytically accepted
within their stated hypotheses. No blocking issue was identified.** This is an
audit of an autonomous, positive, quadratic scalar model and manufactured
preparations, not physical-origin validation or authorization for protected
execution. No retained arrays, observed data, source/trajectory callbacks, or
network resources were consulted.

## 1. Conditions needed for the conclusion

The conclusion uses fixed real parameters with `M>0`, `c>0`, `m²>0`,
`0<g²<m²(c+M)`, one open half-space bulk channel, and no time dependence or
additional interaction. For a global three-dimensional brane packet,
`a(k,p)` is smooth with compact support in `k` and `p`, with common bounds
`|k|<=K` and `0<p0<=p<=p1<infinity`. A single fixed `k` instead defines a
Fourier fiber; it is not a finite-energy plane wave on the whole brane.

The packet is the exact continuum-mode superposition given in equation (13).
Zero initial bound occupation is a spectral preparation condition, not a
consequence of an initially vanishing brane coordinate. The source states the
fiber/global distinction, the compact-support assumptions, and the finite-time
preparation caveat explicitly.

## 2. Independent proof of finite energy and local decay

Write `F=(ip-c)D+g²`, `D=m²-M²-p²`, and
`omega=sqrt(|k|²+M²+p²)`. Since

```
|F|²=(g²-cD)²+p²D²>0  for p>0 and g!=0,
```

`R`, `Q` and every derivative needed below are bounded on the compact packet
support. In particular, the point `D=0` causes no singularity. At every fixed
finite `t`, the amplitudes obtained by multiplying `a` by `R`, `Q`, `omega`,
`k`, `p`, and `exp(-it omega)` are smooth and compactly supported. Extend the
incoming and outgoing transverse Fourier amplitudes by zero to the positive
and negative momentum axes. Their support stays away from zero, so the
extensions are smooth. Their inverse Fourier transforms are Schwartz functions
of `(x,y)` on the full space; restriction to `y>=0` has finite bulk energy.
The brane field and its derivatives are Schwartz in `x`. The trace and mixed
boundary energies are finite as well. The same argument in the `p` variable
gives finite half-line fiber energy at fixed `k`.

On the support,

```
partial_p omega = p/omega >= p0/sqrt(K²+M²+p1²) > 0.
```

For a compactly supported smooth `f(k,p,y)`, integrating by parts gives

```
I(t,k,y) = integral f(k,p,y) exp(-it omega(k,p)) dp
         = (it)^(-N) integral L_k^N f(k,p,y) exp(-it omega(k,p)) dp,
L_k f = partial_p [f/(partial_p omega)].
```

There are no endpoint terms. On `0<=y<=Y` and compact `k` support, the integral
of `|L_k^N f|` is uniformly bounded for the amplitudes of `Phi`, `q`, their
first time and spatial derivatives, and their boundary traces. Thus every such
fiber quantity is `O(|t|^(-N))` as `t` tends to either infinity. Plancherel in
`k` proves the same bound for the corresponding norms in `L²(dx)`; integrating
over finite `Y` proves the finite-layer statement. Squaring gives
`E_loc,abs=O(|t|^(-2N))`. The mixed term obeys

```
|g integral q* Phi_0 dx| <= |g| ||q||_2 ||Phi_0||_2,
```

so it obeys the same estimate. This proves decay of the positive majorant in
equation (15) without assigning an unwarranted positive sign to the interaction
energy alone.

No estimate here is uniform as the support approaches `p=0`, as the parameter
`g` approaches zero, or as the compact momentum bounds are removed. The source
does not claim those extensions.

## 3. Energy escapes the finite layer; global energy does not decay

For the incoming phase `-py-t omega`, on `y>=0` at late positive time its
`p` derivative has magnitude at least `y+v_min t`. The reflected component has
the corresponding property on `y>=0` at early negative time. Repeated
integration by parts in the full phase bounds those wrong-direction components
and their energy derivatives by

```
C_N (1+y+|t|)^(-N),
```

uniformly on the compact `k` support. Derivatives of the phase beyond the first
are `O(|t|)`, which is bounded by a constant times `y+v_min|t|`, so the repeated
integration gives the stated bound. Integrating its square over `y>=0` gives
`O(|t|^(1-2N))`. On the opposite half-line the same argument shows that the
surviving incoming or outgoing component approaches its full-line energy norm.
The unwanted cross term tends to zero by Cauchy–Schwarz, and the boundary terms
decay as proved above.

With the declared `1/sqrt(2pi)` transverse Fourier normalization, full-line
Plancherel gives the complex positive-frequency energy

```
E_in  = integral d³k dp omega² |a(k,p)|²,
E_out = integral d³k dp omega² |a(k,p) R(p)|² = E_in.
```

This verifies the normalization and inference in equation (17). The complex
energy is the sum of the energies of the real and imaginary solutions. Each
real solution obeys its own energy conservation law. The calculation proves
escape from every fixed layer and unit reflected energy, never disappearance
of the full bulk-plus-brane energy. A large but finite observation time can
still contain a transient or long-lived resonance.

## 4. Exact orthogonality, including velocity data

Let `s>0` denote a possible bound inverse length, and use the unnormalized bound
mode

```
b=(q_b,Phi_b)=(1, g exp(-sy)/(c+s)),
m²-z_b=g²/(c+s),  s²=M²-z_b.
```

The Hilbert inner product contains both the brane coordinate and the bulk
integral. Omitting the brane coordinate would give an incorrect overlap. With
`N=(ip+c)D-g²`, `R=N/F`, one has

```
(s-ip)F+(s+ip)N = 2ip[(s+c)D-g²].
```

Consequently the overlap is identically

```
<b,psi_p>
 = Q + g/(c+s) [1/(s+ip)+R/(s-ip)]
 = (2ipg/F) [1 + (D-g²/(c+s))/(s²+p²)]
 = 0,
```

because `D=g²/(c+s)-s²-p²`. This is valid at `D=0` too. The integrals are
absolutely convergent against the bound profile, and the packet amplitudes are
integrable, so superposition preserves the exact zero. The velocity overlap
is also zero, since differentiation only multiplies each continuum mode by
`-i omega(k,p)`. Thus both initial bound position and momentum vanish.

The spectral projection of the fixed self-adjoint spatial operator commutes
with its wave evolution. Zero bound position and momentum remain zero at all
times; a pre-existing nonzero bound component persists as a stable oscillation
in a fixed `k` fiber. The no-capture result therefore follows without a damping
approximation or a claim that every localized initial datum is a continuum
packet. When there is no bound mode the condition is vacuous.

An explicit excluded preparation clarifies the issue. At a finite starting
time, set `q=q_dot=0`, `Phi_dot=0`, and choose a nonzero nonnegative smooth bulk
bump supported strictly away from `y=0`. This satisfies the boundary condition
at that time, but its overlap with `b` is proportional to
`integral exp(-sy) Phi(y) dy`, which is nonzero. Multiplying by a smooth finite
energy profile in `x` gives the global example. Moving the bump farther from
the brane makes this overlap small, not identically zero. An eventual bound
oscillation in these data is pre-existing spectral occupation and does not
contradict the theorem. The source's finite-time preparation qualification
correctly excludes it.

## 5. Bounded scientific and quantum interpretation

For the specified coherent displacement of the full coupled ground state, the
mean field follows the classical solution and the connected Gaussian covariance
does not change. Finite-energy continuum data supported above the positive
continuum gap also have the finite one-particle norm required for this coherent
preparation. No bound coherent displacement is introduced. Residual quantum
vacuum fluctuations are not residual coherent brane excitation, and no
unsmeared vacuum energy has been used.

The theorem does not imply thermalization, entropy production, reheating,
derivation of the incident energy, a primordial source, a Big Bang, or a
dimension-specific empirical signal. Nonlinear interactions, extra accessible
channels, time-dependent couplings/backgrounds, or gravitational dynamics need
a new analysis; the present theorem does not rule those models out. These
limits are explicit in the source and are necessary to its acceptance.

## 6. Editorial observation and review boundary

The harmless typo `Equivalently,+a change` in the initial version was corrected
in the reread source. No source edit was made by this reviewer. Final source hashes and exact
manufactured controls belong to the main review; this document supplies the
independent analytic packet audit only.
