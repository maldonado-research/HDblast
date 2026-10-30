# Current first master, dynamic shell trace, and physical junctions

30 October 2026. This is an exact implementation-to-proof coefficient join.
The older arbitrary-function boundary proof and smooth dynamic-boundary
construction remain analytic inputs. No numerical trace solution or new
existence theorem is produced here.

## Conventions and source scope

On each regular interior, `y` is outward proper distance and inward scaled
depth satisfies `partial_y=-(ell*rho)^-1 partial_r`, `x=t/ell`. One-copy
copy-even physical residuals are

```
R_ab = K_out,ab - sigma(phi)*gamma_ab/6,
P = n_out(phi) + sigma_prime(phi)/2.
```

There is one common tension in the two-copy problem. Dividing the even
equations into one-copy rows gives the displayed `1/6` and `1/2`; it does
not add another tension. Background values satisfy `sigma0=6a` and
`sigma1=-2p`, where `a=rho_y/rho`, `p=phi0_y`.

Before selected AST execution, the checker authenticates every manifest
member of the geometric-source, interior-master, and regular-shell-spectra
capsules. It repeats this authentication after the calculation. No project
module is imported. It executes only the authenticated current
`frame_amplitudes` function, selected background coefficient expressions
`R=alpha/rho^2`, `B=3alpha/(rho^3p)`, `k=rho(5a/2-h_ratio)`, its alpha
recurrence, and the actual tension polynomial expressions. The retained
dynamic trace law is read from authenticated mathematical source. This is
not a claim to have rerun a numerical trace solver or its existence proof.

The actual tension is `sigma=2Wsuper+delta(1+c+c eta)`,
`Wsuper=1/3+eta^2+eta^3/3`, fixed `delta=1/1000` and fixed rational `c`.
Its second derivative is exactly `sigma2=4+4eta`; put `s=sigma2/2` below.
The algebra initially retains `s` free, then binds its actual polynomial
value. The source `B` multiplying the scalar reconstruction is distinct
from the shell coefficient `Bcal` defined below.

## Direct first physical residuals

The actual longitudinal pair has

```
u = alpha*vartheta*W/rho^2,
h_yy=-4uY, h_ya=0, h_ab=2uY*gamma0_ab,
chi = B*vartheta*(W_r/ell+kW)*Y,
vartheta=cosh(t)^(-3/2), Hubble=tanh(t).
```

The current alpha recurrence gives
`alpha_y/alpha=h_ratio-5a/2`, with `h_ratio=V_prime/p`. Consequently,
for `Pmaster=rho^2*uY`, the scalar reconstruction obeys the exact identity

```
chi = -3*Pmaster_y/(rho^2*p).
```

Expand the unit normal on this first longitudinal path:
`N^-1=1+2epsilon*uY+O(epsilon^2)`. The physical metric first residual is

```
R1_ab = (u_y+2a*u+p*chi_amplitude/3)*Y*gamma0_ab = 0.
```

This checks the single isotropic coefficient multiplying the full
Lorentzian shell metric, with the other metric structures zero; it does
not count ten copies of the same scalar calculation as independent
equations. The physical scalar residual includes the lapse/unit-normal
term:

```
P1_amplitude = chi_amplitude,y + 2p*u + s*chi_amplitude.
```

The derivative is generated directly from exact correlated functions using
`rho_y=a rho`, `p_y=p(h_ratio-4a)`,
`a_y=-rho^-2-p^2/3`,
`h_ratio,y=V_second-h_ratio^2+4a h_ratio`, and the actual alpha recurrence.
No independent interval endpoints, midpoint background, or differentiated
unknown Taylor remainder is substituted.

## Master and dynamic boundary factorization

At the physical shell, set

```
Bcal = h_ratio-4a+s = p_y/p+sigma2/2,
mu^-1 = rho_b*Bcal,
gamma = rho_b*(h_ratio-5a/2) = -k_b,
Htrace = W_r-ell*gamma*W,
A0(t) = 5/2-35*Hubble^2/4,
Dtrace = W_xx-(ell/mu)*Htrace+ell^2*A0*W.
```

The authenticated retained dynamic law is
`y''=(ell/mu)Htrace-ell^2*A0*y`, `y(x)=W(x,0)`, precisely `Dtrace=0`.
Its positive `mu` is supplied by the retained actual branch certificate;
the algebraic factorization only divides by `rho,p,ell`, which are nonzero
on the retained collar. It does not infer the dynamic trace from an
arbitrary Robin datum.

The actual master residual is

```
Rmaster = W_xx-W_rr+ell^2*(rho^2*Q-35*Hubble^2/4)*W,
Q=2h_ratio^2-10a*h_ratio+35a^2/4-V_second-p^2/6.
```

Direct coefficient collection of the physical scalar residual gives

```
P1_amplitude = [3alpha*vartheta/(ell^2*rho_b^4*p_b)]
               *(Rmaster-Dtrace).
```

This identity is checked before either residual is set to zero. The
`W_xx` and Hubble-dependent terms cancel in their difference, as they must:
the scalar normal condition does not acquire an independent tangential
acceleration before the bulk equation is used. The separate constant
`5/2` and actual potential coefficient are essential. The shell
correlation `mu^-1-gamma=rho_b(s-3a/2)` is checked exactly as well.

Therefore the retained master equation and retained dynamic shell trace
jointly imply the physical scalar junction. The master equation alone
does not. The checker includes a formal jet witness `W=Wr=0`,
`Wrr=Wxx=1`: it has `Rmaster=0`, `Dtrace=1` and nonzero `P1` equal to
minus the displayed prefactor. This witness demonstrates the missing
boundary premise; it is not an actual solution of the boundary problem.

## Current harmonic representative and free tangent labels

The same current reconstruction AST is evaluated with generic raw gauge
amplitudes `f,b,d`. At the shell impose only `f=0`, `b_r=d_r=0` and
their tangential derivatives. In particular, `b`, `d`, and `f_r` remain
free. This is the actual shell-preserving vector condition, not a
requirement that its tangential coordinate labels vanish.

The scalar pullback is `chi_harmonic=chi_longitudinal+p*f`. Its normal
derivative contributes `p*f_y`, canceled by the `-p*h_yy/2` term from the
unit normal. For the physical metric rows, the normal derivative of the
`a*f` tangent metric term cancels the normal lapse gauge contribution.
Mixed metric traces vanish with the stated boundary conditions and the
remaining offdiagonal/anisotropic physical normal rows vanish. Thus the
same physical residuals hold for the current harmonic representative.
The checker verifies these cancellations with free tangential values,
and explicitly checks that the current reconstructed metric still
contains nonzero coefficients of `b`, `d` and `f_r` before cancellation.
No gauge-wave PDE is needed for this physical boundary covariance join;
the separately retained gauge-wave equation supplies harmonicity.

## Verification and limitations

The checker uses only the standard library and exact collected rational
Laurent polynomials. Normal and optimized executions produce identical
complete JSON report bytes. Twelve genuinely corrupted fresh subprocesses
per mode must fail a mathematical identity before producing an output:
normal orientation, scalar reconstruction sign, tension normalization,
unit-normal term, gamma sign, mu sign, Htrace sign, boundary constant,
ell scale, potential sign, gauge scalar normal term, and the missing
dynamic-boundary premise. Source authentication failures are not counted
as mathematical mutation rejection.

This strengthens the existing first-field/boundary consistency chain. It
does not construct another master solution, improve numerical errors,
prove second-order evolution or subsidiary uniqueness, establish an exact
finite-amplitude family, demonstrate a physical force, or establish a
higher-dimensional origin of the Big Bang. Existing smoothness, corner
compatibility, background/master/gauge existence and uniqueness arguments
remain separately identified analytic inputs.
