# Exact quadratic registered screen magnetic coefficient

This is a kinematic, registered-coordinate diagnostic. Its observer and radius
prescription is fixed below. No geometric clock transport, gauge-invariance
certificate, full second-order Einstein–matter solution, detector response, or
cosmological-origin conclusion follows from it.

The full scalar quadratic term is **not zero** in general. It reduces exactly to
one tensor in the registered TT space. Set `R=rho cosh(t)`, `theta=tanh(t)`;
all subscripts `t,y` below are ordinary physical-coordinate derivatives. For the
registered scalar frame amplitudes `A,B,C,D,E,L,M`, define

```
F = ((D+4M-L) M_t - 2 M L_t)/(4 rho)
    - E(D+L)/(2R)
    + (C/R-M_y/2) [ R(C_t-theta C)/(2rho) - R E_y/2 + B/2 ].       (1)
```

For `g(lambda)=g0+lambda q+lambda^2 k/2`, with the frozen standing quadrupole
`Y=n0^2-n1^2`, the **entire scalar-quadratic screen tensor** is

```
Q_H[q,q] = (4F/R) (Tplus-Tminus).                              (2)
```

For equal registered second-response TT amplitudes `u` in `k`,

```
Pi_TT H''(0) = h (Tplus-Tminus),
h = 3 u_x/(2 ell rho^2 cosh(t)) + 4F/R.                        (3)
```

For unequal TT amplitudes replace the two linear coefficients by
`+3(uplus)_x/(2ell rho^2 c)` and `-3(uminus)_x/(2ell rho^2 c)`.
Equation (3) refers to the second parameter derivative itself: there is no
additional `1/2` on its linear term. The quadratic expression (1) needs only
`B,C,D,E,L,M,C_t,E_y,L_t,M_t,M_y` and background `rho,c,theta`. The cancellation
of `A` and `a=rho_y/rho` is derived below, not assumed by omitting metric slots.
The API in `screen_coefficient.py` consumes precisely the dictionaries returned
by the inherited `actual_first_jets.actual_jet` function.

## Conventions and complete observer prescription

Use signature `(-++++)`, curvature
`R^a_bcd=partial_c Gamma^a_db-partial_d Gamma^a_cb+...`, and positive time
extrinsic curvature `K=Lie_U gamma/2`. Write the background metric

```
g0=dy^2+rho(y)^2[-dt^2+c(t)^2 kappa_AB dOmega^A dOmega^B].
```

At every parameter value, `U` is the future unit normal to `t=constant`.
Within that spatial slice, `n` is the positive unit normal to `y=constant`:
`n=grad_gamma(y)/|grad_gamma(y)|`. This is **not** the full spacetime gradient
of `y` when the time shift is nonzero. The screen is the tangent space to the
coordinate surfaces `t,y=constant`; its metric is `s_AB=g_AB`. Use an oriented
orthonormal screen frame obtained by the positive symmetric inverse square
root of `s` relative to the registered background frame. This fixes the axes
continuously without imposing an additional rotation.

The sphere orientation is the outward ambient orientation,
`epsilon_abc=sum_j EPS[j,a,b,c] n_j`, with `EPS[0,1,2,3]=+1`. This choice fixes
the sign of every displayed magnetic and curl coefficient. Reversing it
reverses all these signs.

The complete scalar metric perturbation in background orthonormal slots is

```
q_yy=A Y,         q_y0=B Y,          q_00=D Y,
q_ya=C v_a,       q_0a=E v_a,
q_ab=L Y kappa_ab + M h_ab,
v_a=D_aY,        h_ab=D_aD_bY.
```

Thus the spatial metric, coordinate covariant shift, and lapse are

```
gamma_yy=1+lambda A Y,       gamma_yA=lambda R C v_A,
s_AB=R^2(kappa_AB+lambda p_AB), p_AB=L Y kappa_AB+M h_AB,
beta_y=lambda rho B Y,      beta_A=lambda rho R E v_A,
N^2=rho^2[1-lambda D Y+lambda^2(B^2Y^2+E^2|v|^2)]+O(lambda^3).
```

In the background orthonormal frame, the first observer variations are

```
U1=(D Y/2) e0 - B Y ey - E v^a ea,
n1=-(A Y/2) ey - C v^a ea,
(ea)1=-(1/2) p_a^b eb.
```

These formulas include the time and radial shifts. Their normalization and
mutual orthogonality hold through first order, including the `q` terms.
The screen epsilon with one lowered and two raised indices has expansion

```
epsilon(s)_A^{CD} = R^-1 [epsilon(kappa)_A^{CD}+lambda e1_A^{CD}]+O(lambda^2),
e1_A^{CD}=p_A^E epsilon_E^{CD}-(tr_kappa p/2)epsilon_A^{CD}.      (4)
```

The second observer/frame variations multiply `C0=0`. First screen-frame
changes applied to the *complete* `H1=0` likewise vanish. Equation (4), however,
cannot be omitted from the quadratic calculation: its contraction with the
unsymmetrized first Codazzi tensor is generally nonzero.

## Why ADM Codazzi gives the full Weyl quantity

Define `H_AB=(1/2) epsilon_(A^{CD} C_{B)UCD}`, with normalized symmetrization.
For three screen slots and one `U`, every Schouten subtraction is a sum of
`s_BC P_UD` and `s_BD P_UC`; the remaining metric contractions with `U` vanish.
The epsilon contraction of each such term is antisymmetric in `A,B`, so its
symmetric part vanishes **at every parameter value**. This includes all Ricci,
scalar-curvature, metric, and observer variations of those subtractions.
Codazzi therefore gives the exact identity

```
H_AB = epsilon_(A^{CD} D^(gamma)_C K_{B)D}.
```

Let `chi_AB=gamma(D_A n,e_B)`, `w_A=K(n,e_A)` and `K^S_AB=K(e_A,e_B)`.
The spatial Gauss formula gives the exact screen split

```
H_AB = epsilon_(A^{CD} [D^(s)_C K^S_{B)D}+chi_{B)C}w_D].       (5)
```

Here `D^(gamma)_C K_BD=D^(s)_C K^S_BD+chi_CB w_D+chi_CD w_B`;
the last term has a symmetric `C,D` pair and vanishes under epsilon. Thus
(5) includes the full observer and radial-screen corrections without a
separate gauge or matter field equation.

The background is conformally flat: after `dz=dy/rho`, it is conformal to
`dz^2+ds_dS4^2`, whose Schouten tensor gives `Riem=g wedge P` directly. Hence
`C0=0`. The following first-order calculation proves the complete defined
`H1=0`; it does not assert that the full five-dimensional magnetic Weyl
three-index tensor vanishes.

## First variations and exact scalar identities

Put `H=theta/rho`, `a=rho_y/rho`. Direct lapse-shift and spatial radial ADM
expansion gives

```
K0^S=H R^2 kappa,        chi0=a R^2 kappa,        w0=0,
K1^S/R^2=alpha Y kappa+beta h,
alpha=(L_t+2theta L)/(2rho)+H D/2-a B,
beta =(M_t+2theta M)/(2rho)-E/R,
w1=eta v,
eta=R(C_t-theta C)/(2rho)-R E_y/2-B/2,
chi1/R^2=(aL+L_y/2-aA/2)Y kappa+(aM+M_y/2-C/R)h.
```

It is useful to subtract tensors whose curl or symmetric epsilon product
vanishes identically at all orders:

```
S=K^S-H s,    J=chi-a s,
S1/R^2=l Y kappa+m h,
l=L_t/(2rho)+H D/2-aB,     m=M_t/(2rho)-E/R,
J1/R^2=(L_y-aA)Y kappa/2+(M_y/2-C/R)h.                        (6)
```

Then (5) is exactly `Htensor=curl_s S+epsilon_(A^{CD}J_{B)C}w_D`.
Because `S0=J0=w0=0`, this form retains the connection variations acting on
`K0` automatically and avoids accidentally dropping a `Gamma2*K0` term.
The first symmetric curl of both `Y kappa` and `h` is zero. The first radial
product is also zero. Therefore the complete `H1` vanishes as a tensor.

For the specific degree-two harmonic, all angular calculations are finite:

```
D_c h_ab=-kappa_ca v_b-kappa_cb v_a-2kappa_ab v_c,
A_ab=epsilon_(a^{cd}v_c h_{b)d},
curl(Y h)=A,       curl(v tensor v)=-A,
epsilon_(a^{cd}h_{b)c}v_d=-A.
```

The first connection variation of `s/R^2=kappa+lambda p` is

```
Gamma1^e_cb=(L/2-M)(delta_b^e v_c+delta_c^e v_b)
            -(L/2)kappa_cb v^e.                              (7)
```

The exact frozen basis calculation, with the orientation above, gives

```
A=2(Tplus-Tminus),  curl Tplus=3Tplus,  curl Tminus=-3Tminus.
```

This is an all-point polynomial identity, not a parity argument or a sampled
fit. The normalized sphere mean norms are `||Tplus||^2=||Tminus||^2=1/6`,
`||A||^2=4/3`; pointwise `|A|^2=8 S0(1-S0)` with `S0=n0^2+n1^2`.

## Quadratic expansion and the indispensable correction terms

In this section `K2` denotes the coefficient of `lambda^2`, **half** its
second derivative. This convention is responsible for the final factor two.
The lapse expansion is

```
N^-1=rho^-1[1+lambda D Y/2
  +lambda^2(3D^2Y^2/8-(B^2Y^2+E^2|v|^2)/2)]+O(lambda^3).
```

Write `z=C/R-aM-M_y/2`. Using
`K_ij=(dot(gamma)_ij-D_i beta_j-D_j beta_i)/(2N)`, its angular second
coefficient has the form

```
K2^S/R^2 = (trace scalar) kappa + f_h Y h + f_v (v tensor v),
f_h=D beta/2+B z,               f_v=(E/R)(L-2M).              (8)
```

For completeness the omitted trace scalar is

```
[H(3D^2/8-B^2/2)+(D/2)(alpha-HD/2)+B(a(A-L)-L_y/2)]Y^2
 +[-H E^2/2+aCE-EL/(2R)]|v|^2.
```

It contributes zero to the symmetric background curl, exactly. Formula (8)
uses the full first spatial connection, including its radial and angular
pieces; no scalar-field or Einstein equation has been substituted.

After dividing by `R^2`, the four contributions to the coefficient of `A` in
(5) at order `lambda^2` are

| Contribution | Coefficient |
|---|---:|
| Curl of `K2^S` | `D beta/2+B z-(E/R)(L-2M)` |
| First screen connection acting on `S1` | `-(L/2-M)m` |
| Epsilon raising/volume correction (4) acting on `D S1` | `M(-l+m)` |
| Radial correction `J1*w1` | `(C/R-M_y/2)eta` |

The sum is

```
F = D beta/2+B z-(E/R)(L-2M)
    +(2M-L/2)m-Ml+(C/R-M_y/2)eta.
```

Substitute (6), `beta=m+HM`, `z=C/R-aM-M_y/2`; the `H D M/2` and `a B M`
pairs cancel exactly, leaving (1). This cancellation establishes that retaining
`A`, the radial shape, lapse, shifts, and screen correction is compatible with
a final result independent of `A,a`; it does not justify dropping them earlier.
Converting coordinate screen components to background orthonormal components
contributes `R^-2`, while (4) contributes `R^-1`. Hence the quadratic second
derivative is `(2F/R)A`, giving (2).

The linear TT Codazzi result is `u_t curl(T)/(2rho^2 c)`, and `u_t=u_x/ell`.
Using the exact eigenvalues gives (3). Scalar-type linear second metric
components have zero complete screen H by the same first-order calculation.
SO(4) equivariance separates vector-derived and TT representations, so a
linear vector sector cannot change this TT coefficient. This statement
isolates the coefficient; it does not supply unknown second fields.

## Reproduction and bounded test strategy

`check_reduction.py` uses exact rational polynomials modulo the unit-sphere
relation. It verifies the third derivative identity, complete TT reconstruction,
scalar curls, both connection variations, their trace and Hessian actions,
all four epsilon correction actions, and the radial shape identity. Normal and
optimized Python must produce the same result. `screen_coefficient.py` is the
small interval-compatible implementation of (1)–(3).

A strong independent check computes ADM curvature at a rational sphere chart
point from arbitrary first metric jets, expands in the parameter, and compares
with (1). Vary all seven amplitudes, all needed first derivatives, `rho,c,theta,a`,
and the chart point. This checks the coordinate/frame conversion independently
of the all-point angular identities. Use separate controls for a diagonal
scalar perturbation (`B=C=E=M=0`, giving `F=0`) and a deliberately nonzero
angular shift; a zero-only test cannot detect lost observer corrections.

Useful exact controls with `rho=c=1`, `theta=0` are: `D=M_t=1` and other
inputs zero gives `F=1/4`, so `4F/R=1`; `E=D=1` gives `F=-1/2`, so `4F/R=-2`;
`C=B=1` gives `F=1/2`, so `4F/R=2`. For a volume-correction mutant, take
`M=L_t=1` with other inputs zero: the exact coefficient is `4F/R=-2`, and
omitting the epsilon correction changes it. Run this as a named rejection
control rather than presuming a passing identity suite detects it implicitly.

For interval propagation use the actual common-gauge amplitude enclosures,
not only comparator values. The first-derivative slots already follow from the
inherited differentiated error certificate; no new fourth metric derivative
is required. Evaluate (1) with outward rational arithmetic, combine it with the
certified `u_x` enclosure, and retain the resulting signed interval even if it
contains zero. Tightness may benefit from propagating correlated master/gauge
jets directly before evaluating the polynomial, but no correlation should be
asserted from separate value boxes. The sign of `u_x` alone cannot establish
the sign of the sum in (3).

The production angular runner authenticates the copied angular source against
an immutable SHA256 before executing those exact bytes. It requires a fresh
nonsymlink report path outside its source directory and any enclosing packaged
research capsule, then creates the output exclusively.
`check_reduction_guards.py` exercises the path and pre-import checks without
rerunning the scientific algebra; write replay reports to external scratch.
