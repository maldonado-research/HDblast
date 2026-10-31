# Registered radial quadratic screen coefficient

Use the retained curvature convention, outward screen orientation, scalar
amplitudes `A,B,C,D,E,L,M`, `R=rho cosh(t)`, `H=tanh(t)/rho`,
`a=rho_y/rho`, and source axial tensor
`Aax_ab=epsilon_(a^{cd} v_c h_{b)d}=2(Tplus-Tminus)`.
All `t,y` derivatives below are ordinary physical-coordinate derivatives.

The complete scalar quadratic second derivative contracted with the
registered outward spatial radial normal is

```
Q_radial[q,q] = (4 F_radial/R) (Tplus-Tminus),

F_radial = ((-A+4M-L) M_y-2M L_y)/4
           + C(A-L)/(2R)
           + (E/R-M_t/(2rho)) eta,

eta = R(C_t-theta C)/(2rho)-R E_y/2-B/2.
```

This is the second parameter derivative itself. The radial linear TT term
is `u_y curl(T)/(2R)`; in particular an equal second TT amplitude contributes
`3 u_y/(2R) (Tplus-Tminus)` without another factor one-half.

## A radial hypersurface calculation

The registered radial normal `n` is the spatial unit normal to `y=constant`
inside `t=constant`. Introduce the unit spacetime normal `n_tilde` to the
timelike hypersurface `y=constant`. Its first variation is

```
(n_tilde)_1 = B Y e0 - (A Y/2) ey - C v^a ea,
n_1         =          -(A Y/2) ey - C v^a ea.
```

Thus their difference is `B Y e0`. The complete first magnetic screen
tensor for `e0` is zero. Since the background Weyl tensor vanishes, this
first normal difference contributes zero to the second screen tensor.
The two prescriptions therefore have the same scalar quadratic second
tensor, although the normals are not identical.

Write `X` for the second fundamental form of the timelike radial
hypersurface with outward unit normal `n_tilde`. On that hypersurface let
`u_tilde` be the future timelike unit normal to its screen, `Ktilde` its
screen second fundamental form, and `w_A=X(u_tilde,e_A)`. Codazzi and the
timelike screen Gauss formula give

```
H_radial,AB = epsilon_(A^{CD}
    [D^s_C X^S_{B)D} - Ktilde_{B)C} w_D].
```

The minus sign is from `g(u_tilde,u_tilde)=-1`. Every Weyl Schouten
subtraction cancels in this symmetric epsilon contraction, as it does for
the time-normal expression. Set `S=X^S-a s` and `J=Ktilde-H s`; the
isotropic terms cancel identically, giving
`H_radial=curl_s S-epsilon J w`, with `S0=J0=w0=0`.

For the radial lapse-shift decomposition the first fields are

```
S1/R^2 = l_r Y kappa + m_r h,
l_r = L_y/2-aA/2+HB,       m_r = M_y/2-C/R,

(J1/R^2)_Hessian = m_t h,  m_t = M_t/(2rho)-E/R,
w1 = eta_r v,
eta_r = R E_y/2-R(C_t-theta C)/(2rho)-B/2 = -eta-B.
```

The un-subtracted Hessian coefficients are
`beta_r=m_r+aM` and `beta_t=m_t+HM`.
Writing `X2` for the coefficient of the squared perturbation parameter
(half its second derivative), its nontrace angular terms are

```
X2^S/R^2 = (trace scalar) kappa + f_h Y h + f_v v tensor v,
f_h = -A beta_r/2+B beta_t,
f_v = (C/R)(L-2M).
```

These follow directly from
`X_ij=(partial_y h_ij-D^h_i beta_j-D^h_j beta_i)/(2Nrad)`.
The radial lapse contributes `Nrad^-1=1-lambda A Y/2+O(lambda^2)`.
The time-shift term has positive `B beta_t`, because the induced radial
hypersurface metric has negative inverse time component. Any omitted
trace coefficient has identically zero symmetric curl for arbitrary
angular scalar coefficient.

Using the source's all-point angular identities, the four coefficients
of `Aax` at quadratic parameter order are

| Contribution | Coefficient |
|---|---:|
| Curl of `X2^S` | `-A beta_r/2+B beta_t-(C/R)(L-2M)` |
| First screen connection acting on `S1` | `-(L/2-M)m_r` |
| Epsilon raising and volume acting on `D S1` | `M(-l_r+m_r)` |
| Timelike screen product `-epsilon J1 w1` | `m_t eta_r` |

Their sum is

```
F_radial = -A beta_r/2+B beta_t-(C/R)(L-2M)
           +(2M-L/2)m_r-M l_r+m_t eta_r.
```

Substitution cancels both `a A M/2` and `H B M` pairs. Collecting terms
then gives the displayed compact formula. Conversion to orthonormal
screen slots contributes `R^-2`, epsilon contributes `R^-1`, and the
second derivative contributes a factor two. Hence the tensor is
`(2F_radial/R) Aax=(4F_radial/R)(Tplus-Tminus)`.

## Direct exact witness and limits

`check_direct_radial.py` is a derived version of the retained independent
coordinate curvature oracle. It constructs the full 5D metric, inverse,
connection, Riemann, Ricci, Schouten and Weyl tensors as exact rational
parameter series at a curved nonnormal stereographic point. All seven
scalar amplitudes and all independent ordinary two-jets are populated.
The registered radial normal is computed from the full spatial inverse
metric, and the physical screen epsilon from its Gram determinant.

The complete result equals the compact formula. Actual mutations omitting
the first radial-normal correction or the first epsilon correction change
the result. With `--combined`, a nonzero equal TT second metric is added
with the required `lambda^2/2` normalization and a nonzero radial jet;
the direct combined tensor tests the radial TT sign and factor.

This finite witness supports the universal algebra above; it is not a
numerical enclosure of the proper-time coefficient. No observer-independent,
detector, radiation, or cosmological-origin conclusion is claimed.
