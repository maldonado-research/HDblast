# Uniform pointwise enclosure for the retained narrow shell state

15 October 2026 UTC. Independent Astra analytic derivation. This note concerns
the actual-background first variation with delta=1/1000 and ell=10^-10. It
retains the frozen source's smooth analytic background, exact-data and
finite-time evolution hypotheses, gauge/reconstruction identification, and
strict causal margin for the whole nonlinear constraint correction. It does
not supply a finite nonlinear amplitude radius or physical-units prediction.

## 1. Inputs, output, and the pointwise result

On the characteristic triangle 0<=t<=T=3ell, 0<=z<=L-t, L=8ell, let

    w_tt - w_zz + W(t,z) w = 0,
    -w_z(t,0) + mu w_tt(t,0) + beta(t) w(t,0) = 0.

The actual positive-inertia collar in the frozen response certificate gives

    mu >= m=1/1000, |W| <= M=10^6,
    |beta| <= b=1000, |gamma| <= g=250, gamma=C_b/2.

All weak inequalities here are safe outer bounds for that source's strict
inequalities. Only these actual coefficient bounds are used in this note;
the coefficient functions are not replaced with constant toy coefficients.
The true coefficients need not be differentiated for the bound below.

The initial data are w(0,z)=q_ell(z), w_t(0,z)=exp(-z)q_ell(z), where

    q_ell(z) = exp(4 - ell^2/((z-ell)(2ell-z))) for ell<z<2ell,
             = 0 otherwise.

Let u be the free Dirichlet wave formed from the smooth odd extensions of
these initial data. Then u(t,0)=0, |u|<=U:=1+T, and

    Dhat(t) = u_z(t,0) = q_ell'(t) + exp(-t) q_ell(t).

The actual observed derivative is D=w_z-gamma*w at z=0. Define x=t/ell,

    Q(x) = q_ell(ell*x),
    H(x) = ell*D(ell*x),
    Hhat(x) = Q'(x) + ell*exp(-ell*x)*Q(x).

The following uniform enclosure is established below:

    sup_{0<=x<=3} |H(x)-Hhat(x)| <= epsilon_star
      = 8000003000800720779892143999973
        / 79999999999992799998920000162000000000
      < 1000001/10^13.

Thus epsilon_star is about 1.00000037510018e-7. This is a pointwise
derivative-output estimate for the actual retained solution, rather than an
inference of pointwise control from an L2 observation estimate.

## 2. Incoming characteristic and the boundary equation

Set e=w-u, q_b(t)=e(t,0), v=q_b', and J=e_t+e_z. Initial e and e_t vanish.
The true error equations are

    e_tt-e_zz = -W(u+e) =: F,
    -e_z + mu*v' + beta*q_b = Dhat.

In particular, the boundary forcing is +Dhat: the free Dirichlet comparator
has boundary residual -Dhat. At the shell e_z=J-v, and therefore

    (partial_t-partial_z) J = F,
    mu*v' + v + beta*q_b = J(t,0) + Dhat(t),
    q_b'=v.

All incoming characteristics lie within the triangle and reach the initial
edge. Consequently

    J(t,z) = integral_0^t F(s,z+t-s) ds.

Put E=sup_triangle |e| and V=sup_[0,T] |v|. These are finite for the retained
smooth solution, and |F|<=K:=M(U+E). Hence |J(t,z)|<=tK.

For K_-=e_t-e_z, an outgoing characteristic reaches either the initial edge
or the physical shell. In the first case |K_-|<=tK. In the second case its
shell foot is (t-z,0), where K_-=2v-J. The incoming source interval has
length t-z and the outgoing interval length z, so

    |K_-(t,z)| <= 2V+tK.

It follows that |e_t(t,z)|<=V+tK. The vertical segment from (0,z) to (t,z)
lies in the triangle, so integration gives

    |e(t,z)| <= tV + t^2 K/2,
    E <= TV + k(U+E),     k=MT^2/2.                         (1)

The reflected characteristic source intervals are disjoint in time. No
artificial boundary condition at the moving right edge has been introduced.

## 3. The cancellation that makes the pointwise bound sharp

Positive mu gives the damped boundary variation-of-constants formula

    mu*v(t) = integral_0^t exp(-(t-s)/mu)
                 [q_ell'(s)+h(s)+J(s,0)-beta(s)q_b(s)] ds,
    h(s)=exp(-s)q_ell(s).

Integration by parts, using q_ell(0)=0, transforms the large derivative:

    integral_0^t exp(-(t-s)/mu)q_ell'(s) ds
      = q_ell(t) - mu^-1 integral_0^t exp(-(t-s)/mu)q_ell(s) ds.       (2)

Both q_ell(t) and the weighted integral on the right belong to [0,1], since
0<=q_ell<=1 and mu^-1 integral_0^t exp(-(t-s)/mu)ds=1-exp(-t/mu)<=1.
Their difference therefore has absolute value at most 1. Also

    integral_0^t exp(-(t-s)/mu)h(s) ds <= ell,

because h is nonnegative, at most 1, and supported in an interval of length
ell. The remaining terms satisfy

    integral_0^t exp(-(t-s)/mu)|J(s,0)| ds <= T^2 K/2,
    integral_0^t exp(-(t-s)/mu)|beta(s)q_b(s)| ds <= b T^2 V/2,

where |q_b(s)|<=sV and exp(-(t-s)/mu)<=1. Dividing by mu>=m yields

    V <= A + c(U+E) + dV,                                  (3)
    A=(1+ell)/m, c=MT^2/(2m), d=bT^2/(2m).

Equation (2) is essential to the sharp number here. It is not the invalid
claim that a large pointwise derivative is small. Its damped integral is
bounded by two quantities between zero and one.

## 4. Closed exact-rational majorant and the observation

Suppose k<1 and

    delta_den=(1-k)(1-d)-cT > 0.

Solving the two scalar inequalities (1) and (3) gives

    V <= V_star := [A(1-k)+cU]/delta_den,
    E <= E_star := [T V_star+kU]/(1-k).

This elimination is valid directly for the finite suprema of the solution;
it does not assume their differentiability or a numerical fixed point.

At the shell, u_b=0, so the exact observed difference is

    D-Dhat = e_z-gamma*q_b = J-v-gamma*q_b.

Using |q_b(t)|<=tV_star gives the uniform result

    |D-Dhat| <= MT(U+E_star)+(1+gT)V_star,
    epsilon_star=ell*[MT(U+E_star)+(1+gT)V_star].             (4)

For the stated rational inputs,

    k=d=9/200000000000000,
    c=9/200000000000,
    delta_den=39999999999996399999460000081/40000000000000000000000000000,
    V_star=40000000004000000000000360000000/39999999999996399999460000081,
    E_star=120000018012005405399191619999757
           /399999999999963999994600000810000000000.

These are exact values. Decimal approximations are V_star=1000.00000010009
and E_star=3.000000450300405e-7. The apparently large unscaled boundary
velocity is consistent with the scaled waveform error of order ell/m.

The file `check_uniform_pointwise.py` checks these rational identities,
positive denominators and claimed bound, the characteristic/boundary signs,
and the integration-by-parts identity. Deliberately wrong signs, output
normalization, and a discarded boundary forcing must be rejected. Its
finite checks audit algebra and arithmetic conditional on the analytic
proof; they do not formalize existence or the underlying physical reduction.

## 5. What is resolved, and the physical normalization

The entire first-variation graph is enclosed by

    Hhat(x)-epsilon_star <= H(x) <= Hhat(x)+epsilon_star,
    0<=x<=3.

Any certified interval evaluation of Hhat over an interval gives a
continuous-domain enclosure of H there by adding this error. This supports
rigorous signs, peak heights, and peak-location exclusion where the margin
exceeds the error. Values from a plot or a grid alone do not supply those
claims. Uniform function error alone does not identify an exact number of
zeros or local extrema of H; an exact waveform can have small extra
oscillations inside the enclosure.

The physical scalar first variation is

    S(t)=3*cosh(t)^(-3/2)*D(t)/(rho_b^3*p_b).

The archived actual branch has p_b<0, so S has the opposite sign from D.
The signed normalized physical waveform

    H_phys(x)=ell*rho_b^3*p_b*S(ell*x)/3
             =cosh(ell*x)^(-3/2)*H(x)

has the same orientation as H. It obeys, for example,

    |H_phys-Hhat| <= epsilon_star +(3/4)T^2(16+ell),

using |Hhat|<=16+ell and 0<=1-cosh(t)^(-3/2)<=3t^2/4. The latter inequality
follows from log(cosh t)<=t^2/2 and exp(-a)>=1-a for a>=0. The additional
clock-conjugation term is about 1.08e-18; it is not a physical scale
calibration. Actual unnormalized S intervals also require the retained
rho_b and p_b intervals.

Finite propagation and the identically zero shell-neighborhood initial data
give H=0 on 0<=x<=1. For x>=2 the comparator is zero, and q_ell(t)=0
removes the leading term in (2). Its convolution term can now be bounded
using its support length rather than the whole positive kernel:

    V_tail = ell/m^2 + ell/m + c(U+E_star) + d V_star,
    |H(x)| <= epsilon_tail
       := ell*[MT(U+E_star)+V_tail+gT V_star]
        = 190040072038996407198648650081
          /3999999999999639999946000008100000000000000
        < 4752/10^17,             2<=x<=3.

The global V_star remains in the beta-memory and gamma*q_b terms; the
post-pulse velocity bound cannot be substituted retroactively for the
earlier history. No exactly vanishing late tail is claimed. Combining the
zero pre-arrival interval, global pulse bound, and this tail bound gives

    ||H-Hhat||_L2(0,3)^2 <= epsilon_star^2+epsilon_tail^2.

The exact first nonlinear arrival and nonlinear waveform remain outside
this first-variation certificate.

## 6. Certified signs and resolved peak geometry

The independent exact-rational verifier in `peak_analysis/verify_peaks.py`
evaluates the explicit comparator and its extrema. Its analytic derivation
is `peak_analysis/peak_analysis.txt`. Put

    a=sqrt(2/sqrt(3)-1),
    x_L=(3-a)/2, x_R=(3+a)/2,
    M_Q=6(2+sqrt(3))*a*exp(1-sqrt(3)).

Then Q' has exactly the two interior extrema +M_Q at x_L and -M_Q at
x_R. Exact-rational evaluation gives

    4.235640422134888218565987509535 < M_Q
      < 4.235640422134888218565987509536.

The small nonzero velocity input is included in Hhat, whose positive and
negative peaks are certified separately in that verifier. The resulting
actual waveform conclusions, using epsilon_star<1000001/10^13, are

    H(5/4)>0, H(7/4)<0,

    every global maximizer of H lies in (1.3033,1.3034),
    every global minimizer of H lies in (1.6966,1.6967),

    4.235640322182880388580255 < max H
      < 4.235640522183080388580257,
   -4.235640522086896048554158 < min H
      < -4.235640322086696048554156.

For the location claim, the exact intervals for M_Q-Q'(1.3033) and
M_Q-Q'(1.3034) lie respectively above 2.635596782936e-7 and
5.903033349921e-7. Both exceed 2*(1000001/10^13)+ell=2.000003e-7.
Monotonicity of Q' away from its unique positive peak excludes every point
outside the stated left interval; reflection gives the right interval.
This is an exclusion based on function-value gaps, so it requires no
derivative bound on the actual waveform error.

There is at least one actual zero between 5/4 and 7/4. The absolute
maximum of |H| lies in the union of the two peak windows. The uniform
error is too large to decide which lobe attains it, and does not establish
unique actual extrema or a unique actual zero. The raw physical S signs
are reversed: S(5ell/4)<0 and S(7ell/4)>0 because p_b<0.

The peak verifier uses only exact Fraction intervals, integer-square-root
bounds, and alternating Taylor enclosures for exponentials. Normal and
optimized runs give identical JSON, and four deliberate wrong formulas
are rejected. Its high precision comparator constants do not imply the
same accuracy for the actual waveform or its physical gain.

## 7. Source scope

Read against the immutable source capsule
`validated-response-20261015/frozen-source/HDBLAST_SHELL_RESPONSE_ERROR_CERTIFICATE_20261015`:

* `SHELL_RESPONSE_ERROR_CERTIFICATE.md`, sections 1, 2, 4 and 5: state,
  physical output, comparator, actual coefficient box.
* `archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md`: actual analytic
  branch and rational first-exit collar bounds, including p_b<0.
* `astra/RESIDUAL_TO_PHYSICAL_RESPONSE.md`: retained characteristic signs,
  regularity/compatibility, reconstruction and exact-family scope.

This is a new pointwise continuation on the same narrow state. It does not
relabel the prior energy/L2 estimate as a waveform computation, and it does
not modify the frozen source.
