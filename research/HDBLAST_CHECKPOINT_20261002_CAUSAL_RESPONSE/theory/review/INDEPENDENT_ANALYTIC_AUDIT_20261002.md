# Independent analytic audit of the prospective causal-response protocol

Scope: analytic review only. No pulse/source samples, mode solutions, response
integrals, quadrature experiments, or new physical predictions were evaluated.
No checkout file was changed. The reviewed proposal is
`/workspace/HDblast/research/HDBLAST_CHECKPOINT_20261002_STATIONARY_QUANTUM/theory/PROSPECTIVE_CAUSAL_RESPONSE_PROTOCOL.md`.

## Response identities

The Kubo sign, Wick multiplicity, momentum measure, and forced-mode response
agree. In particular, multiplying the commutator coefficient
`-i/(2 pi^2)` by the Kubo coefficient `-i/2` gives `-1/(4 pi^2)`.
The physical-field interaction produces exactly the canonical source
`s(eta')=a(eta')^2 delta x(eta')` and the output factor `a(eta)^(-2)`.

At fixed geometry and fixed positive reference r, varying the stated
order-two subtraction gives `delta |v|_sub^2=-s/(4 omega_r^3)`.
The hard-cutoff finite part therefore retains both `gamma_E` and `+1`.
The logarithmic-memory expression follows by integration by parts with
`s(eta_i)=0`. The past-infinite stationary physical-mass source is a separate
state/switching limit, correctly yielding
`Q_x=-(2 gamma_E+ln(r/H^2))/(16 pi^2)` at the plane-wave reference.

The initial Bogoliubov term is correct for
`v_k=[exp(-ik(eta-eta_i))+beta(k)exp(+ik(eta-eta_i))]/sqrt(2k)`:
`delta Q_initial=integral dk k Re[beta exp(2ik(eta-eta_i))]/(2 pi^2 a^2)`.
For an unsqueezed occupation-diagonal homogeneous Gaussian state, the
commutator response acquires the stated `1+2n(k)` factor. A general squeezed
Gaussian state would require its anomalous correlations as well; the stated
occupation diagnostic has the required interpretation.

The current chain rule contains both
`x_phi delta Q/2` and `x_phiphi Q0 delta phi/2`. The displayed mass-law contact
is correct. The linear energy Ward identity contains `Q0 delta x_dot/2`;
the product `delta Q delta x_dot` begins at second order.

Three integrations by parts give
`|I_k-s(eta)/(2k)| <= (|s''(eta)|+integral |s'''|)/(8k^3)`
when the initial source derivatives vanish. Together with the subtraction
bound this gives exactly the proposal's combined omitted-tail coefficient.
These estimates do not replace quadrature or floating-point error bounds.

## Analytic derivative norms for the normalized bump

Let `B(u)=exp(1-q)`, `q=1/(1-u^2)`, for `|u|<1`, and let B and every
derivative vanish outside this interval. Set `F(u)=u B(u)`.
All following statements are global, including the smooth endpoints.

On `0<u<1`, define

    P(q) = 4q^4 - 12q^3 + 6q^2,
    R(q) = 4q^4 - 12q^3 + 2q^2,
    C(q) = 2q^3 - 14q^2 + 21q - 6,
    S(q) = 4q^4 - 32q^3 + 64q^2 - 36q + 3.

Then

    B''(u)  = exp(1-q) P(q),
    B'''(u) = -4u q^3 exp(1-q) C(q),
    F''(u)  = sqrt(1-1/q) exp(1-q) R(q),
    F'''(u) = -2q^2 exp(1-q) S(q).

These are derivatives with respect to u. The normalization here is B(0)=1;
the example bump `exp[-1/(1-u^2)]` in the proposal differs by a constant e.

### Root counts and exact total variation

Under `z=q-1`,

    C = 2z^3 - 8z^2 - z + 3,
    S = 4z^4 - 16z^3 - 8z^2 + 12z + 3.

Each polynomial has exactly two positive z roots: Descartes' rule gives
at most two, and the following disjoint sign-changing brackets supply two.
Write the q roots as alpha_B,beta_B and alpha_F,beta_F, respectively:

    1 < alpha_B < 2,       5 < beta_B < 81/16,
    7/4 < alpha_F < 15/8,  21/4 < beta_F < 85/16.

The exact signs are

    C(1)=3, C(2)=-4, C(5)=-1, C(81/16)=2049/2048,
    S(7/4)=129/64, S(15/8)=-1023/1024,
    S(21/4)=-879/64, S(85/16)=101937/16384.

Let `b(q)=exp(1-q)P(q)` and
`f(q)=sqrt(1-1/q)exp(1-q)R(q)`.
The signs fix every monotonicity interval; there are no additional extrema.
Since `B''(0)=-2`, `F''(0)=0`, and both vanish at u=1,

    integral_{-1}^1 |B'''(u)| du
      = 4[b(beta_B)-b(alpha_B)] - 4,

    integral_{-1}^1 |F'''(u)| du
      = 4[f(beta_F)-f(alpha_F)].

These exact algebraic-root/exponential formulas need no response evaluation.

### Convenient conservative certificates

For the negative B'' extremum,

    -b(alpha_B) <= 12/e < 9/2,

because `6q-2q^2-3 <= 3/2` and `q^2 exp(1-q) <= 4/e`.
For its positive extremum, P is increasing throughout its root bracket, so

    b(beta_B) <= P(81/16) exp(-4)
              = 20056977/(16384 e^4) < 45/2.

The last strict inequality follows from the elementary Taylor lower bound
`e > 163/60`. Thus

    integral |B'''| < 104,       sup |B''| < 45/2 < 23.

For the negative F'' extremum, use that `q^2 exp(1-q)` increases below q=2,
`6q-2q^2-1` decreases above q=3/2, and `sqrt(1-1/q)` increases. The bracket
then gives

    -f(alpha_F) <= (6075/256) sqrt(7/15) exp(-7/8) < 34/5.

For the positive extremum, R and the square root increase on its bracket:

    f(beta_F) <= (23647425/16384) sqrt(69/85) exp(-17/4) < 187/10.

Entirely rational certificates for these two last inequalities are

    sqrt(7/15) < 137/200,
    exp(-7/8) < 167/400,
    sqrt(69/85) < 901/1000,
    exp(-17/4) < (7/19)^4 (384/493).

The exponential bounds follow from Taylor lower sums:
`exp(7/8) > 9429967/3932160 > 400/167`,
`e > 163/60 > 19/7`, and `exp(1/4) > 493/384`.
Consequently

    integral |F'''| < 4(34/5+187/10) = 102,
    sup |F''| < 187/10 < 19.

No simple interval-length-times-supremum bound on the third derivatives is
used: the certificates exploit the exact total variation of the second
derivatives and the complete root count.

## Scaling and direct use in the registered tail bound

For `u=(eta-eta_c)/sigma`, with sigma>0, and either
`s(eta)=epsilon B(u)` or `s(eta)=epsilon F(u)`,

    integral |s'''(eta)| d eta
      <= (|epsilon|/sigma^2) times 104 or 102,

respectively. A global observation-independent bound is therefore

    A_eta <= 127 |epsilon|/sigma^2   for B,
    A_eta <= 121 |epsilon|/sigma^2   for uB.

After the pulse, s''(eta)=0, so 104 and 102 directly replace 127 and 121.
Before the pulse, both s and A_eta are exactly zero. During the pulse one
can retain the actual analytic |s''(eta)| and use the full derivative-norm
certificate, or a separately justified truncated norm.

With dimensionless a, eta and sigma have length dimension; k and M have
mass dimension; epsilon and s have mass-squared dimension. Thus A_eta has
mass dimension four and the stated tail bound has the required mass-squared
dimension of Q. If epsilon is recorded as dimensionless, its chosen mass
scale, such as H^2, must be written explicitly. These envelopes support
freezing the numerical inputs; they establish no numerical accuracy result.

## Read-only review of the staged interval implementation

The sibling `tail_bounds.py` was read, not imported or executed. Its exact
partial-A formulas are correct for the frozen unit half-width pulse centered
at eta=-4. Write `D=maximum-minimum` for the appropriate second-derivative
extrema. At u=-1/2 and u=0, A/|epsilon| equals `2D` for both sources. At
u=+1/2 it equals

    2D - 4 + (832/81) exp(-1/3)   for B,
    2D     + (992/81) exp(-1/3)   for uB.

After the pulse it equals `4D-4` and `4D`, respectively. These expressions
follow from the monotonicity intervals certified above and the exact values

    B''(+/-1/2)=-(416/81) exp(-1/3),  B''(0)=-2,
    F''(+1/2)=-(496/81) exp(-1/3),
    F''(-1/2)=+(496/81) exp(-1/3),    F''(0)=0.

The implementation's wider signed-maximum root bracket (5,6) is also valid:
`S(5)=-77` and `S(6)=363`. The exact root count ensures uniqueness inside
each selected bracket. Exact rational bisection followed by outward-rounded
rational endpoints and interval primitive evaluation encloses the required
extrema. Its tail expressions retain interval arithmetic for amplitude,
scale factor, M squared, pi, and cutoff. Final conversion advances each
upper binary64 endpoint toward positive infinity. This is a static analytic
review of the enclosure construction, not a completed runtime check of its
mpmath dependency or of any numerical experiment.
