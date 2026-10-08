# Prospective entire target and all-node stored-state comparison

Prepared 8 October 2026, before any new physical source construction or saved-array
read/decode. This is a finite calculation of a prescribed linear target and exact
represented stored inputs. Cause remains NOT_ESTABLISHED; novelty NOT_ASSESSED.
The full twelve-case pressure/contact certificate and ultraviolet tail remain
UNRESOLVED. No measured state error is reported in this preparation document.

## Fixed problem and coefficient provenance

The target is evaluated at a = -9/2, with U'=W, W'=2ikW-g,
g=4 L²h-2 Lh'-h'', L=-1/s, and zero data at -6. The exact source vanishes
through -5. The sources are h=B(s+4) and h=(s+4)B(s+4), with the same bump,
represented epsilon and represented Pi as the inherited checkpoint. There is
no new vacuum choice at a. The full target is

    W*(k) = - integral[-5,a] exp(2ik(a-s)) g(s) ds,
    U*(k) = - integral[-5,a] Phi_k(a-s) g(s) ds,
    Phi_k(tau) = integral[0,tau] exp(2ikv) dv.

This definition is entire in k, including Phi_0(tau)=tau. Neither construction
nor target evaluation divides by a small k.

The prior independent BUDGET.json contains SHA256 identities for all 22 vectors
of 113 exact real dyadic source coefficients, but does NOT contain those vectors.
Its continuous postprocessor explicitly records this limit. A digest is not a
coefficient vector. Thus this implementation requires 22 prospectively registered
source-model reconstructions, after the new complete registration and external
public GO. Every reconstructed canonical vector must reproduce the corresponding
old SHA256 before it can be used. The inherited independent constructor uses
exact rational recurrences and directed dyadic exponential bounds; its complete
code and dependencies must be pinned by the new registration. Exporting the
reconstructed vectors and their hashes then makes the new result self-contained.
No reconstruction is permitted during fabricated preparation.

For each source the fixed interior consists of eleven parents, starting at
x=s+5=1/128, with right=min(3*left/2,1/2), until right=1/2. Their centers and
halfwidths are c_j,r_j. Each exact polynomial is

    p_j(s) = sum[l=0,112] d_jl (s-c_j)^l.

Use the prior independently certified parent uniform errors

    e_j = max(exported_uniform_error,
              analytic_tail + exported_coefficient_error),
    E = sum_j (2 r_j) e_j.

The new arithmetic must use this ONE authenticated real polynomial family for
all momenta, all four saved source/resolution capsules, and all twelve cases.
Arithmetic output widths at the old nine probes play no role in the new proof.

## Whole-band entire polynomial from normalized exact moments

Set t(s)=2(a-s), so 0 <= t <= 1 on the full prehistory. For the piecewise exact
polynomial family on the interior define exact real moments

    M_n = sum_j integral[parent j] t(s)^n p_j(s) ds.

The exact finite-source responses are

    W_p(k) = - sum[n=0,infinity] i^n k^n M_n/n!,
    U_p(k) = - (1/2) sum[n=0,infinity] i^n k^n M_(n+1)/(n+1)!.

Termwise integration is justified by uniform absolute convergence on every
bounded complex k set. In particular these formulas include k=0 exactly.

Fix N=832, arithmetic precision 512 bits, and 0 <= k <= 256. Let

    J = sum_j 2 r_j sum[l=0,112] |d_jl| r_j^l.

This is an exact rational upper bound on integral |p|; hence |M_n| <= J for all n.
The truncation tails of the two degree-N polynomials satisfy

    T_W = J * 256^(N+1)/(N+1)! / (1 - 256/(N+2)),
    T_U = J * 256^(N+1)/(2 (N+2)!) / (1 - 256/(N+3)).

Indeed the ratio between consecutive omitted positive W terms is at most
256/(N+2)<1; the U ratio is at most 256/(N+3)<1. These are exact rational
geometric majorants, not an empirical convergence test. Numerically
T_W/J < 1.713e-67. This numerical display is descriptive only; the implementation
uses exact rational factorials and ratios.

For efficient arithmetic set v=(s-c_j)/r_j in [-1,1],

    q_j(v) = sum_l (d_jl r_j^l) v^l,
    t(s) = alpha_j - 2 r_j v,  alpha_j = 2(a-c_j).

Start with q_j and repeatedly multiply by the affine polynomial
alpha_j-2r_j v. Integrate each resulting finite polynomial exactly in form,
using outward Arb arithmetic for its coefficients and endpoint evaluations.
This computes all M_0,...,M_(N+1) in approximately eleven quadratic-size native
polynomial loops. Every exact rational input is converted directly through fmpq,
never through decimal or binary64. All source coefficients, including index112,
are retained.

This coordinate choice also controls conditioning: alpha_j+2r_j equals t at
the parent's left endpoint, which is less than one. The sum of absolute
coefficients of q_j(v)*(alpha_j-2r_j v)^n therefore never exceeds that of q_j.
No large raw power (s-c_j)^(-112) is introduced into the moment integrals.
The final oscillatory k evaluation has cancellation, with a coarse e^256
absolute-term amplification. That cancellation is covered by actual 512-bit Arb
balls; it is not assumed harmless. Each entire polynomial is evaluated by the
native ordinary single-point routine at EVERY exact decoded momentum. No
interpolation, fast unstable multipoint evaluation, or nine-probe reuse occurs.

The expected workload is two moment constructions and 49,152 node evaluations
(8,192+16,384 per source). Exact source coefficients can be shared by source;
capsule identities, nodes, saved U/W inputs, quadrature weights and prefix cases
remain distinct. A shared target value at an exactly equal momentum does not
identify the saved states. The implementation presently makes every evaluation.

## Add all target errors before the gate

The omitted cap is the prescribed real source on [-5,-5+delta], delta=1/128.
With H=(3/8)^63, ell=1/(5-delta), D=delta(2-delta),
pd=2(1-delta)/D², rho=0 for positive_B and 1 for signed_uB, the inherited
integration-by-parts proof supplies, for 0<=k<=256,

    C_W(k) = H [pd+rho+2ell+2k
                    +delta(4k²+4k ell+6ell²)],
    C_U(k) = H [(pd+rho+2ell)/2+1
                    +delta(3ell²+2ell+2k)].

Thus the complete exact-target error beyond the arithmetic polynomial balls is

    b_W(k) = T_W + E + C_W(k),
    b_U(k) = T_U + E/2 + C_U(k).

The factors E and E/2 follow from |exp|=1 and |Phi|<=a-s<=1/2 for real k.
A disk error b safely enlarges BOTH real and imaginary component intervals by b.
At k=0 the exact real-source theorem proves both target imaginary components are
zero; only that mathematical target is assigned the singleton zero interval.
No stored quantity is projected or altered.

Export endpoints on the fixed grid 2^-96. Convert each Arb lower/upper endpoint
to its exact rational value, subtract/add b as an exact rational, then floor/ceil
to integer grid endpoints. The complete exported L1 radius is computed FROM
those four final endpoints:

    radius = [(Re_hi-Re_lo)+(Im_hi-Im_lo)] / 2.

The prospective gate is radius <= 10^-18 for EACH exported U and W rectangle.
It includes moment arithmetic, coefficient division, k evaluation, both analytic
tails, source representation, cap, dimensional inflation and grid export. A
failure aborts the registered result; precision, degree, source, gate or node
universe are not changed in response to physical results. This gate concerns
uncertainty in the TARGET, not accuracy of the saved state. It does not replace
the inherited 2e-8 twelve-case full-integral gate.

A source-polynomial J bound can also be checked without evaluating coefficients:
Cauchy bounds give sum_l |d_jl|r_j^l <= 2 M_j + coefficient_error_j, where M_j is
the old disk majorant. The prior JSON makes J<4.85e6 for either source, while
E<4.67e-28. The fabricated stress fixture uses coefficient conditioning larger
than this bound. These observations are preparation evidence; the actual output
radius gate remains mandatory at every node.

## Exact saved state and normalization discrepancy

After the public GO, decode each finite stored binary80 component as its exact
rational value and divide by the SAME exact represented nonzero epsilon:

    U_s = u1/epsilon, W_s = w1/epsilon,
    [delta U] = U_s - [U*], [delta W] = W_s - [W*].

Interval subtraction reverses the target endpoints. The sign convention is
saved minus target throughout. For each strictly positive exact node,

    c_s = Re U_s - Im W_s/(2k)

is an EXACT rational number. Real g and the entire kernels imply the exact
identity Re U*=Im W*/(2k), so c_error=c_s. This is a theorem about the target,
not a claim of stored normalization. The corresponding first-order Wronskian
variation is 2 epsilon c_s. Retain its sign and magnitude without adjustment.

The infrared-stable amplitude coordinate is

    k A_error = -i delta W/2,
    Re(k A_error) = Im(delta W)/2,
    Im(k A_error) = -Re(delta W)/2.

It requires no division by k. For a rectangular enclosure of kA, an exact rational
L1 bound is valid; a rational outward square-root bound for its squared magnitude
is tighter. If |A| is later required at a positive node, divide this bound by the
exact k. The prescribed target itself has a genuine A pole at zero; no bounded-A
continuum premise is inferred. Report all four delta-U/W components as well as c
and kA, since the imaginary constant in deltaU-A can be invisible to the selected
linear stress observables without being a small mode error.

Every capsule receives its own exact node comparisons and each registered prefix
receives its own aggregates, including signed and absolute normalization data.
The finite quadrature sum of a bound is a bound on that finite weighted sum.
These samples supply no regularity or interpolation theorem for the stored
state between nodes. Any continuum envelope, time/momentum quadrature remainder,
full contacts, subsequent numerical evolution or ultraviolet tail needs a
separate argument. Ward closure does not establish incoming-state accuracy.

## API and prospective execution boundaries

`entire_target_engine.py` is a pure library: no filesystem reads, no source
callbacks, no network, no subprocesses. `EntireTarget(rows)` accepts the eleven
ordered dictionaries of exact center, halfwidth, coefficient list and certified
uniform_error. The guarded wrapper authenticates coefficient hashes and the
associated old-budget certificate before constructing it.

`evaluate(k, source)` takes an exact Fraction and one of the two registered source
labels and returns U/W real/imaginary integer endpoint pairs. The common exponent
is -96 and belongs in the output metadata. The guard must pin all constants,
source labels, case/capsule inventory, exact epsilon/Pi, source builder and
certificate references, decoder, exporter and aggregation implementation.
`source=None` is solely the finite-polynomial fabricated test interface; the
actual wrapper must prohibit it and require a registered source label.

The coefficient reconstruction must journal each of exactly22 source constructions
before calling it. Saved-array decode events and compared node counts must be
journaled independently. Shared source coefficients do not reduce the required
four-capsule/49,152-node coverage. An outer bounded worker, fixed single-thread
settings and pinned interpreter/dependency versions remain separate registration
requirements. A suggested conservative new ceiling is 900 seconds/1GiB per
worker, to be selected from fabricated end-to-end evidence before actual entry.

The companion fabricated tests use only explicitly defined finite polynomials.
They check exact constant-source solutions at zero, tiny and large k, preservation
of coefficient112, nonzero exact normalization defects, the sign/factor in kA,
rejection of altered geometry/degree/precision and enforced exported widths. The
full rehearsal covers all49,152 fabricated nodes and all98,304 complex exports.
It authenticates no real capsule and proves no real stored-state accuracy.
