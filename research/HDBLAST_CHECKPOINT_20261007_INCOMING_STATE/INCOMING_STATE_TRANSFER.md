# Incoming-state transfer for the unchanged finite-band metric model

Ricardo Maldonado · 7 October 2026 UTC · conditional mathematical continuation

This calculation supplies the missing homogeneous-state counterpart to the
5 October forcing-error kernels. It bounds no actual incoming quantum state.
It evaluates no physical source, saved quantum array or physical trajectory,
changes no registered probe, case, gate or frozen publication file, and does
not complete the numerical pressure/contact certificate.

## 1. Fixed premises and direct operators

Retain conformal time `a=-9/2`, `b=-7/2`, `L=-1/t`, `L'=L²`, the normalized
first-order model `U'=W`, `W'=2ikW-g`, and the same real source and complete
action/subtraction contacts for the two compared states. The physical mass,
minimal coupling, fixed subtraction reference and direct action-derived
operators are the inherited model premises. At `k>0`, differences satisfy

\[
\delta U'=\delta W,\qquad \delta W'=2ik\delta W,
\]
\[
\delta R=\frac{(2k^2+3L^2)\operatorname{Re}\delta U
       -k\operatorname{Im}\delta W-L\operatorname{Re}\delta W}{2k},
\]
\[
\delta P=\frac{(2k^2/3-L^2)\operatorname{Re}\delta U
       -k\operatorname{Im}\delta W-L\operatorname{Re}\delta W}{2k}.
\]

The scaling remains the inherited `a0^4/epsilon` convention. Pressure is
computed from its direct operator before conservation is checked. Identical
contacts cancel in these differences; uncertainty in actual contacts does
not cancel against the true target and remains a separate certificate term.

Use `dmu=k² dk/(2 Pi²)` on fixed `0<k<=K`, with one positive measure constant
throughout. The rational bounds use only `Pi>=3` and cover the inherited
represented constant. They do not replace it by mathematical pi.

## 2. Exact homogeneous state decomposition and normalization

Define from the actual incoming differences

\[
c(k)=\operatorname{Re}\delta U(a,k)
          -\frac{\operatorname{Im}\delta W(a,k)}{2k},\qquad
A(k)=\frac{\delta W(a,k)}{2ik},\qquad E=e^{2ik(t-a)}.
\]

Then `c` is real and constant in time, and

\[
\delta W(t,k)=2ik E A,\qquad
\operatorname{Re}\delta U(t,k)=c+\operatorname{Re}(EA).
\]

The constant imaginary part of `delta U(a)-A` does not enter these linear
stress operators. It is not silently normalized away. In the fixed basis
`v0=exp(-ikt)/sqrt(2k)`, the first-order variation of
`i(v* v'-v v'*)` is `2epsilon c`. Thus exact first-order normalization of
both compared directions implies `c=0`; rounding or an unverified incoming
prescription does not imply that condition. A difference of exact nonlinear
normalized modes requires its own perturbative truncation analysis.

Write `p=Re(EA)` and `q=Im(EA)`. Direct substitution gives

\[
\delta R=\left(k+\frac{3L^2}{2k}\right)c
                   +\frac{3L^2}{2k}p+Lq,\tag{1}
\]
\[
\delta P=\left(\frac{k}{3}-\frac{L^2}{2k}\right)c
             -\left(\frac{2k}{3}+\frac{L^2}{2k}\right)p+Lq.\tag{2}
\]

The coupled phase satisfies `p'=-2kq`, `q'=2kp`, so `p²+q²=|A|²`.
If the homogeneous phase reference is moved from `a` to `a_ref`, the same
solution uses `A_ref=exp[2ik(a_ref-a)] A`; then
`exp[2ik(t-a_ref)] A_ref=EA`. Thus the envelope `|A|` is reference-phase
invariant, while real and imaginary components require the consistent phase.
This is a reexpression of the same state, not a change of incoming data.
The density's leading `k p` terms cancel only because the state phase and
canonical normalization correlation are retained. The pressure's leading
term remains `-2k p/3`. Replacing these correlated quantities by independent
amplitude estimates loses this useful distinction.

## 3. State stress envelopes, including the infrared premise

Suppose independent input analysis proves measurable envelopes
`|c(k)|<=chi(k)`, `|A(k)|<=sigma(k)`. The exact phase-amplitude bounds are

\[
e_R(k,t)\le \left(k+\frac{3L^2}{2k}\right)\chi(k)
 +L\sqrt{1+\frac{9L^2}{4k^2}}\,\sigma(k),\tag{3}
\]
\[
e_P(k,t)\le\left|\frac{k}{3}-\frac{L^2}{2k}\right|\chi(k)
 +\sqrt{\frac{4k^2}{9}+\frac{5L^2}{3}
                         +\frac{L^4}{4k^2}}\,\sigma(k).\tag{4}
\]

Integrate these nonnegative envelopes with the same measure to bound the
continuous-band stresses. Sufficient infrared/finite-band assumptions are
`integral_0^K (k³+k)(chi+sigma) dk < infinity`; the intermediate `k²` term
is bounded by `(k³+k)/2`. Uniform envelopes suffice, but finite samples do
not supply them. If a power-law envelope behaves as `k^-p` near zero, `p<2`
is sufficient for these stress bounds. A represented `k=0` mode itself is
not evaluated; these are improper-integral bounds.

For uniform independently proved `chi,sigma`, elementary positive
coefficient bounds give

\[
e_{R,K}(t)\le C_{R,c}(L,K)\chi+C_{R,A}(L,K)\sigma,
\quad e_{P,K}(t)\le C_{P,c}(L,K)\chi+C_{P,A}(L,K)\sigma,\tag{5}
\]
\[
C_{R,c}=\frac{K^4+3L^2K^2}{8\mathrm{Pi}^2},\quad
C_{P,c}=\frac{K^4}{24\mathrm{Pi}^2}
                              +\frac{L^2K^2}{8\mathrm{Pi}^2},
\]
\[
C_{R,A}=\frac{L}{6\mathrm{Pi}^2}
  \left[\left(K^2+\frac{9L^2}{4}\right)^{3/2}
                         -\left(\frac{9L^2}{4}\right)^{3/2}\right]
 \le\frac{LK^3}{6\mathrm{Pi}^2}+\frac{9L^3K}{16\mathrm{Pi}^2},\tag{6}
\]
\[
C_{P,A}\le\frac{K^4}{12\mathrm{Pi}^2}
                                  +\frac{5L^2K^2}{16\mathrm{Pi}^2}.\tag{7}
\]

Equation (6) integrates the exact density phase norm. Its rational majorant
uses `sqrt(1+x)<=1+x/2`. For (7), the positive pointwise majorant
`2k/3+5L²/(4k)` has square exceeding the exact pressure norm squared by
`21L⁴/(16k²)`. The canonical pressure coefficient uses a triangle inequality.
These integrated coefficients are upper bounds, not claims that the phases
can simultaneously attain the pointwise maxima at every momentum and time.

At all times, use `L<=2/7`, `Pi>=3` to obtain the particularly simple
normalized-state bounds (`c=0` independently established)

\[
e_{R,K}\le\left(\frac{K^3}{189}+\frac{K}{686}\right)\sigma,
\qquad
e_{P,K}\le\left(\frac{K^4}{108}+\frac{5K^2}{1764}\right)\sigma.\tag{8}
\]

This improves the density's generic uncorrelated `K^4` sensitivity to `K^3`
for this homogeneous, normalized state contribution. It leaves the pressure
at `K^4`. It is not a ultraviolet bound or a statement about the actual
state's high-momentum decay.

## 4. Direct work and pressure-time primitives

Differentiating (1), with direct pressure (2), proves

\[
\delta R'=L(\delta R-3\delta P),\qquad
\left(-\frac{\delta R}{3L}\right)'=\delta P.\tag{9}
\]

Provided the envelopes justify absolute domination and the interchange of
time and momentum integration, the exact primitives give

\[
\delta I_K:=\int_a^bL(\delta R_K-3\delta P_K)dt
                   =\delta R_K(b)-\delta R_K(a),\tag{10}
\]
\[
\delta J_{P,K}:=\int_a^b\delta P_Kdt
    =\frac{b\delta R_K(b)-a\delta R_K(a)}3.\tag{11}
\]

These are state-difference primitives derived from the direct operators;
they do not define or repair pressure. They do not include source-work or
contact uncertainties in the full certificate.

With `La=2/9`, `Lb=2/7`, `(b-a)=1`, the uniform envelopes imply

\[
|\delta I_K|\le
\frac{3(L_b^2-L_a^2)K^2}{8\mathrm{Pi}^2}\chi
 +\left[C_{R,A}^{up}(L_b,K)+C_{R,A}^{up}(L_a,K)\right]\sigma.\tag{12}
\]

The constant `k c` term cancels from the work difference exactly. Its
canonical coefficient is `16K²/(1323 Pi²)`, at most `16K²/11907`.
The normalized phase coefficient is at most
`16K³/1701+536K/250047`. A tighter input-dependent result retains the norm
of the complex endpoint coefficient

\[
\left(\frac{3L_b^2}{2k}-iL_b\right)e^{2ik}
                   -\left(\frac{3L_a^2}{2k}-iL_a\right),
\]

before momentum integration; the table uses the valid triangle bound.
For direct time-integrated pressure,

\[
|\delta J_{P,K}|\le
\left[\frac{K^4}{24\mathrm{Pi}^2}
       +\frac{(L_b-L_a)K^2}{8\mathrm{Pi}^2}\right]\chi
 +\frac{C_{R,A}^{up}(L_b,K)/L_b
                  +C_{R,A}^{up}(L_a,K)/L_a}{3}\sigma.\tag{13}
\]

The `sigma` contribution here grows as `K³`, despite the pointwise pressure
envelope's `K⁴`; this is proved oscillatory integration, rather than measured
signed cancellation. It does not bound the time integral of the absolute
pressure error by that smaller quantity.

## 5. What these primitives cannot certify

An arbitrary on-shell homogeneous state perturbation obeys both identities
(9), even if its initial amplitudes are physically incorrect. Also, dropping
the paired canonical geometry terms `3L²c/(2k)` from density and `-L²c/(2k)`
from pressure leaves a separately conserved `(kc,kc/3)` sector. That incorrect
paired omission passes a Ward test. The new verifiers explicitly demonstrate
this blind spot and reject it against the direct operators.

Normalization, conservation and nine momentum probes therefore still do not
establish `sigma(k)` for the actual continuum state. The earlier smooth
between-probe counterexample remains applicable. This theorem identifies
which independently justified spectral enclosures would suffice; it supplies
neither those enclosures nor their ultraviolet continuation.

The twelve-case numerical certificate still needs actual incoming-state
analysis, coefficient/phase arithmetic, reconstructed residuals, direct full
contact/source-work evaluation, exact geometry/product terms as needed,
time/momentum quadrature and endpoint serialization. The analytic source
Taylor remainder is a separate already bounded term and is not counted as
state uncertainty. Source/operator arithmetic bounds at nine probes remain
nine-probe results. A finite-band calculation supplies no ultraviolet tail.

The historical metric calibration remains **FAIL**. The original twelve-case
pressure/contact certificate and its unchanged `2e-8` full-integral gate remain
**UNRESOLVED**. Coupled gravity, heating, cosmological agreement and a
higher-dimensional cause of the Big Bang remain **NOT_ESTABLISHED**. External
mathematical novelty remains **NOT_ASSESSED**. Independent internal checks
are not external peer review or a proof-assistant formalization.

## 6. Next concrete calculation

Audit the upstream incoming-state prescription and the provenance of all
retained modes without changing them. Separate the conditional saved-discrete
target (each stored binary value treated as exact) from a continuum/physical
state target, which requires independent state preparation and spectral
regularity. Construct enclosures for `c(k)` and `A(k)` appropriate to each
target and carry their weighted integrals through (3)-(13). The conditional
discrete target can have exact stored inputs without implying zero error in
the physical continuum target.

Only after the complete implementation, fabricated tests, controls, resource
budget and independent review have been frozen and publicly read back should
new physical callbacks or saved-array evaluations begin. Keep the existing
twelve cases and inherited gates. If state provenance supplies no continuum
enclosure, record that missing prerequisite rather than infer it from Ward
closure or replace the state silently.
