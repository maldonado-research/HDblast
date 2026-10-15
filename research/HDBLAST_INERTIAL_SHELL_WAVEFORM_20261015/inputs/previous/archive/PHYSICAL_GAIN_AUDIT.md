# Physical shell gain: correlated exact-rational enclosure

15 October 2026 UTC. The same actual regular branch, exact `delta=1/1000`,
`ell=1e-10`, observation endpoint `T1=3ell`, and inward collar `8ell` are
retained. This audit adds no trajectory, field, physical coefficient or
observational scale from another model. No archived program was executed.

**Result.** Using the prior analytic shell box, and preserving the exact
shell relation between its scalar and radius, the physical gain satisfies

\[
\boxed{-\frac{2750293553}{250000000000}\le G\le
-\frac{2493872451}{250000000000}<0},\qquad
G=\frac3{\rho_b^3p_b}.
\]

The proper-time squared-response prefactor satisfies

\[
\boxed{\frac{9328639849756051}{100000000000}\le K\le
\frac{56727115360439913}{500000000000}},\qquad
K=\frac{3\pi^2}{\rho_b^2p_b^2}.
\]

These are exact outward rational brackets. Their decimal values are
`G in [-0.011001174212,-0.009975489804]` and
`K in [93286.398497560510,113454.230720879826]`.
The JSON supplies finer exact endpoints before the outward decimal rounding.
They are conditional consequences of the retained analytic branch theorem;
they do not certify a numerical background trajectory or waveform.

## 1. Exact shell correlation

Write `x=eta_b`, `f(x)=1+c+c x`, with the retained exact
`c=96714321425692/161839258930731`, and

\[
W(x)=1/3+x^2+x^3/3,\quad U=W'^2/2-2W^2/3,
\quad\sigma=2W+\delta f.
\]

The physical scalar junction and Hamiltonian constraint give

\[
p(x)=-\sigma'(x)/2=-W'(x)-\delta c/2
=-2x-x^2-\delta c/2,
\]
\[
H(x):=\rho_b^{-2}=\frac{\sigma(x)^2}{36}
-\frac{\sigma'(x)^2}{48}+\frac{U(x)}6
=\delta A(x)+\delta^2 B(x),
\]
\[
A=Wf/9-W'c/12,\qquad B=f^2/36-c^2/48.
\]

The new standard-library checker independently constructs these polynomials
with rational coefficients, verifies their exact cancellation, and compares
its `A,B` to the retained certificate coefficients. It uses the exact
`eta_b` endpoints in the prior JSON; it does not substitute the coarse outer
`rho in [129,131]`, `p in [-14e-5,-12e-5]` collar rectangle.

Consequently `G(x)=3H(x)^(3/2)/p(x)` and
`K(x)=3pi²H(x)/p(x)²`. Rational interval evaluation over the whole retained
scalar interval verifies `H>0`, `p<0`, `H'<0`, and

\[
\tfrac32H'p-Hp'>0,\qquad H'p-2Hp'>0.
\]

The first numerator makes `G` strictly increasing. The second, together
with `p³<0`, makes `K` strictly decreasing. Therefore their rigorous extrema
are at scalar interval endpoints. This avoids discarding the exact shell
correlation by independently multiplying all radius and scalar ranges.
The receipt also supplies the wider independent-old-box calculation for
comparison.

Positive square roots are bounded by 128-bit dyadic rationals whose squares
are checked exactly. Pi is enclosed using Machin's identity and 40-term
alternating arctangent series, with exact tangent and angle-range checks.
No floating value controls an endpoint or pass/fail decision.

## 2. Optional refinement from the same retained exact constants

The prior producer rounded its proven constants outward to `K_eta<3.196`.
The retained October 11 exact certificate supplies a slightly smaller exact
`K_eta`. Applying that same proven scalar bound, without new analytic
hypotheses or a new trajectory, gives the optional narrower brackets

\[
-\frac{343782179}{31250000000}\le G\le
-\frac{9975608609}{1000000000000},
\]
\[
\frac{11661077540581211}{125000000000}\le K\le
\frac{14181406361697717}{125000000000}.
\]

They appear separately in `correlated_exact_certificate_eta_box`.
`correlated_retained_eta_box` is the result based directly on the already
registered prior scalar interval. Both concern the same actual background.
The modest refinement does not turn the historical floating profile into a
validated trajectory.

## 3. Coupling and physical-output identities

Let `C=rho(3a+2p_y/p)`, `gamma=C_b/2`, and
`alpha=(rho/rho_b)^(3/2)p/p_b`. Since `dy=-rho dz`,

\[
\alpha_z/\alpha=-\rho(3a/2+p_y/p)=-C/2,
\qquad\alpha_b=1.
\]

For the retained coupled reconstruction `V=alpha cosh(t)^(-3/2)w`,

\[
V_z|_b=\cosh(t)^{-3/2}[w_z-\gamma w]_b.
\]

The full metric-scalar reconstruction, with the shell embedding included
before choosing its straight-shell representative, therefore gives

\[
\boxed{S(t)=G\cosh(t)^{-3/2}D(t)},\qquad
D(t)=[w_z-(C_b/2)w](t,0).
\]

The gain is negative. A signed waveform enclosure must preserve that sign
and the half-factor in `C_b/2`. Squaring the response and multiplying the
retained clock/volume weight cancels `cosh³t` exactly:

\[
\mathcal C_{T1}=\frac{\pi^2\rho_b^4}{3}
\int_0^{T1}\cosh^3t\,S(t)^2dt
=K\int_0^{T1}D(t)^2dt.
\]

An optional background dependency reduction uses `h=U_eta/p`, with the exact
identity `h_z=-rho(U_etaeta-h²+4ah)`. It is an auxiliary representation of
the same background. No approximate value of `h` or trajectory is imported
by the gain calculation.

## 4. Inputs, replay and limits

`GAIN_INPUT_SOURCE_PINS.json` records six exact copied references. The
October 11 analytic proof hash is
`f1e2571a551ec122dd7755a13145bae73bd92bd7f7426bdc3002e2f52f8dbfc8`, its exact
certificate hash is
`22c6cd4e505ee2c7e549420355b94a5b3c079ebb6945314f5b285c13bb67835d`, and the
prior analytic shell-box JSON hash is
`8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8`.
The copied prior producer is explicitly a reference and was not executed.

Run only the newly written producer from this directory:

```sh
python -I -B enclose_physical_gain.py
python -I -O -B enclose_physical_gain.py
```

Both modes passed and wrote identical `PHYSICAL_GAIN_ENCLOSURE.json` bytes;
`GAIN_REPLAY_MODES.json` records the hashes. These finite arithmetic checks
do not formally certify the analytic branch proof. The next waveform step
still needs a validated signed enclosure of `D`, its physical-output error,
and any desired time-coordinate conversion. No effective nonlinear epsilon
or remainder bound, observed scale or cosmological-origin inference follows
from this gain audit.
