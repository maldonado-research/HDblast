# A resolved first-variation shell waveform with continuous error bounds

15 October 2026 UTC. Conditional analytic and verified-arithmetic checkpoint.

For the same narrow aligned quadrupole as the preceding checkpoint, an
explicit comparison wave now resolves the **linear first variation** of the
physical shell scalar. Uniform pointwise errors, a signed physical gain,
continuous time-cell enclosures and an improved leading-response integral are
supplied. The calculation does not require a PDE time-step solver: it uses
an explicit approximation and a rigorous bound on its difference from the
actual coupled solution. It is not an exact nonlinear waveform or a finite
epsilon prediction. Analytic hypotheses remain mathematical inputs; finite
checks do not formally certify the analytic PDE proof or establish novelty.

## 1. Retained state and the physical quantity

Keep the registered near-phi=+1 regular Einstein–scalar branch, identical
copies and one common shell, with exact parameters

\[
\delta=1/1000,\quad c_{\rm model}=96714321425692/161839258930731,
\quad\ell=10^{-10},\quad T=3\ell,\quad L_{\rm energy}=8\ell.
\tag{1}
\]

The retained proof and background box are copied verbatim under
`inputs/previous/`. The present gain audit additionally uses the exact shell
constraint to preserve parameter correlations. Historical floating trajectories
and other scalar branches supply no validation inputs.

Set Y=n_1^2-n_2^2 on the unit S^3 and use the previous normalized bump

\[
q_\ell(z)=Q(z/\ell),\quad
Q(x)=\begin{cases}
\exp\{4-1/[(x-1)(2-x)]\},&1<x<2,\\0,&\text{otherwise}.
\end{cases}
\]

The full coupled master data are P|_0=alpha q_ell Y and
P_t|_0=e^{-z}alpha q_ell Y, where
alpha=(rho/rho_b)^(3/2)p/p_b. Thus both aligned profiles are nonzero,
A=B=diag(1,-1,0,0). With P=VY and V=alpha cosh^(-3/2)t w, the retained
metric–scalar reconstruction and physical junction laws give

\[
w_{tt}-w_{zz}+\mathcal W(t,z)w=0,\qquad
-w_z+\mu w_{tt}+\beta(t)w=0\quad(z=0),
\tag{2}
\]

\[
w|_0=q_\ell,\quad w_t|_0=e^{-z}q_\ell,\quad
\gamma=C_b/2,\quad\mu=(\rho_bB_b)^{-1},\quad
\beta=\mu(35\operatorname{sech}^2t/4-25/4)+\gamma.
\]

On the actual positive collar, the authenticated bounds are

\[
\mu>m=10^{-3},\quad |\mathcal W|<M=10^6,\quad
|\beta|<B=1000,\quad |\gamma|<\Gamma=250.
\tag{3}
\]

Here B is only a bound on beta, not the background B_b; gamma is half C_b,
not c_model. Define the rescaled time and derivative trace

\[
x=t/\ell=\tau/(\rho_b\ell),\quad
D(t)=w_z(t,0)-\gamma w(t,0),\quad H(x)=\ell D(\ell x).
\tag{4}
\]

For the intrinsic proper-time clock of the previous physical-response
construction, the full scalar first variation is

\[
\widehat\Phi_1(\tau,\Omega)=S(\tau/\rho_b)Y(\Omega),\qquad
\ell S(\ell x)=\mathcal G\cosh^{-3/2}(\ell x)H(x),\quad
\mathcal G=\frac3{\rho_b^3p_b}<0.
\tag{5}
\]

This is the quadrupole coefficient, not an angularly averaged or isotropic
signal. Embedding variation and the harmonic gauge bridge are the retained
physical construction; the constant background shell scalar makes first
clock/label shifts drop out. The relation to the exact family remains
conditional on all prior localized-inverse, reflection, gauge, support and
finite-regularity evolution hypotheses.

Retain a strict margin for the entire nonlinear correction. For example,
conditional on the localized-inverse hypotheses for that annulus, choose
A_0=ell/2, L_corr=3ell, T_evol=7ell/2, c_*=5/4 and Z=8ell. Then
L_corr+c_*T_evol=59ell/8<Z and T<T_evol. The background radial speed is one;
shrinking the exact-family parameter interval supplies c_* without an effective
epsilon value. L_energy and L_corr are distinct quantities.

## 2. Uniform waveform theorem

The free halfline Dirichlet solution with the same initial data is

\[
u(t,z)=\tfrac12[q_o(z+t)+q_o(z-t)]
+\tfrac12\int_{z-t}^{z+t}h_o(s)ds,
\quad h(z)=e^{-z}q_\ell(z),
\tag{6}
\]

where odd extensions are smooth because the data vanish near z=0. It has
u_b=0, |u|<=1+T and u_z(t,0)=q_ell'(t)+e^{-t}q_ell(t). The rescaled
comparison output is therefore the explicit function

\[
\widehat H(x)=Q'(x)+\ell e^{-\ell x}Q(x).
\tag{7}
\]

**Conditional waveform enclosure.** On 0<=x<=3, the actual first-variation
trace satisfies

\[
\boxed{|H(x)-\widehat H(x)|<1.0000004\,10^{-7}.}
\tag{8}
\]

More precisely,

\[
H(x)=0\ (0\le x\le1),\qquad
|H(x)|<4.752\,10^{-14}\ (2\le x\le3),
\tag{9}
\]

and

\[
\|H-\widehat H\|_{L^2(0,3)}<1.0000005\,10^{-7}.
\tag{10}
\]

The post-pulse assertion is a small-tail bound, not a zero-tail assertion.
The exact prefront statement concerns the first variation; higher-order
constraint corrections and changing nonlinear cones prevent it from fixing
the exact nonlinear first-arrival time.

For the signed physical coefficient, (5) reverses the sign of H. Certified
witnesses at x=5/4 and x=7/4 give respectively ell S<0 and ell S>0.
Thus the physical first variation has two opposite-signed lobes. Every
global maximizer of the **normalized trace H** lies in (1.3033,1.3034),
and every global minimizer lies in (1.6966,1.6967), with

\[
4.2356403<\max H<4.2356406,\quad
-4.2356406<\min H<-4.2356403.
\tag{11}
\]

The peak corollary does not assert uniqueness of actual critical points,
zeros or extrema, or determine which lobe attains the absolute maximum.
Detailed high-precision critical-point geometry in the auxiliary note refers
to the comparator unless explicitly labeled otherwise.

## 3. Positive boundary inertia and incoming characteristics

Let A=sup|w| on the causal triangle t+z<=L_energy, 0<=t<=T, and let
V_b=sup|v|, with q_b=w(t,0), v=q_b', q_b(0)=v(0)=0. The free Dirichlet
Duhamel kernel is nonnegative and has total spatial mass at most t-s.
The boundary lift is q_b(t-z) when z<=t and zero otherwise. Consequently

\[
(1-r)A-TV_b\le U,\quad U=1+T,\quad r=MT^2/2.
\tag{12}
\]

The incoming characteristic J=w_t+w_z obeys
(partial_t-partial_z)J=-Ww. Thus at the shell

\[
J(t,0)=q_\ell'(t)+e^{-t}q_\ell(t)+R_J(t),\quad
|R_J(t)|\le MtA.
\]

The complete dynamic boundary equation is

\[
\mu v'+v+\beta q_b=J(t,0).
\tag{13}
\]

Its positive exponential kernel allows the derivative input to be integrated
by parts. The contribution of q_ell' is

\[
\frac{q_\ell(t)}\mu-
\frac1{\mu^2}\int_0^te^{-(t-s)/\mu}q_\ell(s)ds.
\]

Both terms lie in [0,1/mu], so their difference has absolute value at most
1/mu, not a bound proportional to sup|q_ell'|. The other initial input,
the potential correction and beta q_b then give

\[
(1-s)V_b-\frac rm A\le\frac Um,
\quad s=BT^2/(2m).
\tag{14}
\]

The matrix has positive diagonal and positive determinant

\[
d=(1-r)(1-s)-Tr/m>0.
\]

Its positive inverse yields

\[
A\le A_*:=\frac{U(1-s+T/m)}d,\qquad
V_b\le V_*:=\frac U{md}.
\tag{15}
\]

Because D=J-v-gamma q_b and |q_b|<=TV_*,

\[
\eta_*:=\ell[MT A_*+(1+\Gamma T)V_*]
<1.0000004\,10^{-7}
\tag{16}
\]

proves (8). Incoming characteristic control is essential: no pointwise
normal trace is inferred from the preceding L2 energy estimate.

For t>=2ell the leading q_ell(t) term vanishes. Since both initial profiles
are supported on an interval of length ell, the remaining memory is bounded
by ell/m^2+ell/m. Therefore

\[
V_{\rm tail}=\ell/m^2+\ell/m+(r/m)A_*+sV_*,\quad
\eta_{\rm tail}=\ell[MT A_*+V_{\rm tail}+\Gamma TV_*]
<4.752\,10^{-14}.
\tag{17}
\]

Finite propagation gives exact zero before x=1. Integrating the global error
only on [1,2], and the tail error on [2,3], gives
eta_*^2+eta_tail^2<(1.0000005e-7)^2 and proves (10).
All data are smooth and compatible at the shell corner; no auxiliary deep
boundary affects this causal triangle. Independent error-field and direct
solution derivations are included in `referee/`, `astra/` and `sol/`.

## 4. Continuous-cell arithmetic and the physical gain

The correlated retained shell box gives rigorous outward rational bounds

\[
-0.011001174212\le\mathcal G\le-0.009975489804,
\]
\[
93286.398497560510\le\mathcal K:=
\frac{3\pi^2}{\rho_b^2p_b^2}\le113454.230720879826.
\tag{18}
\]

They use p_b=-sigma'_b/2 and the exact constraint
rho_b^(-2)=sigma_b^2/36-sigma_b'^2/48+U_b/6. Derivative signs justify
correlated endpoint extrema. A separately labeled tighter eta-box option
is not needed for this theorem. No rounded historical trajectory is used.

For a real cell [a,b] inside [1,2], put v=1/[(x-1)(2-x)]>=4. Then
Q=exp(4-v), Q'=(3-2x)v^2 exp(4-v), and v^2 exp(4-v) decreases with v.
Exact reciprocal ranges therefore enclose Q and Q' throughout the cell.
At support endpoints the limits are zero, avoiding division by zero or an
unjustified complex-analytic quadrature assumption.

`waveform_intervals.py` uses 96-bit dyadic outward floor/ceil operations.
Exponentials use range reduction to argument<=1/8, a 24-term Taylor sum
with explicit remainder, repeated interval squaring and a separately proved
large-argument tail. Independent positive-Taylor reciprocal checks verify
selected values without importing that implementation.

`compute_waveform.py` exports 4096 closed support cells plus one prefront
cell and one tail cell. Every plotted band covers the whole cell. The
physical factor uses the safe inequality
1-t^2<=cosh^(-3/2)t<=1 for 0<=t<=1/2; together with negative G it is
multiplied outward without dropping its sign. Plotting arithmetic is display
only; the CSV integer endpoints and rational JSON are the certificate data.

The graphic's center curve is the comparator at a midpoint gain, for visual
orientation. The actual first variation is enclosed by the full band. Gain
uncertainty dominates the visible physical-amplitude uncertainty. The tail
inset shows a bound rather than a resolved tail shape.

![Conditional first-variation waveform envelope](figures/SHELL_FIRST_VARIATION_WAVEFORM.png)

## 5. Sharper leading proper-time response coefficient

The exact cancellation of the proper-time cosh weight remains

\[
\mathcal C_{3\ell}=\frac{\mathcal K}{\ell}\int_0^3H(x)^2dx.
\tag{19}
\]

Continuous real interval sums enclose int Hhat^2 on the 4096 support cells.
The support endpoints are included. Outward square roots and the norm
triangle inequality with (10) then give the actual H norm. Combining with
(18) and the Jacobian 1/ell yields

\[
\boxed{6.06\,10^{15}<\mathcal C_{3\ell}<7.44\,10^{15}.}
\tag{20}
\]

Full exact endpoints are in `computed/WAVEFORM_ENCLOSURE.json`; the main
computed interval is approximately [6.069818375e15,7.427495591e15]. This
tightens the preceding [1e15,4e17] enclosure for the **same** state and
normalization. It is not a bound transferred to wider or differently scaled
perturbations. An independent checker reconstructs all 4098 exported cells
and the physical/norm arithmetic without importing the production evaluator.

## 6. Scientific limits and next target

The exact-family proper-time response still has
R_epsilon=epsilon^2 C+o(epsilon^2), conditionally on the retained hypotheses,
and its leading coefficient is independent of the second initial correction.
There is no effective nonlinear parameter/remainder bound. No finite-amplitude
nonlinear waveform, physical time calibration, exact nonlinear arrival time,
gravitational energy, cosmological fit or higher-dimensional Big Bang origin
is established. Internal independent AI derivations/review and finite exact
checks are distinct from external peer review and formal proof verification.

The next step is a wider, practically chosen state with a validated background
tube and continuous coupled characteristic/shell reconstruction, followed by
an effective nonlinear remainder estimate. The present small-width analytic
waveform is a verified benchmark for that implementation. No new Zenodo
version or APS submission accompanies this checkpoint.
