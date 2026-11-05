# A wider shell waveform with a certified signed inertial tail

15 October 2026 UTC. Conditional analytic and verified-arithmetic checkpoint.

For the newly specified width ell=10^-6, the actual coupled **linear first
variation** now has a resolved positive normalized post-pulse tail:

\[
\boxed{3.2\,10^{-8}<H(x)<5.3\,10^{-8}\qquad(2\le x\le3).}
\tag{1}
\]

The physical quadrupole coefficient has the opposite, negative sign.
Continuous sign witnesses also establish at least one zero in each of
(1.49998,1.49999) and (1.96,1.97). A fresh response-integral calculation gives
6.06e11<C_(3ell)<7.44e11. These are results for the wider initial tangent,
obtained by retaining shell inertia in an explicit comparison wave and
bounding its difference from the actual coupled solution. They do not
describe an exact finite-amplitude nonlinear waveform. The retained
analytic and exact-family hypotheses remain mathematical inputs.

## 1. The wider state and the physical observation

Keep the registered near-phi=+1 regular Einstein–scalar background branch,
identical copies and one common shell, at the exact parameters

\[
\delta=1/1000,\qquad c_{\rm model}=96714321425692/161839258930731.
\]

The perturbation width is now ell=1/1000000, with T=3ell and
Z=L_energy=8ell. The previous completed waveform used ell=10^-10. Its
smaller absolute comparison error and response coefficient belong to that
different initial tangent; neither is transferred here. The present
calculation keeps the same background and rechecks its collar at the wider
depth by a strict interval first-exit argument.

Let Y=n_1^2-n_2^2 on the unit S^3, and define

\[
Q(x)=\begin{cases}
\exp\{4-1/[(x-1)(2-x)]\},&1<x<2,\\
0,&\text{otherwise},
\end{cases}\qquad q_\ell(z)=Q(z/\ell).
\]

The retained master data are P|_0=alpha q_ell Y and
P_t|_0=exp(-z)alpha q_ell Y, with
alpha=(rho/rho_b)^(3/2)p/p_b. Both aligned profiles are nonzero:
A=B=diag(1,-1,0,0). The retained reduction
P=alpha cosh^(-3/2)t wY gives

\[
w_{tt}-w_{zz}+\mathcal W(t,z)w=0,\qquad
-w_z+\mu w_{tt}+\beta(t)w=0\quad(z=0),
\tag{2}
\]
\[
w(0,z)=q_\ell(z),\qquad w_t(0,z)=e^{-z}q_\ell(z),
\]
\[
\mu=(\rho_bB_b)^{-1},\quad\gamma=C_b/2,\quad
\beta(t)=\gamma+\mu(35\operatorname{sech}^2t/4-25/4).
\]

Here B_b is a background coefficient; the later constant B bounds beta.
The derivative trace and its physical reconstruction are

\[
D(t)=w_z(t,0)-\gamma w(t,0),\quad H(x)=\ell D(\ell x),
\quad x=\frac t\ell=\frac\tau{\rho_b\ell},
\tag{3}
\]
\[
\widehat\Phi_1(\tau,\Omega)=S(\tau/\rho_b)Y(\Omega),\qquad
\ell S(\ell x)=\mathcal G\cosh^{-3/2}(\ell x)H(x),
\quad\mathcal G=\frac3{\rho_b^3p_b}<0.
\tag{4}
\]

The intrinsic proper-time clock, embedding variation and harmonic gauge
bridge are the retained construction. The constant background shell scalar
removes first-order clock/label shifts from this observable. A sign assertion
for S concerns its quadrupole coefficient; Y itself changes sign with angle.
This is not an angular average, isotropic signal or gravitational energy.

The same compact outer background rectangle
rho in [129,131], eta in [-9/100000,-8/100000],
p in [-14/100000,-12/100000], a in [111/1000,112/1000]
has positive first-exit margins through Z=8e-6. Independent differentiation
and exact interval evaluation give the safe constants

\[
\mu>m=1/500,\quad |\mathcal W|<M=102000,\quad
|\beta|<B=250,\quad |\gamma|<\Gamma=250,
\quad |\mathcal W_z|<K=3\,10^7.
\tag{5}
\]

The sharper retained parameter intervals are approximately
0.0020597712<mu<0.0022716618 and
202.08605<gamma<247.36781. Decimal displays do not decide any gate.
The derivative audit also proves
|W_t|<21/400000 and |beta'|<21/160000000 on [0,T].
See [the background audit](archive/BACKGROUND_DERIVATIVE_AUDIT.md).

For the conditional nonlinear-family interpretation, the newly specified
tangent must separately satisfy the retained localized-inverse, reflection,
gauge, whole-correction support and finite-regularity evolution hypotheses.
A strictly separated example is A_0=ell/2, L_corr=3ell,
T_evol=7ell/2 and c_*=5/4, since
L_corr+c_*T_evol=59ell/8<Z and T<T_evol. Shrinking an unspecified family
parameter interval may supply the speed bound; this is not an effective
epsilon radius. L_corr and L_energy have different roles.

## 2. Exact free comparator with boundary inertia

Set beta0=gamma+5mu/2 and a0=35/4. Then
beta(t)=beta0-a0 mu tanh^2t. Let u solve the free bulk equation with the
same initial data and constant dynamic boundary
-u_z+mu u_tt+beta0 u=0. Its boundary value b(t)=u(t,0) satisfies exactly

\[
\mu b''+b'+\beta_0b=j_0(t),\qquad
j_0=q_\ell'+e^{-t}q_\ell,\qquad b(0)=b'(0)=0.
\tag{6}
\]

The term b' comes from the incoming characteristic. With impulse kernel k,

\[
\mu k''+k'+\beta_0k=0,\qquad k(0)=0,\quad k'(0)=1/\mu,
\]

one has b=k*j0. Integration by parts of the known q_ell' input gives

\[
b(t)=\int_0^t[k'(t-s)+e^{-s}k(t-s)]q_\ell(s)\,ds,
\]
\[
b'(t)=q_\ell(t)/\mu+
\int_0^t[k''(t-s)+e^{-s}k'(t-s)]q_\ell(s)\,ds.
\]

Put a=1/mu, c=beta0/mu and

\[
P(r,s)=(a-\gamma-e^{-s})k'(r)+(c-\gamma e^{-s})k(r).
\]

The exact comparison trace is

\[
H_0(x)=Q'(x)+\ell(e^{-\ell x}-1/\mu)Q(x)
+\ell^2\int_1^{\min(x,2)}P(\ell(x-y),\ell y)Q(y)\,dy,
\tag{7}
\]

where the integral is zero for x<=1. This expression retains the prescribed
initial velocity in both its local term and memory kernel. For fixed x>=2,
its small-width expansion is

\[
H_0(x)=\frac{\ell^2I_Q}{\mu}
\left(\frac1\mu-\gamma-1\right)+O(\ell^3),
\qquad I_Q=\int_1^2Q(y)\,dy.
\tag{8}
\]

The minus one in (8) is the initial-velocity contribution. The sign theorem
does not rely on this formal expansion: the complete finite-width kernel
and its arithmetic remainder are enclosed below. The comparator solves
its stated free linear problem; the difference from the actual potential
and variable boundary coefficient remains in the error proof.

## 3. Actual-minus-comparator error and cancellation

Let w_b=w(t,0) and
I(t)=integral_0^t W(s,t-s)w(s,t-s)ds. Incoming transport gives
w_t+w_z=j0-I at the boundary. For the boundary error q=w_b-b, v=q',

\[
\mu v'+v+\beta_0q=-I+a_0\mu\tanh^2t\,w_b,
\qquad q(0)=v(0)=0,
\tag{9}
\]
\[
H-H_0=\ell(-I-v-\gamma q).
\tag{10}
\]

Use the positive matrix majorants

\[
U_0=1+T,\ r=MT^2/2,\ s=BT^2/(2m),\
d=(1-r)(1-s)-Tr/m>0,
\]
\[
A_*=U_0(1-s+T/m)/d,\qquad V_*=U_0/(md),
\tag{11}
\]

which bound |w| and |w_b'|. They follow from the free Dirichlet Duhamel
representation and the positive boundary exponential kernel; the derivative
input is integrated by parts before estimation.

Write W(t,z)=V(z)-a0 tanh^2t. On the incoming ray the free Dirichlet field
has displacement part [q_ell(t)+q_odd(t-2s)]/2. Its unweighted integral is
t q_ell(t)/2, exactly zero once t>=2ell. Pairing positive and negative odd
arguments before estimating the weighted integral leaves a coefficient
difference bounded by (K+a0 t)r. The velocity part is nonnegative on this
ray, with integral (1/2)integral_0^t r exp(-r)q_ell(r)dr. Thus, with

\[
R_0=\frac38(K+a_0T)\ell^2+\frac34M\ell^2,
\]

the exact boundary lift and remaining bulk Duhamel term give

\[
|I(t)|\le \frac{Mtq_\ell(t)}2+R_0+
\frac{MV_*t^2}{4}+\frac{M^2A_*t^3}{6}.
\tag{12}
\]

The leading term disappears after the pulse, but its previous history still
forces the error boundary equation. Retaining that history yields

\[
I_{\rm tail}=R_0+MV_*T^2/4+M^2A_*T^3/6,
\]
\[
J_1=3M\ell^2/4+TR_0+MV_*T^3/12+M^2A_*T^4/24,
\]
\[
E_v=\frac{J_1/m+a_0V_*T^4/4}{1-s},\qquad
|v|\le E_v,\quad |q|\le TE_v.
\tag{13}
\]

Consequently

\[
\sup_{0\le x\le3}|H-H_0|<1.1\,10^{-7},\qquad
\sup_{2\le x\le3}|H-H_0|<1.7\,10^{-10}.
\tag{14}
\]

Exact evaluation gives stronger errors below 1.022e-7 and 1.645e-10.
The first variation and comparator both vanish for 0<=x<=1. The support
aware norm bound is ||H-H0||_(L2(0,3))<1.100002e-7.
No pointwise normal trace is inferred from an L2 energy estimate.

For a closed time cell [l,r] with Q<=Q_max throughout, the exported error
uses the minimum of the global bound and

\[
\ell\left[\frac{M\ell rQ_{\max}}2+R_0+
\frac{MV_*(\ell r)^2}{4}+\frac{M^2A_*(\ell r)^3}{6}
+(1+\Gamma T)E_v\right].
\tag{15}
\]

This refinement is needed to resolve the late zero witnesses where Q is
small. The proof is detailed in the
[independent residual derivation](referee/INDEPENDENT_INERTIAL_BOUND.md).
An independent energy-kernel route in the
[Astra derivation](astra/INERTIAL_COMPARATOR_AND_POSITIVE_TAIL.md)
obtains a looser tail error below 2.97e-10 and independently proves H>9.5e-9.

## 4. Continuous arithmetic, signed tail and zero existence

The comparator evaluator uses 112-bit outward integer intervals. A
degree-eight matrix Taylor enclosure for (1000k,k') has the explicit
k' remainder bound 2a_upper(1000r_upper)^9/9!. The factor two follows
from 1000T=3/1000 and exp(1000T)<2. Separate short-time inequalities prove
k'>0 and k>=0 on the complete needed interval. Exponentials include their
remainders; support endpoint limits of Q and Q' are handled explicitly.

The bump primitive is enclosed by 4096 continuous interval rectangles.
Its monotonicity and the positive kernel enclose the integral in (7) over
whole time cells. The export contains one exact prefront cell, 4096 support
cells and 256 post-pulse cells: 4353 closed cells covering [0,3]. Integer
endpoints in the [CSV](computed/INERTIAL_WAVEFORM_CELLS.csv) certify full
cells, not sampled time points. The physical factor is enclosed by
1-t^2<=cosh^(-3/2)t<=1 and multiplied outward with the negative gain.

The full-interval comparator tail is approximately
[3.23575026180407e-8,5.27014430788258e-8]. Enlarging it by the rigorous
1.7e-10 error proves (1) for every real x in [2,3]. This establishes a
strictly nonzero signed tail of the actual coupled first variation.

The actual-H signs at x=5/4 and 7/4 are positive and negative, respectively,
so the physical coefficient retains opposite-signed lobes. Further
strict witnesses have signs

| x | Sign of actual H | Sign of physical coefficient S |
| --- | --- | --- |
| 1.49998 | positive | negative |
| 1.49999 | negative | positive |
| 1.96 | negative | positive |
| 1.97 | positive | negative |

Continuity proves at least one zero in each stated open bracket. It does
not prove uniqueness within a bracket, a count of all zeros, monotonicity,
or uniqueness of extrema. The zero interval [0,1] is a separate exact
prefront statement for this linear variation. Higher-order corrections and
changing nonlinear characteristic cones prevent it from fixing an exact
nonlinear arrival time.

## 5. Newly computed proper-time response coefficient

The correlated retained background gain and norm prefactor have outward bounds

\[
-0.011001174212\le\mathcal G\le-0.009975489804,
\]
\[
93286.398497560510\le\mathcal K:=\frac{3\pi^2}{\rho_b^2p_b^2}
\le113454.230720879826.
\]

These background-only quantities apply to the unchanged shell. In contrast,
the waveform integral is recomputed for the new width. The exact proper-time
cosh cancellation and the t=ell x Jacobian give

\[
\mathcal C_{3\ell}=\frac{\mathcal K}{\ell}\int_0^3H(x)^2\,dx.
\tag{16}
\]

In particular, the factor is 1/ell, not ell. Outward squares of the actual
whole-cell intervals, integrated over all 4353 cells including the tail,
give approximately

\[
6.506603243394882\le\int_0^3H^2\,dx\le6.546734713874441.
\]

The exact rational endpoints in the
[computed report](computed/INERTIAL_WAVEFORM_ENCLOSURE.json) give
approximately [6.069775830288547e11,7.427547506963038e11] for C_(3ell), hence

\[
\boxed{6.06\,10^{11}<\mathcal C_{3\ell}<7.44\,10^{11}.}
\tag{17}
\]

This is a new integral certificate, not a width rescaling of the old
coefficient. Conditional on the exact-family hypotheses,
R_epsilon=epsilon^2 C+o(epsilon^2); the leading coefficient is independent
of the second initial correction. There is no effective epsilon or nonlinear
remainder estimate in this checkpoint.

## 6. Residual-energy method and literature limits

The accompanying [public methods note](millennium/HDBLAST_METHODS_TRANSFER_PUBLIC_SAFE.md)
derives a complementary route for a future approximate coupled waveform.
For an error e with bulk residual r_e and boundary residual q_e, its scalar
energy includes the boundary kinetic term mu|e_t(0)|^2/2. With no far-end
flux the exact energy identity contains both coefficient terms
(1/2)integral W_t|e|^2 and beta'|e(0)|^2/2, in addition to bulk and boundary
residual work. Initial error and both residuals must be controlled.

Under the audited bounds W>=10000, beta>=200 and the stated derivative
bounds, that note derives

\[
A_e(t)\le(1+10^{-14})\left[A_e(0)+
\int_0^t(\|r_e\|_2+23|q_e|)\,ds\right],\qquad
|e(t,0)|\le A_e(t)/\sqrt{300},
\]

where A_e=sqrt(2E_e). This controls a boundary value under its regularity,
coercivity and support premises. It does not by itself control the normal
derivative observable in (3), prove nonlinear evolution, or close the
missing nonlinear defect estimate. The present derivative-waveform proof
uses the incoming-characteristic estimate (10) instead.

The [bounded primary-source review](literature/INERTIAL_WAVEFORM_PRIMARY_LITERATURE.md)
discusses causal Green kernels, characteristic-path expansions and nonlinear
Einstein boundary theories. It records exactly which passages were read,
which external hypotheses do not match this shell model, and which search
failed. Neither those methods nor this finite review establishes novelty.
No external result supplies the unresolved model-specific inverse,
support, gauge or nonlinear-evolution hypotheses.

## 7. Scope and reproduction

The [README](README.md) links the production calculation, independent proofs,
arithmetic checkers and provenance. Inputs are authenticated; normal and
optimized Python modes use explicit checks that survive assertion removal.
Meaningful fault controls include boundary/kernel sign changes, omitted
local inertia, a wrong bump derivative and an incorrect response Jacobian.
Finite exact checks supplement written analytic arguments. Internal
independent AI derivation and review are not external peer review or formal
verification of the PDE theorem.

The result is a conditional first-variation waveform on one specified
background and initial tangent. It supplies no finite-amplitude nonlinear
waveform, effective nonlinear parameter range, observational time calibration,
cosmological fit, gravitational energy or higher-dimensional Big Bang origin.
The earlier full-second-order-contact checkpoint already constructs a
conditional complete second Taylor coefficient through shell contact; its
existence statement is preserved. The next quantitative target is an explicit
validated second-order observable for this wider state and a nonlinear
defect-to-error map compatible with the actual constraints, junction laws,
gauge and support.
