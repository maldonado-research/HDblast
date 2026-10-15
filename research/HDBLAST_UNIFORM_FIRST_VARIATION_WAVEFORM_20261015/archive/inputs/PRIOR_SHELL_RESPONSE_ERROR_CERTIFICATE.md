# A quantitative leading shell-response enclosure from a free-wave comparator

15 October 2026 UTC. Conditional analytic checkpoint, with exact-rational
arithmetic verification. This note extends the preceding physical-response
checkpoint. It supplies a coarse numerical interval for one newly specified
leading response coefficient. It supplies no numerical waveform, effective
nonlinear amplitude radius, observed time scale or cosmological-origin test.
The analytic inputs remain assumptions/theorems to be assessed mathematically;
the finite programs do not formally certify those inputs or the PDE proof.

## 1. The state, result and retained hypotheses

Use the same five-dimensional Einstein–scalar model, identical interior
copies, one common shell, scalar potential and shell tension as the retained
near-phi=+1 regular branch. The exact parameters are

\[
\delta=1/1000,\qquad c_{\mathrm{model}}=
96714321425692/161839258930731.
\]

The analytic branch and its boundary box are authenticated in
`archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md`. Historical floating
trajectories and the different S85/negative-phi branch are not inputs.

Retain the earlier localized exact-data and finite-time evolution hypotheses:
reflection symmetry, the specified adjoint/weighted inverse, support control,
the completed initial gauge and uniqueness, and the strict causal margin for
the **whole** nonlinear constraint correction. The background collar used
below does not replace those hypotheses. The evolution interval must extend
strictly beyond the observation endpoint. No effective interval in the
nonlinear parameter epsilon is supplied.

Select

\[
\ell=10^{-10},\qquad T_1=3\ell,\qquad L=8\ell.
\tag{1}
\]

Here z and t are the dimensionless coordinates of the preceding checkpoint.
These numbers specify a very narrow mathematical perturbation, not an
inferred physical length or cosmological time. Let

\[
q(z)=\begin{cases}
\exp\{4-\ell^2/[(z-\ell)(2\ell-z)]\},&\ell<z<2\ell,\\
0,&\text{otherwise},
\end{cases}\qquad Y=n_1^2-n_2^2.
\]

The symbol L in this note is the **energy triangle's initial length**,
not the endpoint L_corr of the nonlinear correction annulus. A sufficient
causal geometry, conditional on the retained localization/inverse hypotheses
for that annulus, is A_0=ell/2, L_corr=3ell, T_evol=7ell/2,
c_*=5/4 and Z=8ell: then L_corr+c_*T_evol=59ell/8<Z and T_1<T_evol.
The background radial speed in (t,z) is one; a bound c_*=5/4 for the exact
nearby metrics follows after shrinking their parameter interval. This is
an existence-smallness statement with no effective epsilon value. Other
admissible geometries satisfying the retained strict margins work as well.

As previously, with \(\alpha=(\rho/\rho_b)^{3/2}p/p_b\), prescribe the
**coupled** master data

\[
P|_0=\alpha qY,\qquad P_t|_0=e^{-z}\alpha qY.
\tag{2}
\]

Both displacement and velocity are nonzero; A=B=diag(1,-1,0,0). Full
metric–scalar reconstruction and the previous harmonic gauge bridge identify
this master with the first variation of the conditional exact family.
This is a member with the explicit new value (1), not a numerical evaluation
of an earlier unspecified profile width.

**Conditional quantitative result.** For the previous intrinsic-proper-time
physical scalar response functional,

\[
\mathcal R_\varepsilon(\rho_b T_1)
=\varepsilon^2\mathcal C_{T_1}+o(\varepsilon^2),\qquad
\boxed{10^{15}<\mathcal C_{T_1}<4\,10^{17}.}
\tag{3}
\]

The interval uses the normalization (2) and the fixed model units. It is
deliberately coarse. Large coefficients for narrow normalized data give no
finite nonlinear prediction without an effective epsilon/remainder bound.
The leading coefficient remains independent of the second initial correction
among exact curves with the same first tangent. This result does not compute
the waveform, gravitational energy or the exact nonlinear first-arrival time.

## 2. Physical output and the exact cancellation of the clock weight

Set \(P=VY\) and \(V=\alpha\cosh^{-3/2}t\,w\). On the actual positive
collar the preceding full reconstruction gives

\[
w_{tt}-w_{zz}+\mathcal Ww=0,\qquad
-w_z+\mu w_{tt}+\beta(t)w=0\quad(z=0),
\tag{4}
\]

\[
\mu=(\rho_bB_b)^{-1}>0,\quad
\gamma=C_b/2,\quad
\beta=\mu(35\operatorname{sech}^2t/4-25/4)+\gamma.
\]

The symbol gamma is **half** the published coefficient C_b. It is unrelated
to c_model above. The physical shell scalar, including embedding variation
before taking the straight-shell representative, is

\[
S(t)=\frac{3\cosh^{-3/2}t}{\rho_b^3p_b}
[w_z(t,0)-\gamma w(t,0)].
\]

Consequently the cosh weight of the proper-time observable cancels exactly:

\[
\mathcal C_{T_1}
=\frac{\pi^2\rho_b^4}{3}\int_0^{T_1}\cosh^3t\,S^2dt
=\underbrace{\frac{3\pi^2}{\rho_b^2p_b^2}}_{K}
\|D\|_{L^2(0,T_1)}^2,\quad D=w_z-\gamma w|_{z=0}.
\tag{5}
\]

Thus an error certificate must control the derivative trace D, rather than
only the bulk displacement. The clock, gauge and exact-family statements
used in (5) are the retained physical-response result, not consequences of
the arithmetic checks here.

## 3. An energy estimate that controls the derivative observation

For an error e on the shrinking triangle \(0\le t\le T<L\),
\(0\le z\le L-t\), suppose

\[
e_{tt}-e_{zz}+\mathcal We=F,\qquad
-e_z+\mu e_{tt}+\beta e=G\quad(z=0).
\tag{6}
\]

Assume a smooth solution, or sufficient energy/trace regularity obtained by
smooth approximation with compatible corner data. Coefficients are bounded,
mu is constant and positive, and residual norms below are integrable. The
present comparator and retained solution are smooth, so no weak trace
extension is needed for (3). Fix lambda,nu>0 and define

\[
X^2=\int_0^{L-t}(e_t^2+e_z^2+\lambda e^2)dz
+\mu v^2+\nu q_b^2,\quad q_b=e(t,0),\quad v=q_b'.
\]

Integration by parts and the boundary law give

\[
\tfrac12(X^2)'=
\int[(\lambda-\mathcal W)ee_t+Fe_t]dz
+(\nu-\beta)q_bv+Gv
-\tfrac12[(e_t-e_z)^2+\lambda e^2]_{z=L-t}.
\tag{7}
\]

The positive boundary inertia is essential. The moving right edge has
nonpositive flux: no artificial boundary condition is imposed there.
Let mu>=m>0, |W|<=M, |beta|<=b, and uniformly ||F(t)||_2<=f, |G(t)|<=g.
Then a safe constant and forcing bound are

\[
a=\max\{(\lambda+M)/\sqrt\lambda,
(\nu+b)/\sqrt{m\nu}\},\qquad h=\sqrt{f^2+g^2/m}.
\]

From (7), \(X'\le aX/2+h\), interpreted by a regularized square root
at X=0. In particular

\[
X(t)\le B_*:=e^{aT/2}(X(0)+Th).
\tag{8}
\]

This by itself does not bound e_z at the shell. Instead put J=e_t+e_z.
Its incoming characteristic law is
\((\partial_t-\partial_z)J=F-\mathcal We\). All its characteristics
ending on the shell lie in the shrinking triangle since L>T. Minkowski's
integral inequality therefore gives

\[
\|J(\cdot,0)\|_2\le I+T(f+MB_*/\sqrt\lambda),\quad
I=\|e_t(0,\cdot)+e_z(0,\cdot)\|_{L^2(0,T)}.
\]

Since e_z=J-v and
\(|v+\gamma q_b|\le\sqrt{1/m+\gamma_*^2/\nu}\,X\),

\[
\boxed{\|e_z-\gamma e_b\|_2\le
I+T(f+MB_*/\sqrt\lambda)
+\sqrt{T(1/m+\gamma_*^2/\nu)}\,B_*=:\Delta.}
\tag{9}
\]

Here gamma_*>=|gamma|. There is no differentiation of a residual or time
derivative of W or beta in this estimate. General coefficient uncertainty
is included by computing the true residuals of an approximate solution u:

\[
F=-r-(\mathcal W-\widehat{\mathcal W})u,\quad
G=-r_b-(\mu-\widehat\mu)u_{tt,b}
-(\beta-\widehat\beta)u_b.
\]

If an approximate output uses gamma_hat, add
\(|\gamma-\widehat\gamma|\|u_b\|_2\) to (9). Neither coefficient errors
nor boundary acceleration may be silently omitted in that general case.
For (3), u_b=0 and the coefficients used in (6) are the true ones.

## 4. An explicit comparison wave and elementary residual bounds

Let q_o and v_{0,o} be the smooth odd extensions of q and v_0(z)=e^{-z}q(z).
Define the free Dirichlet solution

\[
u(t,z)=\tfrac12[q_o(z+t)+q_o(z-t)]
+\tfrac12\int_{z-t}^{z+t}v_{0,o}(s)ds.
\tag{10}
\]

It has exactly the initial values in (2) after conjugation, u_b=0, and

\[
u_z(t,0)=q'(t)+e^{-t}q(t)=:\widehat D(t),\quad
|u|\le1+T.
\tag{11}
\]

The error e=w-u has X(0)=I=0. Its true residual sources are
F=-Wu and G=D_hat. No boundary residual is discarded: the free Dirichlet
wave is only a comparator, not an asserted solution of the shell law.

The normalized bump obeys 0<=q<=1 and |q'|<=16/ell. To see the latter,
write x=(z-ell)/ell and v=1/[x(1-x)]>=4. Then
|ell q'|<=v^2 exp(4-v)<=16. It follows that

\[
f=M(1+T)\sqrt L,\qquad g=16/\ell+1.
\tag{12}
\]

On the rising half, q increases from 0 to 1 over length ell/2 and v_0>=0.
Cauchy–Schwarz on that interval and (11) give the conservative bounds

\[
2/\ell\le\|\widehat D\|_2^2
\le\ell(16/\ell+1)^2\qquad(T=3\ell).
\tag{13}
\]

The difference of normal traces is exactly the output in (9), since u_b=0.
Thus

\[
\max\{0,\sqrt{2/\ell}-\Delta\}\le\|D\|_2
\le\sqrt\ell(16/\ell+1)+\Delta.
\tag{14}
\]

These are actual-solution estimates. They do not use the contradictory
zero-extension lower bound of the preceding qualitative activation proof.

## 5. Rigorous background box and exact arithmetic

The separately proved first-exit collar bound, obtained from the retained
analytic branch theorem, supplies on 0<=z<=8e-10, for all real t:

| Quantity | Sufficient rigorous bound |
| --- | --- |
| mu | 1/1000 < mu < 1/400 |
| abs(W) | < 10^6 |
| abs(beta) | < 1000 |
| abs(gamma) | < 250 |
| rho_b squared | 16641 < rho_b^2 < 17161 |
| p_b squared | 9/625000000 < p_b^2 < 49/2500000000 |

Its inward ODE is evaluated in a rational outer box with strict first-exit
margins. The potential uses the exact derivative-free background identity

\[
\mathcal W=-\tfrac{35}4\tanh^2t+
\rho^2[2h_{\rm bg}^2-10ah_{\rm bg}+\tfrac{35}4a^2-U_{\eta\eta}-p^2/6],
\quad h_{\rm bg}=U_\eta/p.
\]

The root certificate takes lambda=nu=1. Outward exact-rational arithmetic
in (8)–(14) gives Delta<8. With the positive prefactor interval from (5),
it proves the strict coarse enclosure (3). `RESPONSE_CERTIFICATE_NORMAL.json`
retains the full rational endpoints, not rounded display values.

`certify_response.py` uses integer-square-root dyadic enclosures and the
rigorous rational bound exp(x)<=(1-x/n)^(-n), n>=max(1,2x). Pi is enclosed
by Machin's identity and alternating arctangent remainders. Floating-point
inputs are rejected. These checks certify arithmetic conditional on the
stated analytic inequalities. They neither enclose a numerical trajectory
nor constitute formal verification of the analytic PDE argument.

Independent Astra energy/trace derivation, Sol characteristic majorants,
critical review, source audit and finite wrong-formula controls are included
separately. Their scope and counts are recorded by the fresh replay receipt.

## 6. Remaining scientific target

The next useful calculation is a resolved shell waveform on a practically
chosen scale with a substantially sharper observation-norm enclosure.
One must propagate validated coefficient/background errors and all bulk,
boundary, initial and reconstruction residuals. Resolution convergence alone
is insufficient. To predict a finite exact nonlinear response, also obtain
an effective epsilon/remainder bound and connect model normalization to a
physical scale. No external peer review, novelty certification, global
stability, universal-profile observability or evidence for a higher-dimensional
origin of the Big Bang is asserted here.
