# A residual-to-physical-response theorem for a massive dynamic shell

15 October 2026 UTC. Independent analytic continuation of the registered
HDBLAST coupled master reduction. This note supplies a computable a posteriori
bound for its normal-derivative output. It does not assign numerical values to
an unvalidated background or infer an effective nonlinear amplitude radius.

The retained physical reduction is [PHYSICAL_SHELL_RESPONSE.md](../archive/sources/current/PHYSICAL_SHELL_RESPONSE.md),
equations (9), (11a), (13), and (20). Its identification of the master solution
with the first variation of an exact Einstein–scalar shell family remains a
conditional input. The estimates below concern that master solution itself.

## 1. Model and the exact output normalization

Let

\[
w_{tt}-w_{zz}+W(t,z)w=0,\qquad
-w_z(t,0)+\mu w_{tt}(t,0)+\beta(t)w(t,0)=0,\qquad \mu>0.
\tag{1}
\]

For the registered quadrupole,

\[
W=\frac{35}{4}\operatorname{sech}^2t-\frac{25}{4}
 +\frac12 C_z+\frac14C^2+\frac23\rho^2p^2,
\quad
\mu=(\rho_b B_b)^{-1},
\quad
\beta=\mu\left(\frac{35}{4}\operatorname{sech}^2t-\frac{25}{4}\right)
 +\frac{C_b}{2}.
\tag{2}
\]

Put \(c=C_b/2\), \(h(t)=\cosh^{-3/2}t\), and
\(A=3/(\rho_b^3p_b)\). Since \(\alpha(0)=1\) and
\(\alpha_z(0)=-c\), the exact physical first variation is

\[
S(t)=Ah(t)D(t),\qquad D(t)=w_z(t,0)-c w(t,0).
\tag{3}
\]

The weighted output norm and the quadratic response coefficient therefore obey

\[
\|S\|_{L^2(\cosh^3t\,dt)}=|A|\|D\|_{L^2(0,T)},\qquad
\boxed{\mathcal C_T=\frac{3\pi^2}{\rho_b^2p_b^2}\|D\|_{L^2(0,T)}^2.}
\tag{4}
\]

No estimate of a growing \(\cosh^3t\) weight is needed: the cancellation is
exact. Here the physical clock endpoint is the exact background value
\(\tau_T=\rho_bT\). Comparing outputs at a separately fixed dimensional
endpoint while varying \(\rho_b\) requires an additional endpoint estimate.

## 2. Domain, regularity, and corner compatibility

Fix \(L>T>0\) inside the actual positive-cone collar, and use the shrinking
domain

\[
\mathcal D_{T,L}=\{(t,z):0\le t\le T,\ 0\le z\le L-t\}.
\tag{5}
\]

The shell's backward characteristic triangle \(0\le z\le T-t\) lies inside
this domain with the strict spatial gap \(L-T\). No condition at the
artificial curve \(z=L-t\) is assumed. For computation on a fixed interval
\([0,L]\), a boundary imposed at \(L\) cannot influence \(z=0\) before
time \(L\), provided the numerical trial is a globally smooth reconstruction
on (5), with its residual actually bounded there.

All identities below hold for a classical error \(e\in C^2(\mathcal D_{T,L})\)
and continuous \(W,\beta,F,G\) satisfying

\[
e_{tt}-e_{zz}+We=F,\qquad -e_z+\mu e_{tt}+\beta e=G\quad(z=0).
\tag{6}
\]

One convenient sufficient finite-regularity input envelope is
\(W,F\in C^2\), \(\beta,G\in C^2\), initial displacement
\(e_0\in C^3([0,L])\), initial velocity \(e_1\in C^2([0,L])\), and the
following two corner conditions. Define

\[
a_0=e_0''-W(0,\cdot)e_0+F(0,\cdot),\qquad
a_1=e_1''-W_t(0,\cdot)e_0-W(0,\cdot)e_1+F_t(0,\cdot).
\tag{7}
\]

Then require

\[
-e_0'(0)+\mu a_0(0)+\beta(0)e_0(0)=G(0),
\tag{8a}
\]
\[
-e_1'(0)+\mu a_1(0)+\beta(0)e_1(0)+\beta'(0)e_0(0)=G'(0).
\tag{8b}
\]

The characteristic equations in section 5 and the shell ODE construct the
solution by contraction on short slabs. Their first and second differentiated
equations construct a \(C^3\) solution under this envelope: (8a) matches the
first derivatives of the outgoing characteristic at the corner; (8b) matches
its second derivatives. The undifferentiated outgoing values already match
when \(q(0)=e_0(0)\) and \(v(0)=e_1(0)\). The coefficient bounds are finite
and \(\mu>0\), so finitely many slabs cover (5). Equivalently, extend the
data and coefficients smoothly past \(L\), solve there, and use finite
propagation to remove that extension. This is a sufficient envelope, not a
claim of sharp Sobolev regularity.

For the exact homogeneous solution, set \(F=G=0\) in (7)–(8). The registered
initial bump and velocity vanish on an open neighborhood of the shell, so all
orders of homogeneous compatibility hold there. A \(C^4\) trial and these
exact data produce \(F,G\in C^2\); their error corner conditions then hold
identically by subtraction. They must not be separately set to zero when a
trial has nonzero residual. A grid function by itself is insufficient:
the trial used here must be a reconstruction with the stated regularity.

The master characteristic speed is exactly one in (1). The independent
strict support margin for identifying it with the exact nonlinear family
remains \(L_{\rm corr}+c_*T<Z\), including the entire constraint-correction
annulus; \(L>T\) alone does not replace that retained hypothesis.

## 3. Residual and coefficient errors, with signs fixed

Let \(\widehat w\) be any trial, with approximate coefficients
\(\widehat W,\widehat\mu,\widehat\beta\). Define its displayed residuals by

\[
R=\widehat w_{tt}-\widehat w_{zz}+\widehat W\widehat w,
\qquad
r=-\widehat w_z+\widehat\mu\widehat w_{tt}
       +\widehat\beta\widehat w\quad(z=0).
\tag{9}
\]

For \(e=w-\widehat w\), the exact sources in (6) are

\[
\boxed{F=-R-(W-\widehat W)\widehat w,}
\tag{10a}
\]
\[
\boxed{G=-r-(\mu-\widehat\mu)\widehat w_{tt}(t,0)
          -(\beta-\widehat\beta)\widehat w(t,0).}
\tag{10b}
\]

In particular, a true-coefficient boundary residual \(r\) gives \(G=-r\).
The energy source is \(+Gq'\), whereas the incoming boundary ODE has source
\(+G\). These signs are mutually consistent.

Assume validated enclosures

\[
0<m_-\le\mu\le m_+,\quad |W|\le M(t),\quad
|W-\widehat W|\le\delta_W(t,z),\quad
|\mu-\widehat\mu|\le\delta_\mu,\quad
|\beta-\widehat\beta|\le\delta_\beta(t).
\tag{11}
\]

Every norm below is evaluated over the indicated domain, not only at grid
nodes. Computable nonnegative source envelopes are, for example,

\[
f(t)=\|R(t)\|_{L^2(0,L-t)}
       +\|\delta_W(t)\widehat w(t)\|_{L^2(0,L-t)},
\tag{12a}
\]
\[
g(t)=|r(t)|+\delta_\mu|\widehat w_{tt}(t,0)|
                 +\delta_\beta(t)|\widehat w(t,0)|.
\tag{12b}
\]

The right sides may be replaced by any certified upper enclosures. Interval
enclosures of (2), including \(C_z\), \(\rho\), and \(p\), can supply (11).
No \(W_t\), \(\beta_t\), \(R_t\), or \(r_t\) is needed in the numerical
bound. The regularity envelope and the construction of the trial remain
separate from the quantities appearing in that bound.

## 4. Positive shifted energy and its explicit majorant

Choose arbitrary fixed \(\lambda,\nu>0\). Write
\(q=e(t,0)\), \(v=e_t(t,0)\), and define

\[
E(t)=\frac12\int_0^{L-t}(e_t^2+e_z^2+\lambda e^2)\,dz
       +\frac\mu2v^2+\frac\nu2q^2,\qquad X=\sqrt{2E}.
\tag{13}
\]

The artificial positive weights avoid any assumption that \(W\) or
\(\beta\) is positive. Direct integration by parts gives the exact identity

\[
\begin{split}
E'={}&-\frac12[(e_t-e_z)^2+\lambda e^2]_{z=L-t}\\
 &+\int_0^{L-t}[(\lambda-W)ee_t+Fe_t]dz
 +(\nu-\beta)qv+Gv.
\end{split}
\tag{14}
\]

At \(z=0\), the bulk flux is \(-ve_z\). Substitution of
\(e_z=\mu v'+\beta q-G\) cancels the derivative of \(\mu v^2/2\)
and proves the sign of the last term. At the moving right end the flux is
nonpositive, so it can be discarded without imposing a boundary condition.

Let validated functions \(a,d\ge0\) satisfy

\[
a(t)\ge\max\left\{
\frac{\sup_{0\le z\le L-t}|\lambda-W(t,z)|}{\sqrt\lambda},
\frac{|\nu-\beta(t)|}{\sqrt{m_-\nu}}
\right\},\qquad
d(t)\ge\sqrt{f(t)^2+g(t)^2/m_-}.
\tag{15}
\]

The two potential cross terms are bounded by \(aE\); Cauchy–Schwarz on
the combined bulk and shell velocities bounds the forcing by
\(\sqrt{2E}\,d\). Consequently

\[
X'\le\tfrac12aX+d.
\tag{16}
\]

At zeros of \(X\), use \(\sqrt{2E+\epsilon^2}\) and pass to the limit.
If \(x_0\) encloses

\[
X(0)^2=\int_0^L[e_1^2+(e_0')^2+\lambda e_0^2]dz
                 +\mu e_1(0)^2+\nu e_0(0)^2,
\tag{17}
\]

in the sense \(x_0\ge X(0)\), then

\[
\boxed{X(t)\le B(t):=
 e^{\frac12\int_0^t a}\left[x_0+
 \int_0^t e^{-\frac12\int_0^s a}d(s)ds\right].}
\tag{18}
\]

For initial bounds, using \(m_+\) in the shell velocity term of (17) is
safe. Any explicit supersolution \(B'\ge aB/2+d\), \(B(0)\ge x_0\)
can replace (18), including a rational interval enclosure.

## 5. The derivative observation: a characteristic estimate

The bulk \(H^1\) energy alone does not bound the normal derivative trace.
The required extra identity is elementary and does not differentiate the
source. Put \(J=e_t+e_z\). Then

\[
(\partial_t-\partial_z)J=F-We,
\tag{19}
\]
\[
J(t,0)=e_1(t)+e_0'(t)
       +\int_0^t(F-We)(s,t-s)ds.
\tag{20}
\]

Minkowski's inequality, followed by the substitution \(z=t-s\), yields

\[
\begin{split}
\left\|\int_0^t H(s,t-s)ds\right\|_{L^2_t(0,T)}
&\le\int_0^T\left[\int_s^T|H(s,t-s)|^2dt\right]^{1/2}ds\\
&=\int_0^T\|H(s,\cdot)\|_{L^2(0,T-s)}ds.
\end{split}
\tag{21}
\]

The latter spatial interval lies inside the energy domain. Therefore, if
\(I\ge\|e_1+e_0'\|_{L^2(0,T)}\),

\[
\|J(\cdot,0)\|_2\le
I+\int_0^T\left[f(t)+\frac{M(t)B(t)}{\sqrt\lambda}\right]dt.
\tag{22}
\]

Since \(e_z=J-v\), the desired shell quantity is
\(e_z-cq=J-(v+cq)\). From the positive shell part of (13),

\[
|v+cq|\le\sqrt{1/\mu+c^2/\nu}\,X.
\tag{23}
\]

If \(|c|\le c_*\), we obtain the principal theorem:

\[
\boxed{\|e_z-c e\|_{L^2(0,T;z=0)}\le E_D^{\rm sol}:=
I+\int_0^T\left[f+\frac{MB}{\sqrt\lambda}\right]dt
 +\sqrt{1/m_-+c_*^2/\nu}\left(\int_0^T B^2dt\right)^{1/2}.}
\tag{24}
\]

All quantities on the right can be computed from coefficient enclosures,
initial error enclosures, and trial residual enclosures. There is no unknown
normal derivative or unknown boundary acceleration on that side.

For reference, the boundary characteristic system itself is

\[
\mu v'+v+\beta q=J+G,\qquad q'=v,\qquad
(e_t-e_z)|_{z=0}=2v-J.
\tag{25}
\]

This verifies independently both the incoming sign and the positive boundary
mass. It also underlies the construction invoked in section 2.

## 6. Physical error and a certified response interval

Suppose the trial uses \(\widehat c\), with
\(|c-\widehat c|\le\delta_c\), and let
\(\widehat D=\widehat w_z-\widehat c\widehat w\) at the shell.
Then

\[
\|D-\widehat D\|_2\le
E_D:=E_D^{\rm sol}+\delta_c\|\widehat w(\cdot,0)\|_2.
\tag{26}
\]

For a physical trial \(\widehat S=\widehat A h\widehat D\), if
\(|A|\le A_*\) and \(|A-\widehat A|\le\delta_A\),

\[
\boxed{\|S-\widehat S\|_{L^2(\cosh^3t\,dt)}
\le A_*E_D+\delta_A\|\widehat D\|_2.}
\tag{27}
\]

Alternatively, directly enclose the physical coefficient. Suppose

\[
0<r_-\le\rho_b\le r_+,
\qquad 0<p_-\le|p_b|\le p_+,
\qquad d_-\le\|\widehat D\|_2\le d_+.
\tag{28}
\]

The triangle and reverse-triangle inequalities give the rigorous interval

\[
\boxed{
\frac{3\pi^2}{r_+^2p_+^2}\,[\max(0,d_--E_D)]^2
\ \le\ \mathcal C_T\ \le\
\frac{3\pi^2}{r_-^2p_-^2}\,(d_++E_D)^2.}
\tag{29}
\]

A strictly positive lower endpoint follows when \(d_->E_D\).
Validated quadrature of \(\widehat D^2\), interval square roots, and outward
rounding of (29) suffice; sampled waveforms or ordinary floating-point
quadrature alone do not certify these inputs. A positive interval in (29)
does not itself give a numerical \(\varepsilon\)-radius for the nonlinear
remainder \(\mathcal R_\varepsilon=\varepsilon^2\mathcal C_T+o(\varepsilon^2)\).

## 7. An explicit short-pulse comparator and positive quantitative bound

The theorem can be used before obtaining a numerical waveform. For a
registered \(\ell>0\), define

\[
q_\ell(z)=\begin{cases}
\exp\left(4-\dfrac{\ell^2}{(z-\ell)(2\ell-z)}\right),&\ell<z<2\ell,\\
0,&\text{otherwise},
\end{cases}\qquad h_\ell(z)=e^{-z}q_\ell(z).
\tag{30}
\]

Take the smooth odd extensions \(q_o,h_o\) to the real line. Their smoothness
at zero follows from vanishing on an open shell collar. The free Dirichlet
trial on \(z\ge0\) is

\[
\widehat w(t,z)=\tfrac12[q_o(z+t)+q_o(z-t)]
             +\tfrac12\int_{z-t}^{z+t}h_o(s)ds.
\tag{31}
\]

It has the exact specified initial data, satisfies
\(\widehat w_{tt}-\widehat w_{zz}=0\), and obeys

\[
\widehat w(t,0)=\widehat w_{tt}(t,0)=0,
\quad
\widehat D(t)=\widehat w_z(t,0)=q_\ell'(t)+e^{-t}q_\ell(t).
\tag{32}
\]

Thus \(x_0=I=0\), \(F=-W\widehat w\), and
\(G=\widehat D\). In particular the boundary forcing has the positive sign
in (32), since the trial boundary residual is \(-\widehat D\).

For \(T\ge2\ell\), integration by parts gives the exact identity

\[
\|\widehat D\|_2^2
=\int_\ell^{2\ell}(q_\ell')^2dt
 +\int_\ell^{2\ell}(e^{-t}+e^{-2t})q_\ell^2dt.
\tag{33}
\]

Each half of the bump changes by one over a length \(\ell/2\). Applying
Cauchy–Schwarz on the rising and falling halves separately proves

\[
\boxed{\frac4\ell\le\|\widehat D\|_2^2
\le\ell(16/\ell+1)^2.}
\tag{34}
\]

For the upper bound, set \(u=(z-\ell)/\ell\) and
\(x=[u(1-u)]^{-1}\ge4\). Then
\(\ell|q_\ell'|=|1-2u|x^2e^{4-x}\le16\), because
\(x^2e^{4-x}\) decreases for \(x\ge4\). Also \(0\le q_\ell\le1\).
There is no sign assumption about the falling-half value of \(\widehat D\).

Choose \(T=3\ell\), \(L=8\ell\), and \(\lambda=\nu=1\). Suppose the
actual coefficient box on (5) supplies constants

\[
|W|\le M,\quad |\beta|\le B_\beta,\quad
\mu\ge m>0,\quad |c|\le c_*.
\tag{35}
\]

The explicit safe constants

\[
f_*=M(1+T)\sqrt L,\qquad g_*=16/\ell+1,\qquad
a_* =\max\{1+M,(1+B_\beta)/\sqrt m\},
\quad d_* =\sqrt{f_*^2+g_*^2/m}
\tag{36}
\]

follow from \(|\widehat w|\le1+t\). Equation (18) gives
\(X(t)\le d_*t e^{a_*T/2}\). Hence (24) has the fully explicit enclosure

\[
\boxed{E_D\le E_*:=T f_*
 +\frac12M d_* e^{a_*T/2}T^2
 +\sqrt{1/m+c_*^2}\,d_*e^{a_*T/2}\frac{T^{3/2}}{\sqrt3}.}
\tag{37}
\]

Since the trial shell value is zero, there is no additional \(\delta_c\)
term. If \(E_*<2/\sqrt\ell\), substitute
\(d_-=2/\sqrt\ell\) and \(d_+=\sqrt\ell(16/\ell+1)\) into (29) to get
a strictly positive, computable physical response interval.

For fixed bounds in (35), \(E_*=O(\sqrt\ell)\) as \(\ell\downarrow0\),
while \(d_-\) grows like \(\ell^{-1/2}\); positivity therefore holds for
all sufficiently small \(\ell\). This is a family of explicitly registered
short pulses, not a numerical computation for an earlier fixed pulse width.
Its full constraint completion still needs an enclosing annulus and the
retained strict nonlinear causal margin. The initial higher norms grow as
the pulse narrows, so this result supplies no uniform exact-family amplitude
radius as \(\ell\downarrow0\).

If a conservative rational criterion is preferred, let
\(b=\sqrt{1/m+c_*^2}\) and

\[
K=12\sqrt8M+54\sqrt8M^2+\frac{459M}{2\sqrt m}
          +36\sqrt8Mb+\frac{153b}{\sqrt m}.
\tag{38}
\]

When \(0<\ell\le\min\{1,2/(3a_*)\}\), use
\(e^{a_*T/2}\le e<3\),
\(f_*\le4\sqrt8M\sqrt\ell\), and
\(d_*\le4\sqrt8M\sqrt\ell+17/(\sqrt m\ell)\)
in (37). Each remaining power of \(\ell\) is at least \(1/2\), giving
\(E_*\le K\sqrt\ell\). In particular \(K\ell\le1\) proves
\(\|D\|_2\ge1/\sqrt\ell\), and consequently

\[
\boxed{\mathcal C_{3\ell}\ge\frac{3\pi^2}{r_+^2p_+^2\ell}>0.}
\tag{39}
\]

Equation (39) is conditional on the specified coefficient enclosure and all
retained physical-family hypotheses. It is not an asserted numerical value.

## 8. Executable checks and scope

`check_residual_response.py` checks the source signs, moving-boundary flux,
energy cross-term normalization, characteristic identity, physical weight
cancellation, comparator cross term, and explicit constant algebra. Its
injected wrong-sign and wrong-normalization controls must fail in ordinary
and optimized Python. The script uses explicit exceptions, not `assert`.
The checks do not certify coefficient enclosures, numerical integration, or
the analytic existence and trace arguments. Those roles are stated separately
above so a downstream interval certificate can name its actual inputs.
