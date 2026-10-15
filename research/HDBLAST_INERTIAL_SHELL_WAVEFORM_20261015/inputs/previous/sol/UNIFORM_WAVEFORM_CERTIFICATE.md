# A uniform characteristic enclosure of the specified shell waveform

15 October 2026 UTC. Independent Sol derivation for the retained coupled
scalar master. All bounds are conditional on the retained background,
reconstruction, exact-data and evolution hypotheses. The finite checks verify
the rational deductions and bookkeeping; they are not a formal PDE proof.

The specified narrow state admits the following **uniform waveform enclosure**:

\[
\boxed{\quad |H(x)-\widehat H(x)|\le\eta
 <\frac{10000004}{10^{14}}=1.0000004\,10^{-7},
 \qquad 0\le x\le3.\quad}                                      \tag{1}
\]

Here \(\ell=10^{-10}\), \(H(x)=\ell D(\ell x)\), and the explicit comparison
waveform is

\[
\widehat H(x)=Q'(x)+\ell e^{-\ell x}Q(x),\qquad
Q(x)=\begin{cases}
e^{4-1/[(x-1)(2-x)]},&1<x<2,\\0,&\text{otherwise}.
\end{cases}                                                     \tag{2}
\]

The exact rational \(\eta\) in (17) is approximately
\(1.000000375300181\,10^{-7}\); its rational value, not this decimal
display, makes the decision. The certificate encloses the actual coupled
master's shell output throughout the interval. It does not depend on
sampled residuals, a numerical background trajectory, an energy-to-trace
inference, or a putative lower bound for an auxiliary zero extension.

The response is exactly zero on \(0\le x\le1\). On the post-pulse interval
\(2\le x\le3\), a sharper enclosure is available:

\[
\boxed{|H(x)|\le\eta_{\rm tail}<4752/10^{17}
 =4.752\,10^{-14}.}                                             \tag{3}
\]

No identification of model coordinates with an observed physical time,
effective nonlinear amplitude radius, or finite nonlinear waveform is made.
The physical leading scalar is recovered using the retained factor
\(S(t)=3\cosh^{-3/2}(t)D(t)/(\rho_b^3p_b)\).

## 1. Actual equation, initial state and domain

The frozen source supplies the actual regular near-\(\phi=+1\) branch at
\(\delta=1/1000\), exact
\(c=96714321425692/161839258930731\), and the coupled reconstruction

\[
w_{tt}-w_{zz}+W(t,z)w=0,\qquad
-w_z+\mu w_{tt}+\beta(t)w=0\quad(z=0),                           \tag{4}
\]
\[
\mu=(\rho_bB_b)^{-1},\quad\gamma=C_b/2,\quad
\beta=\mu(35\operatorname{sech}^2t/4-25/4)+\gamma,
\quad D=w_z(t,0)-\gamma w(t,0).
\]

The conjugated data are exactly
\(w(0,z)=q_\ell(z)=Q(z/\ell)\),
\(w_t(0,z)=e^{-z}q_\ell(z)\). Both components are retained.
The initial bump is smooth, flat at both endpoints, has support contained in
\([\ell,2\ell]\), and satisfies \(0\le q_\ell\le1\).
It has value one at \(3\ell/2\).

Use

\[
T=3\ell,\quad L=8\ell,\quad
\Omega=\{(t,z):0\le t\le T,\ 0\le z,\ t+z\le L\}.
\]

All incoming backward characteristics from a point of \(\Omega\) reach the
initial line inside \([0,L]\). Outgoing ones reach either that initial line
or the physical shell; reflected source points remain inside \(\Omega\).
There is no artificial boundary condition at \(z=L-t\). The frozen analytic
collar supplies the bounds on its full depth \(0\le z\le8\ell\):

\[
\mu\ge m=1/1000>0,\quad |W|\le M=10^6,\quad
|\beta|\le B=1000,\quad |\gamma|\le\Gamma=250.                   \tag{5}
\]

Only these true-coefficient bounds are used below. In particular no
approximate coefficient is substituted into a boundary residual.
Smoothness/existence is inherited from the retained compatible evolution
hypothesis. The bump vanishes in a neighborhood of the shell corner, so
the zero shell data \(q_b(0)=q_b'(0)=0\) are compatible to all orders.

The linear domain only requires \(T<L\), here with margin \(L-T=5\ell\).
The retained example of a nonlinear causal geometry has
\(A_0=\ell/2\), \(L_{\rm corr}=3\ell\),
\(T_{\rm evol}=7\ell/2\), \(c_*=5/4\), \(Z=8\ell\). Its exact strict
margins are \(T_{\rm evol}-T=\ell/2\) and
\(Z-(L_{\rm corr}+c_*T_{\rm evol})=5\ell/8\). These inequalities do
not construct the nonlinear localized inverse or replace its hypotheses.

## 2. Incoming characteristic and the positive boundary ODE

Put \(q_b(t)=w(t,0)\), \(v(t)=q_b'(t)\), and \(J=w_t+w_z\). Then

\[
(\partial_t-\partial_z)J=-Ww,
\qquad J(t,0)=j_0(t)+R_J(t),                                   \tag{6}
\]
\[
j_0(t)=q_\ell'(t)+e^{-t}q_\ell(t),\qquad
R_J(t)=-\int_0^tW(s,t-s)w(s,t-s)\,ds.
\]

Define the finite suprema
\(A=\sup_\Omega|w|\), \(V=\sup_{[0,T]}|v|\). Thus

\[
|R_J(t)|\le MtA,\qquad |q_b(t)|\le tV.                         \tag{7}
\]

The boundary sign in (4) gives the **damped** equation

\[
\mu v'+v+\beta q_b=J(t,0),\qquad q_b'=v,\qquad v(0)=0.         \tag{8}
\]

Its exact variation-of-constants formula is

\[
v(t)=\frac1\mu\int_0^t e^{-(t-s)/\mu}
 [q_\ell'(s)+e^{-s}q_\ell(s)+R_J(s)-\beta(s)q_b(s)]\,ds.
                                                                    \tag{9}
\]

Positive inertia and the positive damping kernel matter quantitatively.
Integrating **only the known derivative input** by parts yields

\[
\frac1\mu\int_0^t e^{-(t-s)/\mu}q_\ell'(s)ds
=\frac{q_\ell(t)}\mu-
 \frac1{\mu^2}\int_0^t e^{-(t-s)/\mu}q_\ell(s)ds.               \tag{10}
\]

Both terms on the right belong to \([0,1/\mu]\), because
\(0\le q_\ell\le1\) and
\(\int_0^t e^{-(t-s)/\mu}ds=\mu(1-e^{-t/\mu})\le\mu\).
Consequently their **difference** has absolute value at most \(1/\mu\).
This is where the large derivative amplitude \(O(1/\ell)\) is avoided.
The remaining terms in (9), using (7), satisfy

\[
V\le\frac{1+T}{m}+\frac{MT^2}{2m}A+
                     \frac{BT^2}{2m}V.                       \tag{11}
\]

No derivative of a coefficient, unknown trace, residual, or sampled field
has been taken. Positivity of the kernel is used for an inequality, not to
assert positivity of the signed waveform.

## 3. A displacement bound that closes the characteristic estimate

Let \(u\) be the exact free Dirichlet wave with these initial data:

\[
u(t,z)=\tfrac12[q_o(z+t)+q_o(z-t)]
 +\tfrac12\int_{z-t}^{z+t}(e^{-\cdot}q_\ell)_o(s)ds,\qquad u_b=0.
                                                                    \tag{12}
\]

Subscript \(o\) means the smooth odd extension of the **whole positive-side
function**. In particular \((e^{-\cdot}q_\ell)_o(s)\) on the negative
side is \(-e^{s}q_\ell(-s)\). It is not \(e^{-s}q_o(s)\).
The data have magnitude at most one, so \(|u(t,z)|\le1+t\).

The half-line d'Alembert formula for (4), with the true source
\(F=-Ww\) and physical boundary value \(q_b\), is

\[
w(t,z)=u(t,z)+\mathbf1_{t\ge z}q_b(t-z)
 +\tfrac12\int_0^t\int_{z-(t-s)}^{z+t-s}F_o(s,\xi)\,d\xi\,ds.
                                                                    \tag{13}
\]

Here \(F_o\) is the odd reflection in space. Each source argument obeys
\(s+|\xi|\le t+z\le L\), so (5) applies. The integration region has
absolute-weight area at most \(t^2/2\). Formula (13) is justified by
the standard zero-Dirichlet source formula plus the boundary traveling
wave; continuity also suffices where the reflected source has a jump.
It is an identity for the retained smooth solution, not an assertion that
the odd extension solves the physical boundary problem. Hence

\[
A\le1+T+TV+\frac{MT^2}{2}A.                                   \tag{14}
\]

Set

\[
a=1+T,\quad r=MT^2/2,\quad s=BT^2/(2m),\quad
\Delta=(1-r)(1-s)-Tr/m.
\]

For the specified constants,
\(r=s=9/200000000000000\), and

\[
\Delta=
\frac{39999999999996399999460000081}
     {40000000000000000000000000000}>0.
\]

The matrix of (14), (11) is
\(\left(\begin{smallmatrix}1-r&-T\\-r/m&1-s\end{smallmatrix}\right)\).
Its inverse is entrywise positive. Applying it to the componentwise
inequalities gives

\[
A\le A_*:=\frac{a(1-s+T/m)}\Delta
 =\frac{40000012011998203599999460000}
        {39999999999996399999460000081}<1.00000031,
\]
\[
V\le V_*:=\frac a{m\Delta}
 =\frac{40000000012000000000000000000000}
        {39999999999996399999460000081}<1000.000001.             \tag{15}
\]

These are upper bounds on the true solution; the equality solution of the
majorant system is not a numerical approximation to the master field.

## 4. Uniform output error and support information

Differentiate the explicit comparator at the shell. The odd extensions
give exactly

\[
\widehat D(t)=u_z(t,0)=q_\ell'(t)+e^{-t}q_\ell(t)=j_0(t).
\]

Since \(w_z=J-v\), equations (6) and (7) imply the derivative-trace
bound directly:

\[
|D(t)-\widehat D(t)|
 =|R_J(t)-v(t)-\gamma q_b(t)|
 \le MT A_*+(1+\Gamma T)V_*.                                  \tag{16}
\]

Multiplication by \(\ell\), with \(t=\ell x\), gives (1), where

\[
\eta=\ell[MT A_*+(1+\Gamma T)V_*]
 =\frac{20000007506001802249730539999919}
        {199999999999981999997300000405000000000}.
                                                                    \tag{17}
\]

This has the expected leading scale \(\ell/m\) as
\(\ell\to0\), for fixed coefficient bounds and \(T=3\ell\).
For example the next potential term is \(3M\ell^2\).
The proof controls the derivative observation through its incoming
characteristic and the boundary velocity. It does not infer (16) from a
small bulk displacement error.

Initial data vanish on \([0,\ell]\). The equation and boundary system are
homogeneous, and the backward characteristic region \(t+z\le\ell\)
has only those zero data. The usual characteristic Volterra uniqueness,
or the same inequalities (11), (14) with zero free input on that region,
gives \(w=0\) there. Continuity includes its edge. Hence
\(H=\widehat H=0\) on \([0,1]\). This proves the pre-arrival zero interval;
it does not identify an exact nonlinear arrival time.

For \(2\ell\le t\le T\), \(q_\ell(t)=0\). The memory integral in
(10) remains. Its integration support has length \(\ell\), so (9) yields

\[
|v(t)|\le V_{\rm tail}:=
\frac\ell{m^2}+\frac\ell m+\frac rm A_*+sV_*.
\]

Since \(\widehat D=0\) on this interval,

\[
\eta_{\rm tail}=\ell[MT A_*+V_{\rm tail}+\Gamma T V_*]
 =\frac{95020036022498208999324730040419}
        {1999999999999819999973000004050000000000000000}
 <4752/10^{17}.                                               \tag{18}
\]

The tail is not asserted to vanish. Replacing the dynamic shell by a
Dirichlet shell would lose precisely this boundary memory and the true
boundary velocity.

The support information also sharpens propagation to a response norm:

\[
\|H-\widehat H\|_{L^2(0,3)}^2\le\eta^2+\eta_{\rm tail}^2,
\qquad
\|H-\widehat H\|_{L^2(0,3)}
 <\zeta:=10000005/10^{14}=1.0000005\,10^{-7}.                   \tag{19}
\]

There is no contribution on \([0,1]\), and each of \([1,2]\),
\([2,3]\) has length one. A certified comparator norm
\(\|\widehat H\|_2\in[N_-,N_+]\) therefore yields
\(\|H\|_2\in[\max(0,N_--\zeta),N_++\zeta]\).
The retained observable identity becomes
\(\mathcal C_T=(K/\ell)\|H\|_2^2\),
\(K=3\pi^2/(\rho_b^2p_b^2)>0\). Its background interval and any
comparator quadrature uncertainty must also be propagated.

## 5. Arithmetic checks and independent failure controls

`check_uniform_waveform.py` reads the frozen proof, prior Sol derivation and
background bounds as bytes, authenticates their SHA-256 pins, and independently
evaluates (14)--(18) with `fractions.Fraction`. It imports no archived
producer. Its decisions use only integers and fractions and remain active
under `python -O`. It records exact fractions and checks the full-domain
and causal margins. Inputs reside at the capsule-relative path
`inputs/previous`, with unchanged source hashes. The script's output records
which inputs were checked.

The following finite controls detect specified wrong formulas:

* A polynomial boundary jet checks the damping sign in (8), the sign of
  \(\beta q_b\), the identity \(w_z=J-v\), and both the sign and the factor
  \(\gamma=C_b/2\) in the physical observation.
* Removing the coupling term \(Tr/m\) from the determinant produces
  numbers that fail the supersolution inequalities (11), (14).
* Setting the velocity-input contribution to zero changes (2) at its peak
  by a strictly positive amount, certified using
  \(e^{-3\ell/2}\ge1-3\ell/2>0\).
* The memory term in (10) cannot be deleted when \(q_\ell(t)=0\).
  For the diagnostic input \(f(s)=s(1-s)\) on \([0,1]\), zero outside,
  and \(\mu=1,t=2\), it is strictly negative, with magnitude greater
  than \(1/48\). The lower bound follows from \(e^2<8\) and
  \(\int_0^1f=1/6\). This control concerns the boundary ODE, not a
  replacement physical state.

These checks support the finite algebra and rational inequalities. The
analytic argument above supplies the universal characteristic/kernel
estimate. Neither these checks nor agreement with another derivation
certifies the retained nonlinear hypotheses.

Five explicit failure injections are available through `--inject-error`:
`drop-determinant-coupling`, `flip-boundary-damping`, `use-full-C`,
`omit-velocity-input`, and `delete-boundary-memory`. Each is required to
exit unsuccessfully with a stable named failed check. Injected runs write no
passing output artifact.

Reproduce using the new standard-library checker:

```sh
python -I -B sol/check_uniform_waveform.py
python -I -O -B sol/check_uniform_waveform.py --output UNIFORM_WAVEFORM_CHECKS_OPTIMIZED.json
```

The specified leading response coefficient can be tightened separately by
validated integration of \(\widehat H^2\) and propagation of (19).
Those quadrature results are not inputs to this uniform certificate.
