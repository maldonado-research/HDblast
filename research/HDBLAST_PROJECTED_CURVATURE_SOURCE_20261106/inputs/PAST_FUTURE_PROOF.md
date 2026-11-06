# A past–future estimate for a registered timelike-shell packet

4 November 2026, America/Los_Angeles. Conditional mathematical working note.

For the registered quadrupole with width `ell=10^-6`, let the embedding-included
first scalar on the shell be

\[
\widehat\Phi_1(\tau,\Omega)=S(\tau/\rho_b)(n_1^2-n_2^2).
\]

Under the retained branch, shell, reconstruction and smooth linear development
premises described below, this note establishes

\[
\boxed{\sup_{0\le t\le3\ell}|S(-t)-S(t)|
 < B:=\frac{264000000001}{13738809600}<20.}                 \tag{1}
\]

Combining (1) with two **actual**, same-state future enclosures gives

\[
S(-5\ell/4)<-37376,\qquad S(-7\ell/4)>37379.               \tag{2}
\]

Consequently the continuous first scalar coefficient has at least one zero in
`(-7ell/4,-5ell/4)`. This is a coarse zero bracket; neither uniqueness nor the
older narrow future zero brackets or late-tail sign follows from this estimate.
The number 20 bounds a difference in the fixed model units. It is not an
individual scalar-amplitude bound, relative error or finite nonlinear signal.

## 1. State, observation and prior premises

Both sides concern one fixed Einstein–scalar background: the registered regular
near-phi=+1 branch at `delta=1/1000` and
`c_model=96714321425692/161839258930731`. Set

\[
Q(r)=\begin{cases}\exp(4-1/[(r-1)(2-r)]),&1<r<2,\\0,&\text{otherwise}.
\end{cases}
\]

The longitudinal master is the same actual solution `w` with

\[
w(0,z)=Q(z/\ell),\qquad w_t(0,z)=e^{-z}Q(z/\ell),
\]
\[
w_{tt}-w_{zz}+\mathcal W(t,z)w=0,\qquad
-w_z+\mu_b w_{tt}+\beta(t)w=0\quad(z=0),                  \tag{3}
\]
\[
\beta(t)=\gamma+\mu_b(35\operatorname{sech}^2t/4-25/4),
\quad\gamma=C_b/2,\quad \mu_b=(\rho_b B_b)^{-1}.
\]

The earlier [past-contact proof](inputs/CONTACT_PROOF.md) supplies the positive
endpoint-mass operator, all-power-domain regularity, local finite propagation
and endpoint physical trace reconstruction. Its quantitative headline concerns
`-3ell<=t<=0`. Section 2 below uses the even equation to construct the positive
elapsed developments with both velocity signs and identify the supplied future
state. Branch existence, validity of the analytic enclosures, coupled
reconstruction and linear uniqueness remain premises. They are not reproved by
source authentication or finite arithmetic.

The older [future waveform proof](inputs/OLD_FUTURE_PROOF.md), sections 1–3,
specifies exactly (3), the same `ell=10^-6`, the same nonzero initial velocity,
the quadrupole `Y=n1^2-n2^2`, and the same embedding-included observation:

\[
S(t)=\frac{3\cosh^{-3/2}t}{\rho_b^3p_b}
       [w_z(t,0)-\gamma w(t,0)],\qquad \tau=\rho_b t.      \tag{4}
\]

The constant background shell scalar removes first clock/label pullback
amplitude terms. The spatial normal and the straight-shell reconstruction are
the retained ones, with no time-reflection sign change in (4). A sign of `S`
is a sign of its coefficient: the angular factor itself changes sign and has
zero angular mean.

The selected old [actual report](inputs/OLD_FUTURE_REPORT.json) pins the same
raw shell box and wider background documents included here. Its actual-minus-
comparator enclosure proof is accepted as a prior conditional certificate.
This note uses its pointwise `ell_S` witnesses, which already include the actual
error enlargement. An `H0` comparator enclosure alone would not suffice. The
older `ell=10^-10` state is a different tangent and is not used.

## 2. Time reflection with the complete dynamic boundary

Use scaled positive elapsed time `s` and radius `r`, and define
`U(s,r)=w(ell*s,ell*r)`, `q(s)=U(s,0)`. On a common auxiliary half-line, the
retained operator is

\[
\mathcal H=L^2(0,\infty)\oplus\mathbb R,
\qquad\|(u,q)\|_{\mathcal H}^2=\int u^2dr+\Lambda q^2,
\quad\Lambda=\mu_b/\ell>1000,
\]
\[
A_0(u,q)=(-u_{rr}+V_0u,\ [-u_r(0)+b_0q]/\Lambda),
\quad q=u(0),\quad A_0\ge\kappa I,\quad\kappa=5500\ell^2,
\]
\[
V_0=\ell^2\rho^2\mathcal Q_{\rm bg}\ \ (0\le r\le6),\quad
b_0=\ell(\gamma+5\mu_b/2),\quad
\omega(s)=-\frac{35}{4}\ell^2\tanh^2(\ell s).
\]

Its closed form domain is the linked `H1` trace domain; its operator domain is
the linked `H2` domain, with no additional static Robin condition. The master
evolution and its **full** boundary row are

\[
Z_{ss}+(A_0+\omega(s)I)Z=0,
\qquad -U_r+\Lambda q_{ss}+[b_0+\Lambda\omega(s)]q=0.     \tag{5}
\]

Beyond `r=6`, use the same retained smooth positive extension, interpolated on
`6<=r<=7` and constant `6000ell^2` for `r>=7`. It preserves the static coercivity;
it is not a statement that the actual coefficients are continued globally by
that extension.

The temporal coefficient, boundary spring and `cosh` factor are even. Thus the
actual past solution `w(-ell*s,ell*r)` is, by the same linear uniqueness, the
forward solution with data `(Q,-g)`, where `g=ell*e^(-ell*r)Q`. The actual future
has data `(Q,+g)`. This reverses the supplied nonzero velocity; it does not
silently set that velocity to zero.

Both developments have support in `[0,2+s]` for `0<=s<=3`, and in
`[1-s,2+s]` before `s=1`. Their actual fields therefore stay within `r<=5`,
where the common auxiliary operator equals the retained actual coefficients
through `r<=6`, inside the validated collar `r<=8`. The exterior extension used for global
operator estimates does not alter these local solutions.

Write `Delta Z=Z_minus-Z_plus` and `delta q=q_minus-q_plus`. Then

\[
\Delta Z(0)=0,\qquad \Delta Z_s(0)=(-2g,0).               \tag{6}
\]

The factor two in (6) is essential. Smooth flat initial data satisfy the linked
corner conditions and every required power-domain condition, as in the retained
operator construction.

## 3. Complete basic and graph energies

The complete positive form energy includes the endpoint mass and spring:

\[
E_\Delta=\tfrac12\|\Delta Z_s\|_{\mathcal H}^2
 +\tfrac12\langle A_0\Delta Z,\Delta Z\rangle_{\mathcal H}
 +\tfrac12\omega\|\Delta Z\|_{\mathcal H}^2.
\]

Since `omega_s<=0` on positive elapsed time, `E_Delta` is nonincreasing and

\[
E_\Delta(s)\le E_\Delta(0)=2\|g\|_2^2\le2\ell^2.         \tag{7}
\]

The mass gives `|delta q_s|^2<=E_Delta(0)/500`. There is no boundary response
before `s=1`, so for `0<=s<=3`

\[
|\delta q(s)|^2\le[(s-1)_+]^2E_\Delta(0)/500
 \le E_\Delta(0)/125\le2\ell^2/125<\ell^2.              \tag{8}
\]

For `W=A0^(1/2)Delta Z`, use the abstract graph energy

\[
E_{1/2}=\tfrac12\|A_0^{1/2}\Delta Z_s\|_{\mathcal H}^2
 +\tfrac12\|A_0\Delta Z\|_{\mathcal H}^2
 +\tfrac12\omega\|A_0^{1/2}\Delta Z\|_{\mathcal H}^2.
\]

The scalar temporal coefficient commutes with `A0`. Graph regularity and the
spectral form identity justify differentiation and retain every endpoint term.
There is no compact-support claim for the generally nonlocal vector `W`.
Again `E_1/2` is nonincreasing. Since `tanh(ell*s)<=ell*s` and `s<=3`,

\[
E_{1/2}(s)\ge\frac{1-\epsilon}{2}\|A_0\Delta Z\|^2,
\qquad\epsilon=\frac{63\ell^2}{4400}<1.                  \tag{9}
\]

The initial displacement in (6) vanishes exactly. The retained unit-bump
estimates `||Q||2<=1`, `||Q_r||2<=16` and the bound `V0<103000ell^2` give

\[
2E_{1/2}(0)=4a_0(g,g)
\le4\ell^2[(16+\ell)^2+103000\ell^2],                   \tag{10}
\]

because `g_r=ell*e^(-ell*r)(Q_r-ell Q)` and its initial endpoint trace is zero.
The identity `||A0^(1/2)g||^2=a0(g,g)` includes all endpoint components; no
endpoint of a half-power vector has been discarded. From (9)–(10), exact
rational arithmetic yields

\[
\|A_0\Delta Z\|_{\mathcal H}^2
\le\frac{4\ell^2[(16+\ell)^2+103000\ell^2]}{1-\epsilon}
 <(33\ell)^2.                                          \tag{11}
\]

These positive energies are estimating norms. No physical Noether or injected
energy identification, thermal arrow or entropy law is asserted.

## 4. Physical difference bound

The actual endpoint/conjugation cancellation in the retained reconstruction is

\[
S(-\ell s)-S(\ell s)=
\frac{3\mu_b\cosh^{-3/2}(\ell s)}{\rho_b^3p_b}
\left[-\frac{(A_0\Delta Z)_b}{\ell^2}+\frac52\delta q\right]. \tag{12}
\]

No extra minus sign appears: the derivative in (4) is spatial. From the
endpoint mass in `H`,
`mu_b*|(A0 Delta Z)_b|<=sqrt(mu_b*ell)*||A0 Delta Z||_H`.
The retained shell guards are `mu_b<1/400`, `sqrt(mu_b*ell)<1/20000`,
`rho_b>=129` and `|p_b|>=3/25000`. Using (8) and (11) in (12) gives

\[
|S(-\ell s)-S(\ell s)|
 <\frac{25000}{129^3}\left[\frac{33}{20000\ell}+\frac\ell{160}\right]
 =\frac{264000000001}{13738809600}<20,                    \tag{13}
\]

which proves (1). It pairs proper times `-rho_b*t` and `+rho_b*t` at one anchor.
The bound also holds for the angular field difference since `|Y|<=1`.

## 5. Transport two retained actual future witnesses

The older report exports **ell times S**, not `S`. Divide its `ell_S` endpoint
intervals by `ell` before widening them by `[-B,B]`. For example, at `x=5/4`
the future input is exactly

\[
\ell S(5\ell/4)\in
\left[
\frac{-214137780846555809997409519084451}{5192296858534827628530496329220096},
\frac{-194172220795747890795288722661291}{5192296858534827628530496329220096}
\right].
\]

At `x=7/4` it is exactly

\[
\ell S(7\ell/4)\in
\left[
\frac{194184829173661869937336484226623}{5192296858534827628530496329220096},
\frac{53537921869770098025334130218115}{1298074214633706907132624082305024}
\right].
\]

For any retained actual future interval `[L,U]`, (1) implies the conservative
past interval `[L-B,U+B]`. Applying this to those divided endpoints proves (2).
Readable outward enclosures are

| Coordinate time | Conservative enclosure for past S |
|---|---:|
| `-5ell/4` | `[-41261,-37376]` |
| `-7ell/4` | `[37379,41264]` |

The full Fraction endpoints and the strict integer sign margins are computed
in the canonical report. Continuity then gives a zero between these two past
times, or between proper times `-7rho_b*ell/4` and `-5rho_b*ell/4`.

The old central witnesses at `x=1.49998,1.49999`, and late witnesses at
`x=1.96,1.97`, each widen to an interval containing zero. Their finer old zero
brackets therefore do not transfer by (1). Neither does the old negative
continuous tail sign, whose physical-S magnitude is below this difference
ceiling. The older response integral is quadratic, so it needs a separate norm
and cross-term estimate; a pointwise transport alone does not certify it.

## 6. Contribution, scope and reproduction

Time reflection of an even linear equation, positive graph-energy estimates and
interval widening are standard methods. The positive endpoint-mass operator,
physical trace reconstruction and future waveform enclosures are earlier
results credited above. This note combines them for one registered packet to
give an explicit absolute past–future difference bound and a coarse past sign
change. It claims no new general principle or exhaustive novelty assessment.

Earlier work already constructs a **conditional formal complete second Taylor
coefficient on a future contact interval**; see
[the prior second-order contact theorem](https://github.com/maldonado-research/HDblast/blob/f4001e5d94e6388bfff4f82dee63c174f2dfd5eb/research/HDBLAST_FULL_SECOND_ORDER_CONTACT_20261014/FULL_FINITE_ANGULAR_SECOND_ORDER_CONTACT.md).
That statement retains localized constraint/gauge/corner hypotheses. The
present first-coefficient estimate supplies neither a quantitative particular
past second-order baseline/observable nor an exact nonlinear development with
an effective parameter interval. It does not identify a remote source,
radiation, SI scale or a higher-dimensional cosmological origin.

Seven selected public input documents are included and byte/Git-blob pinned in
`SOURCE_PINS.json`. The older pointwise enclosure theorem is a retained premise:
its producer, complete campaign, gain proof and continuous-cell CSV are not
rerun or claimed newly closed here. The CSV is unnecessary for the two selected
report witnesses used in this note. Authentication does not prove the inherited
analytic theorems or their scientific significance.

`verify.py --output /absolute/fresh/external/report.json` uses only the Python
standard library and exact `Fraction` arithmetic. Run it once normally and once
with `-O`; their report bytes must agree. It authenticates the selected inputs,
checks state and normalization joins, evaluates (7)–(13), transports the selected
actual witnesses and checks the specified non-transfers. It invokes no old
producer, project scientific import, background/field/source/bump oracle or
numerical PDE solver. The continuous proof and its explicit premises remain
necessary.
