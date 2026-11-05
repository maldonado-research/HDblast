# Scalar endpoint realization, closed ground-state identity, and a fixed-background gap

13 October 2026. New mathematical continuation for internal review. Historical
files and the frozen public capsule are unchanged. This note proves an abstract
operator theorem for the scalar Sturm–Liouville expression supplied by the
October 11 prospective bridge. The separate identification of this operator
with all physical perturbations, including exceptional harmonics and the
kinetic norm, is **not** an assumption silently discharged here.

## 1. Exact hypotheses and notation

Fix one finite background endpoint `T>0`, hence one positive detuning when the
registered branch is used. Write the background momentum as `p=eta'`, and write

\[
 \alpha=\frac1{\rho^2p^2},\qquad
 q=\frac2{3\rho^2},\qquad
 w=\frac1{\rho^4p^2},\qquad
 \tau u=w^{-1}[-(\alpha u')'+qu].                 \tag{1}
\]

Assume:

1. `rho,p` are real and `C^2` on `(0,T]`, `rho>0`, and `p` never vanishes
   there. All coefficients are regular at `T`.
2. For some real `A_c != 0`, the regular-cone asymptotics are
   `rho=y(1+O(y^2))`, `p=A_c y(1+O(y^2))`. For the Liouville/essential-spectrum
   assertions only, require these remainders through two logarithmic
   derivatives: applying `(y d/dy)^j`, `j=0,1,2`, to each relative remainder
   still gives `O(y^2)`. Smooth regular-cone series imply this requirement.
3. A real strictly positive function `P_0` solves `tau P_0=4P_0` and
   `P_0=C_0 y^4(1+O(y^2))`, with `C_0>0`. Its leading coefficient may be
   arbitrarily small. The displayed cone expansion of the Jacobi reference
   has `C_0=-U'(u)/60` on the registered negative-momentum branch.
4. The boundary number `beta>0` is fixed and finite, and
   \[
    r_0=\frac{\alpha(T)P_0'(T)}{P_0(T)}-\frac4\beta>0.       \tag{2}
   \]
   For the supplied scalar boundary problem,
   `beta=rho(T)^4 p(T)^2 B_partial`.

All differentiability is classical away from the cone. Local absolute
continuity is sufficient in the domains below. The leading cone asymptotics,
without differentiated remainder bounds, suffice for the form, gap, and
limit-point assertions. The additional differentiated bounds establish the
stated Liouville potential limit and essential spectrum.

The October 11 public Jacobi signs and the private bridge supply these signs
on every nonconstant registered interior: `p<0`, `P_0>0`,
`B_partial>J>=w_model/2>0`, and `r_0>0`. Here `w_model=2` is a model constant,
not the weight `w` in (1). This application remains conditional on the bridge's
perturbation reduction; the analysis below starts at (1)–(2).

## 2. The Hilbert space and closed form

Use the positive Hilbert space

\[
 \mathscr H=L^2((0,T),w\,dy)\oplus\mathbb C,\qquad
 \langle(u,z),(v,t)\rangle_{
 \mathscr H}=\int_0^T w u\bar v\,dy+\beta^{-1}z\bar t.       \tag{3}
\]

The two components of a general Hilbert-space vector are independent. Trace
identification is imposed on the **form domain**, not on every vector in the
Hilbert space. Define

\[
 \begin{split}
 V=\{(u,u(T)):\;&u\in AC_{\rm loc}((0,T]),\quad
       \int_0^T w|u|^2\,dy<\infty,\quad
       \int_0^T\alpha|u'|^2\,dy<\infty\},\\
 \mathfrak e[(u,u(T))]&=\int_0^T(\alpha|u'|^2+q|u|^2)\,dy.  \tag{4}
 \end{split}
\]

The integral with `q` is automatically finite: `q/w=(2/3)rho^2 p^2` is bounded
on `(0,T]` and tends to zero at the cone. The trace at `T` is well defined by
ordinary `H^1` regularity in a neighborhood of the regular endpoint.

**Density.** Smooth functions supported inside `(0,T)` approximate any bulk
`L^2(w)` vector. Given any desired independent endpoint component, add a smooth
bump supported in `[T-epsilon,T]` with that endpoint value. Its bulk norm tends
to zero with `epsilon`. Each bump has finite form energy, which need not remain
bounded for this Hilbert-space density argument. Thus `V` is dense in
`mathscr H`; the trace condition does not contradict density.

**Closedness.** A sequence Cauchy in `mathfrak e+||.||_mathscr H^2` converges in
`L^2(w)` and its weighted derivatives converge in `L^2(alpha)`. On each compact
interval `[epsilon,T]` it converges in ordinary `H^1`; these local limits agree
with the bulk limit. Continuity of the local trace identifies the limiting
endpoint component with `u(T)`. This proves completeness and closedness.

**A cone core.** Let `V_c` consist of smooth functions on `(0,T]` which vanish
on a neighborhood of zero, with their endpoint traces and with no condition
at `T`. This is a form core. Indeed choose a cutoff `chi_epsilon=0` below
`epsilon`, `chi_epsilon=1` above `2epsilon`, and
`|chi_epsilon'|<=C/epsilon`. The only new derivative term is bounded by

\[
 \int_\epsilon^{2\epsilon}
       \alpha|\chi_\epsilon'|^2|u|^2
 \le C_1\int_\epsilon^{2\epsilon}w|u|^2\longrightarrow0,   \tag{5}
\]

because `alpha/w=rho^2=O(y^2)`. All other errors are tails of integrable
quantities. Ordinary smooth approximation away from zero finishes the proof.
This establishes the endpoint closure rather than extrapolating it from
Frobenius eigenfunctions.

## 3. The self-adjoint operator and the eigenparameter boundary condition

The representation theorem for densely defined closed nonnegative forms gives
one nonnegative self-adjoint operator `A` associated with (4). Explicitly,

\[
 \begin{split}
 D(A)=\{(u,u(T))\in V:\;&\alpha u'\in AC_{\rm loc}((0,T]),
                         \quad \tau u\in L^2(w\,dy)\},\\
 A(u,u(T))&=(\tau u,\ \beta\alpha(T)u'(T)).                \tag{6}
 \end{split}
\]

To verify (6), test the representation identity first against compactly
supported functions to obtain the bulk distributional equation and local
regularity. Testing near `T` identifies the second component as the positive
outward flux in (6). Conversely, integration by parts for core tests, followed
by the form-core approximation (5), proves the representation identity for
all `V`. No evaluation of an unjustified cone boundary product is needed.

This is also a literal Friedrichs construction. On `V_c` define the symmetric
operator `S(u,u(T))=(tau u,beta alpha(T)u'(T))`; its quadratic form is (4) on
that core. Its Friedrichs extension is `A`. Taking compact support at both
ends instead would erroneously remove the allowed shell trace and would not
produce this endpoint-mass problem.

The eigenvalue equation `A(u,z)=Lambda(u,z)` now gives

\[
 \tau u=\Lambda u,\qquad z=u(T),\qquad
 \alpha(T)u'(T)=\Lambda u(T)/\beta.                       \tag{7}
\]

Thus the eigenvalue-dependent boundary condition is encoded by an ordinary
self-adjoint operator in the enlarged Hilbert space. A vector `(0,z)` with
`z!=0` is a legitimate Hilbert vector, but it is not in the form domain and
cannot be an independent boundary eigenmode. Complex spectral values are
excluded by self-adjointness, not merely by a calculation on real trial modes.

**Limit-point cone.** At `Lambda=4`, a second independent bulk solution is

\[
 P_1(y)=P_0(y)\int_{y_*}^{y}\frac{dt}{\alpha(t)P_0(t)^2}.
                                                               \tag{8}
\]

Since `alpha P_0^2 ~ (C_0^2/A_c^2)y^4`, this solution behaves as a nonzero
multiple of `y`, while `P_0~C_0 y^4`. Only `P_0` is in `L^2(w dy)`:
`w P_0^2~(C_0^2/A_c^2)y^2`, whereas `w P_1^2` is proportional to `y^-4`.
The Weyl limit-point/limit-circle alternative therefore classifies zero as
limit-point. No independent self-adjoint cone boundary parameter is available.
The form construction specifies precisely the finite-energy realization.

## 4. Ground-state transform on the full closed form domain

Put

\[
 a_g=\alpha P_0^2,\qquad m_g=wP_0^2,\qquad
 b_g=P_0(T)^2/\beta,\qquad k_g=r_0P_0(T)^2>0.              \tag{9}
\]

For `u=P_0 f`, the Hilbert norm and shifted energy are

\[
 \begin{split}
 N[f]&=\int_0^T m_g|f|^2\,dy+b_g|f(T)|^2,\\
 D[f]&=\int_0^T a_g|f'|^2\,dy+k_g|f(T)|^2.               \tag{10}
 \end{split}
\]

For every core function, integration of `tau P_0=4P_0` gives the exact complex
sesquilinear/real quadratic Picone identity

\[
 \boxed{\mathfrak e[u]-4\|(u,u(T))\|^2_{\mathscr H}=D[u/P_0].} \tag{11}
\]

The finite-endpoint term is
`[alpha P_0'/P_0-4/beta]|u(T)|^2`, with precisely the sign in (2).
There is no cone term for a core function.

Equation (11) passes to the **entire** closed form domain. A core sequence
Cauchy in `mathfrak e+N` is, by the core identity, Cauchy in `D+5N` after
multiplication by `1/P_0`. Local convergence identifies its weighted derivative
with `(u/P_0)'`. Conversely, the natural transformed domain

\[
 W=\{f\in AC_{\rm loc}((0,T]):
       \int m_g|f|^2+\int a_g|f'|^2<\infty\}            \tag{12}
\]

has the same cone-core approximation: `a_g/m_g=alpha/w=rho^2` gives (5) again.
Multiplication by `P_0` maps `W` onto `V`, preserves the norm (10), and
identifies the closed forms by (11). This proves closure of the identity
without assuming that arbitrary form vectors have eigenfunction expansions
at the cone. No lower bound on `|A_c|` is used.

Strict positivity of `D[f]` for nonzero `f` follows immediately: its vanishing
would force `f` constant and then `f(T)=0`. Alone this observation would rule
out only a zero eigenvector; it would not exclude zero from continuous
spectrum. The next estimate supplies the missing coercivity.

## 5. An explicit positive gap without an eigenvalue scan

The cone behavior gives

\[
 a_g\sim c_g y^4,\quad m_g\sim c_g y^2,\quad
 c_g=C_0^2/A_c^2>0.                                    \tag{13}
\]

In the Jacobi normalization `C_0=-U'(u)/60`, `A_c=U'(u)/5`, one has exactly
`c_g=1/144`. Both the small central slope and the chosen normalization cancel.
Define background-dependent constants

\[
 \alpha_*:=\inf_{0<y\le T}\frac{a_g(y)}{y^4}>0,\qquad
 m_*:=\sup_{0<y\le T}\frac{m_g(y)}{y^2}<\infty,\qquad
 M_g:=\int_0^T m_g(y)\,dy<\infty.                       \tag{14}
\]

Their finiteness and positivity follow from (13), interior positivity and
compactness after adjoining the limiting values at zero. They are defined by
one background, not by a detuning-uniform estimate.

Here is an elementary weighted inequality, valid for complex `f`. Let
`g=f-f(T)`. On `[epsilon,T]`, integration by parts gives

\[
 3\int_\epsilon^T y^2|g|^2
 =-\epsilon^3|g(\epsilon)|^2
    -2\Re\int_\epsilon^T y^3\bar g g'
 \le 2\left(\int_\epsilon^T y^2|g|^2\right)^{1/2}
          \left(\int_\epsilon^T y^4|g'|^2\right)^{1/2}.
\]

Drop the nonpositive boundary term, divide when the first integral is nonzero,
and let `epsilon` decrease to zero. Using (14),

\[
 \int_0^T m_g|f-f(T)|^2\,dy
 \le C_*\int_0^T a_g|f'|^2\,dy,
 \qquad C_*:=\frac{4m_*}{9\alpha_*}.                    \tag{15}
\]

No boundedness or trace condition on `f` at zero was used. The inequality holds
directly on (12), and also follows by closure from the core. In
`L^2(m_g) oplus C_(b_g)`, the triangle inequality and Cauchy–Schwarz imply

\[
 \sqrt{N[f]}\le\sqrt{C_*\int a_g|f'|^2}
                   +\sqrt{M_g+b_g}|f(T)|
 \le\sqrt{\left(C_*+\frac{M_g+b_g}{k_g}\right)D[f]}.
\]

Consequently,

\[
 \boxed{\quad A-4I\ \ge\ \varepsilon_\delta I,
 \qquad
 \varepsilon_\delta=
 \left[\frac{4m_*}{9\alpha_*}+\frac{M_g+b_g}{k_g}\right]^{-1}>0.\quad} \tag{16}
\]

An entirely sufficient variant replaces `M_g` by `m_* T^3/3` and uses any
proved lower bound on `alpha_*` and upper bound on `m_*`. Formula (16) is an
explicit background-functional gap bound. Evaluating a numerical rational
lower bound for a chosen background would require corresponding analytic or
validated coefficient bounds; no such numerical value is claimed here. The
existence of a strictly positive bound at every fixed background is proved.
The expression is unchanged by rescaling `P_0` by a nonzero constant.

A sharper one-dimensional option is available without a spectral scan. Put

\[
 H_g=\sup_{0<x<T}
   \left(\int_0^x m_g\right)\left(\int_x^T a_g^{-1}\right).
                                                               \tag{17}
\]

One may replace `C_*` by `4H_g`; the product tends to `1/9` at zero and to zero
at `T`, hence is finite. Also `H_g<=m_*/(9 alpha_*)`. For completeness, an
elementary proof of its Hardy inequality appears below; the power estimate
(15) already suffices for every conclusion.

**Spectral conclusion.** The spectral theorem applied to (16) gives

\[
 \sigma(A-4I)\subset[\varepsilon_\delta,\infty).          \tag{18}
\]

Thus complex spectral values and **all** spectrum at `mu^2<=0` are excluded,
including eigenvalues, continuous spectrum and residual spectrum. For a
self-adjoint operator residual spectrum is absent in any event. This is
stronger than a nodeless-mode comparison restricted to eigenfunctions.

No uniform lower bound on (16) as `delta -> 0` has been proved. The formula
makes the missing uniform estimates explicit. Endpoint asymptotics and
`r_0>0` alone are insufficient for such a claim; a small positive endpoint
penalty or an interior low-conductance region can shrink this estimate.
The divided operator at `p identically 0` is outside the hypotheses.

## 6. Liouville half-line and essential threshold

Under the differentiated cone hypotheses, define

\[
 x(y)=\int_y^T\frac{dt}{\rho(t)},\quad 0\le x<\infty,
 \qquad h(y)=(\alpha w)^{1/4}
             =\frac1{\rho(y)^{3/2}|p(y)|},\qquad \psi=h u. \tag{19}
\]

The shell is `x=0`, the cone is `x=infinity`, and `dx/dy=-1/rho`. The unitary
norm is

\[
 \|(u,u(T))\|^2_{\mathscr H}
 =\int_0^\infty|\psi|^2\,dx+d_b|\psi(0)|^2,
 \quad d_b=\frac1{\beta h(T)^2}
          =\frac1{\rho(T)B_{\partial}}>0.               \tag{20}
\]

Set

\[
 s(x)=\frac{h_x}{h}
      =\rho\left(\frac32\frac{\rho'}\rho+\frac{p'}p\right).
\]

The transformed differential expression and boundary condition are

\[
 -\psi_{xx}+V_L\psi=\Lambda\psi,\qquad
 V_L=\frac23\rho^2p^2+s_x+s^2,
 \qquad
 -\psi_x(0)+s(0)\psi(0)=\Lambda d_b\psi(0).             \tag{21}
\]

In particular, the endpoint mass remains positive. A useful coefficient-only
expression is

\[
 V_L=\rho^2\left[
 \frac34 a^2+2a\ell+\ell^2-\frac32 a'-\ell'
                  +\frac23p^2\right],
 \qquad a=\rho'/\rho,\quad\ell=p'/p.                    \tag{22}
\]

No background equation is required for (22). If the Einstein–scalar background
equations are used, it is also

\[
 V_L=-\frac52+\rho^2\left[
 \frac34a^2+6a\ell+2\ell^2-U''(\eta)-\frac16p^2\right]. \tag{23}
\]

The transformed form and operator are

\[
 \begin{split}
 \mathfrak e_L[\psi]&=\int_0^\infty
       (|\psi_x|^2+V_L|\psi|^2)\,dx+s(0)|\psi(0)|^2,\\
 A_L(\psi,\psi(0))&=
   (-\psi_{xx}+V_L\psi,\ d_b^{-1}[-\psi_x(0)+s(0)\psi(0)]).
                                                               \tag{24}
 \end{split}
\]

Its form domain is `H^1(0,infinity)` with trace identification; its operator
domain is `H^2(0,infinity)` with trace identification. To see the form-domain
statement directly, the original derivative energy is
`integral |psi_x-s psi|^2`, and `q/w` is bounded. The hypotheses make `s` and
`V_L` bounded, so this norm is equivalent to the ordinary `H^1` norm after
adding the bulk norm. Local regularity and the bounded potential then give
`H^2` for the operator domain. For an operator vector, `psi` and `psi_x` tend
to zero at infinity; hence its energy Green boundary product vanishes at the
cone. For arbitrary form vectors the closed identity (11) is the preferred
statement.

The cone expansion implies

\[
 x=-\log y+C+O(y^2),\qquad s=\frac52+O(y^2),\qquad
 V_L=\frac{25}{4}+O(y^2)
     =\frac{25}{4}+O(e^{-2x}).                          \tag{25}
\]

The bulk constant-potential Dirichlet half-line operator has essential spectrum
`[25/4,infinity)`. A bounded potential tending to zero is relatively compact
against its resolvent: compactness on each finite interval is the Rellich
embedding, and the tail operator norm tends to zero. Thus replacing its
potential by `V_L` preserves this essential spectrum.

The dynamic endpoint also does not change it. Let `D_V` be the Dirichlet
half-line operator with the **same actual potential** `V_L`. Compare (24)
with `D_V oplus c`, where `c` is any fixed real operator on the one-dimensional
endpoint space. At a common nonreal resolvent point, the difference of the
two bulk resolvent outputs solves `(-d^2/dx^2+V_L-z)psi=0` and belongs to
`L^2`; its space has dimension at most one by the limit-point property.
This finite-rank boundary comparison is distinct from the compact potential
comparison of `D_V` with the constant-potential Dirichlet operator above. The independent endpoint output contributes
at most one more dimension. The resolvent difference therefore has finite
rank, which preserves essential spectrum. It follows that

\[
 \boxed{\sigma_{\rm ess}(A)=[25/4,\infty),\qquad
        \sigma_{\rm ess}(A-4I)=[9/4,\infty).}           \tag{26}
\]

Below `mu^2=9/4`, any spectrum consists of isolated eigenvalues of finite
multiplicity, with possible accumulation only at `9/4`. Bulk limit-point
uniqueness makes each eigenvalue simple. This note does not need a stronger
claim about spectral type at or above the threshold. In particular, (18)
does not rely on an assertion that the spectrum is purely discrete.

The best dimensionless gap cannot exceed `9/4`, because (26) supplies spectrum
at that threshold. Under the archived physical convention
`m_4^2=H_b^2 mu^2`, (16) implies a fixed-background lower bound
`m_4^2>=H_b^2 epsilon_delta`, and the essential threshold is
`(9/4)H_b^2`. Since the static branch has `H_b^2 -> 0` as `delta -> 0`, these
physical essential thresholds approach zero. A positive lower bound on the
entire **physical** mass spectrum independent of detuning is therefore
incompatible with this family of abstract operators; the dimensionless
uniform-gap question is separate and is not settled here.

## 7. Conditional ordinary-Robin companion

For a separately derived scalar sector with no endpoint mass, take
`mathscr H_o=L^2(w dy)` and the ordinary Robin condition

\[
 \alpha(T)u'(T)+c_b u(T)=0.                              \tag{27}
\]

Assume `c_b>=0` and
`r_o=alpha(T)P_0'(T)/P_0(T)+c_b>0`. Use the form

\[
 \mathfrak e_o[u]=\int(\alpha|u'|^2+q|u|^2)+c_b|u(T)|^2
\]

on the bulk version of `V`. The same density, closedness and cone-core
argument yields its self-adjoint Robin realization. Its full ground-state
identity is

\[
 \mathfrak e_o[u]-4\|u\|_{L^2(w)}^2
 =\int a_g|(u/P_0)'|^2+k_o|(u/P_0)(T)|^2,
 \quad k_o=r_oP_0(T)^2>0.
\]

Consequently the entire shifted spectrum lies above

\[
 \varepsilon_o=
       \left[\frac{4m_*}{9\alpha_*}+\frac{M_g}{k_o}\right]^{-1}>0,
 \qquad \sigma_{\rm ess}(A_o-4I)=[9/4,\infty).          \tag{28}
\]

The candidate odd scalar boundary reduction supplied during this continuation
is `c_b=1/(3 a(T) rho(T)^2)`. If that reduction is established independently,
`a(T)>0` and `P_0'=-rho^2 p X_0/3>0` on the registered branch imply all the
Robin signs used here. This note treats that physical boundary reduction as
conditional, not as a result of the operator-domain argument.

## 8. Elementary proof for the sharper Hardy option

Let `R(x)=integral_x^T a_g^-1`, `M(x)=integral_0^x m_g`, and
`H=sup M(x)R(x)`. For `g=f-f(T)`, weighted Cauchy–Schwarz gives

\[
 |g(y)|^2\le2\sqrt{R(y)}
        \int_y^T a_g(s)|f'(s)|^2\sqrt{R(s)}\,ds,
\]

because `integral_y^T [a_g(s)sqrt(R(s))]^-1 ds=2sqrt(R(y))`.
Integration by parts and `M R<=H` give

\[
 \begin{split}
 \int_0^s m_g(y)\sqrt{R(y)}\,dy
 &=M(s)\sqrt{R(s)}+\frac12\int_0^s
             \frac{M(y)}{a_g(y)\sqrt{R(y)}}\,dy\\
 &\le\frac{H}{\sqrt{R(s)}}+\frac H2\int_0^s
             \frac{dy}{a_g(y)R(y)^{3/2}}
 \le\frac{2H}{\sqrt{R(s)}}.
 \end{split}
\]

The omitted lower boundary product is zero since (13) gives
`M(y)sqrt(R(y))=O(y^(3/2))`. Tonelli's theorem now yields
`integral m_g|g|^2<=4H integral a_g|f'|^2`. Apply this first on positive
truncations or to the core, and pass by monotone convergence/closure.
Thus the optional constant in (17) also has a proof internal to this note.

## 9. Standard functional-analysis inputs and scope

The imported standard results are the representation theorem/Friedrichs
extension for a densely defined closed semibounded form, the Weyl
limit-point/limit-circle alternative, the spectral theorem, local Rellich
compactness, and invariance of essential spectrum under compact resolvent
perturbations. Each is invoked with its hypotheses checked above. No external
citation has been invented and no numerical root scan is used as a substitute
for these statements.

This closes the operator/domain and full-spectrum gap step for the stated
scalar differential expression and boundary condition. It does not identify
this Hilbert norm with the physical kinetic form, establish that every
physical mode lies in this realization, resolve exceptional harmonics,
prove positivity of other physical sectors, or establish nonlinear stability.
The Robin companion remains conditional on its independently derived boundary
condition. All statements concern fixed positive detuning unless explicitly
identified as essential-threshold behavior of the family.
