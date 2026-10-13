# Radial closure and characteristic cone flux: independent audit

13 October 2026 UTC. This note audits the mathematical inference from the
reviewed radial realization to a physical cone boundary condition. It does not
alter the frozen checkpoint or prove completeness of the Einstein–scalar
phase space. All new files are confined to `domain-audit/`.

**Finding.** Limit-point/Friedrichs closure supplies no independently chosen
radial endpoint parameter and a vanishing operator Green form at fixed finite
de Sitter time. It does **not** imply vanishing symplectic or radiation flux
through a finite patch of the characteristic cone. A smooth, finite ambient
energy test-scalar counterexample establishes the failure of that implication
exactly. For the registered gravitational scalar reduction, the independently
derived current has a potentially nonzero characteristic limit under explicitly
stated asymptotic hypotheses; existence and completeness of such registered
solutions remain separate questions.

## 1. Inputs and meanings of the limits

The audited inputs are the checkpoint manuscript,
`domain/SCALAR_SELF_ADJOINT_DOMAIN_AND_GAP.md`,
`operator/PHYSICAL_OPERATOR_AUDIT_20261013.md`,
`tensor/TENSOR_PARTNER_GAP.md`, and `NEXT_CALCULATIONS.md` in
`../basis/prior-checkpoint/`.
The checker records their hashes. Write

\[
 g=dy^2+\rho(y)^2[-dt^2+\cosh^2t\,d\Omega_3^2],\quad
 \rho=y(1+O(y^2)),\quad p=\kappa y(1+O(y^2)),\quad \kappa\ne0.
\]

At fixed finite `t`, `y -> 0` approaches the cone **vertex**. Even in the
exact flat metric `rho=y`, it does not approach an arbitrary null cone point.
The flat ambient coordinates and null coordinates are

\[
 T_{\rm amb}=y\sinh t,\quad R=y\cosh t,\quad
 u=R-T_{\rm amb}=ye^{-t},\quad v=R+T_{\rm amb}=ye^t,
\]
\[
 g=du\,dv+\frac{(u+v)^2}{4}d\Omega_3^2.
\]

The future null cone is `u=0,v>0`. Its limit keeps `v` in a finite positive
interval while `y -> 0` and `t=log(v/y) -> infinity`. Thus a surface
`y=epsilon` samples that patch over a de Sitter time interval translated to
infinity. A finite interval of `t` only samples the vertex.

For the curved regular-cone background use the exact conformal radius
`r=exp integral(dy/rho)`, normalized so `r/y -> 1`. Then
`g=Omega(r^2)^2[dr^2+r^2 ds_dS4^2]` locally when the regular-cone series is
smooth in `y^2`, with `Omega=rho/r`, and set `u=r exp(-t),v=r exp(t)`.
Here `r=y[1-(r_2/2)y^2+O(y^4)]`, so the leading null formulas below agree
with the flat ones. The `y`-based null chart is also smooth under these
even-series hypotheses: writing `f=rho^2/y^2`, the exact metric there is
`g=f(uv) eta+[1-f(uv)](v du+u dv)^2/(4uv)`. Its last coefficient is smooth
because `1-f=O(uv)`. This chart has a finite curvature correction to `g_uu`
on `u=0`; it is not an exactly flat metric there. The conformal radius
is needed for the exact `Omega^2 eta` expression, rather than for smoothness
alone. The null measure `sqrt(-g)g^(uv)=v^3/8` remains unchanged on the cone.

## 2. Scalar Liouville transformation and endpoint form

For `P=rho^2 Psi`, use the checkpoint coefficients

\[
 \alpha=(\rho^2p^2)^{-1},\quad w=(\rho^4p^2)^{-1},\quad
 q=2/(3\rho^2),\quad \tau P=w^{-1}[-(\alpha P')'+qP].
\]

Let `x=integral_y^Tb ds/rho(s)`, `h=(alpha w)^(1/4)`, and `psi=hP`.
Then `dx/dy=-1/rho`, `h=1/(rho^(3/2)|p|)`, and

\[
 \int w|P|^2dy=\int |\psi|^2dx,\qquad
 \int\alpha|P'|^2dy=\int|\psi_x-s\psi|^2dx,
 \quad s=h_x/h.
\]

The exact bulk operator is `-d_x^2+V_L`, where

\[
 V_L=s_x+s^2+\frac23\rho^2p^2,\qquad
 s=\rho\left(\frac32\frac{\rho'}\rho+\frac{p'}p\right).
\]

The differentiated regular-cone hypotheses are necessary to control this
potential from the stated relative remainders. If, more specifically,
`rho=y(1+r_2 y^2+O(y^4))` and
`p=kappa y(1+p_2 y^2+O(y^4))`, with differentiated remainders, then

\[
 \begin{split}
 x&=-\log y+C+(r_2/2)y^2+O(y^4),\\
 h&=|\kappa|^{-1}y^{-5/2}
       [1-(3r_2/2+p_2)y^2+O(y^4)],\\
 s&=5/2+(11r_2/2+2p_2)y^2+O(y^4),\\
 V_L&=25/4+(33r_2/2+6p_2)y^2+O(y^4).
 \end{split}
\]

The `q/w` term begins at order `y^4`; it is not an order `y^2` term.
These coefficients are checked independently by symbolic expansion. At the
level of only the checkpoint hypotheses, the safe conclusion is
`V_L=25/4+O(exp(-2x))` and essential `mu^2=Lambda-4 >= 9/4`.

Define the sesquilinear radial Green form by

\[
 B_y(P,Q)=\alpha(P\overline{Q'}-P'\overline Q).
\]

Its exact transformed expression is

\[
 B_y(P,Q)=-[\psi\overline{\chi_x}-\psi_x\overline\chi],
 \qquad \chi=hQ.                                      \tag{1}
\]

The sign follows from the orientation `x_y<0`. For the reviewed even
realization the form domain is `H^1(0,infinity)` with the shell trace attached
to the positive endpoint component. Its operator domain is `H^2` with that
trace attached; the endpoint output, rather than an imposed fixed Robin
condition, encodes the eigenparameter boundary condition. The odd realization
has the stated regular Robin shell condition. With the bounded potentials
assumed in the checkpoint, all operator-domain `psi,psi_x` tend to zero as
`x -> infinity`. Equation (1) therefore tends to zero for pairs in the
operator domain. This is a radial operator statement.

For reference the leading Frobenius powers are

\[
 P\sim y^{5/2\pm\nu},\quad \nu^2=9/4-\mu^2.
\]

At `mu^2=0` these give `y^4,y`, as in the nodeless-reference argument.
For `mu^2>9/4`, put `nu=ik`: the generalized modes have
`P~y^(5/2 +/- ik)` and `psi~exp(+/- ikx)`. Their radial norms diverge
logarithmically. They are spectral generalized eigenfunctions, not vectors
to which the `H^2` endpoint assertion applies. Suitable spectral packets can
be genuine finite-energy vectors.

The tensor radial Liouville field is instead `psi_T=rho^(3/2) h_T`.
Its Green form `rho^4(h_1 h_2'-h_1' h_2)` also becomes minus the
Liouville Wronskian. Its potential is
`9/4+(45r_2/2)y^2+O(y^4)` for the above series. The different scalar
weight and power must not be exchanged with the tensor ones.

## 3. A form vector need not have a pointwise Green-form trace

`H^1` implies that `psi` is continuous and tends to zero; it does not imply
that `psi_x` has a boundary value. This matters when "finite energy" means
the closed form domain alone. Here is a constructive smooth counterexample.
Choose disjoint centers `x_n=3n`, numbers `a_n=2^-n`, and compact smooth
bumps `chi,eta` supported in `[-1,1]` satisfying
`chi(0)=1,chi'(0)=0,eta(0)=0,eta'(0)=1`. Set

\[
 \psi_1=\sum_n a_n\chi(x-x_n),\qquad
 \psi_2=\sum_n a_n^3\eta((x-x_n)/a_n^4).
\]

Both are smooth on each finite interval and in `H^1`:
`||psi_1||_H1^2` is bounded by a constant times `sum a_n^2`, while
`||psi_2||_L2^2` and `||psi_2'||_L2^2` are bounded respectively by
constants times `sum a_n^10` and `sum a_n^2`. Yet their Wronskian is `1`
at every `x_n` and zero between supports. Thus no pointwise endpoint
Wronskian limit exists. The closed variational identity remains valid;
it does not secretly assert a pointwise boundary product for all form
vectors. A weak finite-energy flux must be defined by conservation or
approximation, with the chosen topology stated.

## 4. What a finite-time vanishing statement can prove

For a compact de Sitter time interval `I`, a sufficient strong radial
condition is continuity into `H^2_x(L^2(S^3))` for each master field. The
Wronskian then tends to zero **uniformly on I** after angular integration.
Indeed compact subsets of `H^1` have uniformly small tails; the one-dimensional
Sobolev estimate applied to those tails makes the point values uniformly
small. Apply this to both `psi` and `psi_x`.

The physical gravitational scalar current involves more than the bare
Wronskian for mixed masses. The covariant-current calculation in the sibling
derivation gives, for one interior copy on the generic longitudinal solution
slice and a fixed orientation,

\[
 j_y=9\sqrt{-\gamma}\,\alpha\{[(\Box_\gamma+4)P_1]P_2'
                       -[(\Box_\gamma+4)P_2]P_1'\}.    \tag{2}
\]

For equal separated `Lambda`, (2) reduces to `9Lambda sqrt(-gamma)`
times the real radial Wronskian. In general it does not. Set
`zeta_i=h(Box+4)P_i`. Equation (2) transforms to

\[
 -9\sqrt{-\gamma}\,[\zeta_1(\psi_{2x}-s\psi_2)
                   -\zeta_2(\psi_{1x}-s\psi_1)].        \tag{3}
\]

A sufficient finite-time condition for (3) is
`psi_i in C(I;H^2_x(L^2Omega))` and
`zeta_i in C(I;H^1_x(L^2Omega))`, with bounded `s` and enough angular
regularity for the fields and current to be defined. Then its angular integral
tends to zero uniformly on `I`; integrating over a fixed time interval
also tends to zero since `cosh^3 t` is bounded there. Smooth solutions
continuous into a sufficiently high graph domain of the radial operator
are one way to satisfy these conditions. The bare form domain alone is
insufficient for this pointwise proof. Equation (2) is independently derived
by the covariant agent; this audit checks its boundary implications, not
the complete Einstein–scalar reduction or gauge quotient.

None of these bounds is uniform for time intervals translated to infinity.
An abstract self-adjoint wave evolution can conserve its finite-time
symplectic form while having nontrivial scattering data at the infinite
endpoint over infinite time. Essential spectrum is a warning that such
propagation channels exist, not proof that they are absent physically.

## 5. Exact smooth finite-energy counterexample in the flat cone

Consider a minimally coupled massless test scalar in five-dimensional
Minkowski space, restricted to the exterior cone wedge. Its radial Hilbert
weight is `y^2`, radial form is `integral y^4 |f'|^2dy`, and radial
operator is `-y^-2(y^4 f')'`. At zero spectral value the two radial
solutions are `1,y^-3`; only the first is in `L^2(y^2dy)` near zero.
Thus the endpoint is limit-point, with the finite-energy/Friedrichs
realization and Liouville potential `9/4`.

The exact local solutions

\[
 \varphi_1=1,\qquad \varphi_2=T_{\rm amb}=y\sinh t
\]

are smooth at the vertex and at every null cone point. At every finite `t`
their radial Liouville profiles `y^(3/2)` and
`y^(5/2) sinh t` lie in `H^2` near the radial endpoint. They therefore
obey the fixed-time radial endpoint condition.

These local solutions are realized by actual smooth globally finite-energy
solutions: choose compact smooth initial data on `T_amb=0` that agree with
`(varphi_1,dot varphi_1)=(1,0)` and
`(varphi_2,dot varphi_2)=(0,1)` on a ball of radius `R_0`. By finite
propagation speed they equal the displayed solutions wherever
`|T_amb|+R<R_0`. Take a cone segment `0<a<=v<=b<R_0`; there
`T_amb=R=v/2` and the equality holds on a neighborhood. The initial
energies are finite since all data are compact smooth. If a remote timelike
shell is added outside this causal diamond, it does not affect this local
counterexample. This is a test-scalar comparator, not an assertion that
these two fields solve the registered coupled Einstein–scalar system.

Using `omega^a=varphi_1 grad^a varphi_2-varphi_2 grad^a varphi_1`, the
coordinate-density flux on `y=epsilon` is

\[
 j_y\,dt=\epsilon^4\cosh^3t\sinh t\,dt\,d\Omega_3.
\]

On any fixed interval `[t_1,t_2]`, the integrated flux per solid angle is
`epsilon^4[cosh^4(t_2)-cosh^4(t_1)]/4 -> 0`. On the interval
`t in [log(a/epsilon),log(b/epsilon)]`, however,

\[
 j_y\,dt=
 \frac{(v+\epsilon^2/v)^3(v-\epsilon^2/v)}{16v}
 \,dv\,d\Omega_3
 \longrightarrow \frac{v^3}{16}\,dv\,d\Omega_3.
\]

The limiting integrated flux is `(b^4-a^4)/64` per solid angle, nonzero.
It also follows directly from the null current:
`sqrt(-g) g^(uv)=R^3`, `partial_v T_amb=1/2`, and `R=v/2` on `u=0`.
The same local solution has nonzero null derivative and hence radiative
stress-energy flux. The limit is uniform on `[a,b]` because `a>0`.

This exact example settles the logical question: radial limit-point closure,
finite radial energy at each finite de Sitter time, smooth cone extension,
and even finite ambient energy do not jointly imply zero characteristic
symplectic flux. They are compatible with flux into the region beyond the
wedge. It is not a contradiction of the radial self-adjoint theorem.

### 5.1 The radial and inertial energy topologies are distinct

There is also a separate exact initial-data comparator. At `t=0` in the
flat cone, the spatial volume is `y^3 dy dOmega_3` and the lapse of the
de Sitter time flow is `y`, since
`partial_t=y partial_Tamb` on that slice. For a test scalar with zero
initial velocity, its instantaneous energy for this time flow has radial
gradient part proportional to `integral y^4 |f_y|^2dy`. Its inertial
Minkowski energy instead has gradient part
`integral y^3 |f_y|^2dy`. The former is a lapse weighted energy. The
radial Hilbert norm is `integral y^2 |f|^2dy`.

Take an outer smooth cutoff equal to one near the vertex, and let
`f(y)=cutoff(y) y^(-a)`, `1<=a<3/2`, with zero initial velocity.
Only the powers near zero affect convergence. There,

\[
 y^2|f|^2=y^{2-2a},\qquad
 y^4|f_y|^2=a^2y^{2-2a},\qquad
 y^3|f_y|^2=a^2y^{1-2a}.
\]

The radial norm and closed form are finite because `2-2a>-1`.
The inertial gradient energy is infinite because `1-2a<=-1`
(logarithmically at `a=1`). For `a=5/4`, the exponents are respectively
`-1/2,-1/2,-3/2`. The field is not a smooth classical field at the vertex;
it is a member of the radial form completion. Cutting it off below
`epsilon` approximates it in radial norm plus form, at squared cost
`O(epsilon^(3-2a))`, including the cutoff derivative term. The same
approximants do not converge in inertial gradient energy: that energy
diverges as the cutoff is removed. This is a strict distinction of the
two natural initial-data topologies even in the simplest exact geometry.

This example disproves identifying those energies without specifying the
lapse and reconstruction. It does **not** prove that the registered
Einstein–scalar radial completion differs from a chosen physical gravity
completion: its constraints, physical kinetic bilinear, gauge quotient,
shell conditions and reconstruction must be analyzed separately. The
smooth compact-data flux example above has finite energies in both senses
and therefore does not rely on this topology mismatch.

## 6. Registered scalar: conditional characteristic asymptotics

The exact global-de-Sitter master equation is

\[
 \psi_{tt}+3\tanh t\,\psi_t-\operatorname{sech}^2t\,
 \Delta_{S^3}\psi+[-\partial_x^2+V_L-4]\psi=0.
\]

Writing `Phi=(cosh t)^(3/2) psi` gives exactly

\[
 \Phi_{tt}-\Phi_{xx}+
 [V_L-25/4+(3/4-\Delta_{S^3})\operatorname{sech}^2t]\Phi=0. \tag{4}
\]

The asymptotic `9/4` scalar mass threshold cancels the de Sitter damping
term in this normalization. Equation (4) approaches a massless radial wave
equation along `x,t -> infinity`, `t-x` fixed. A profile `F(t-x)` is
therefore an asymptotic radiation ansatz, **not an exact solution** of (4)
unless its remaining potential terms are controlled or vanish. No scattering
theorem for the registered background is silently invoked here.

Suppose actual constrained solutions satisfy, uniformly on a compact positive
`v` interval and with the differentiated remainders needed in (2),

\[
 P_i=y^4 A_i(v,\Omega)+O(y^6),\quad
 P_i'=y^3 D A_i+O(y^5),\quad
 (\Box+4)P_i=y^4 L A_i+O(y^6),
\]
\[
 D=4+v\partial_v,\qquad L=4-4v\partial_v-v^2\partial_v^2.
\]

One may use `r` for the exact conformal-flat background formula; replacing
`r` by `y` changes these leading formulas only by the controlled subleading
terms. For angular
dependence `sech^2t Delta P` is subleading of order `y^6` along this limit.
The reconstruction then has

\[
 \Psi_i=y^2 A_i+O(y^4),\qquad
 X_i=-\frac3\kappa D A_i+O(y^2).
\]

These leading fields are compatible with smooth adapted null coordinates.
In the flat leading metric,

\[
 \delta g_i=2uv A_i\eta-\frac32 A_i(v\,du+u\,dv)^2,
\]

which is smooth at `u=0,v>0` for smooth `A_i`. This compatibility does not
prove a smooth registered solution, a smooth global longitudinal gauge, or
its extension through the vertex.

For the **actual current** (2), the assumed expansions give one-copy flux

\[
 \begin{split}
 j_y dt&\longrightarrow
 \frac{9v^2}{8\kappa^2}[(LA_1)DA_2-(LA_2)DA_1]\,dv\,d\Omega_3\\
 &=\frac{9v^3}{8\kappa^2}
       [DA_1\,\partial_v(DA_2)-DA_2\,\partial_v(DA_1)]\,dv\,d\Omega_3\\
 &=\frac{v^3}{8}[X_1\partial_vX_2-X_2\partial_vX_1]\,dv\,d\Omega_3.
                                                               \tag{5}
 \end{split}
\]

The key exact identity is `L A=D A-v partial_v(D A)`.
The remainder in (5) is `O(y^2)` on `[a,b]` under the stated differentiated
asymptotic bounds, so integrating the limit is justified there. Identical
doubled copies contribute twice this current when the two-copy convention
applies. Reversing the chosen hypersurface orientation reverses the sign.

Arbitrary smooth `A_i` need not produce zero (5). For example the local
leading data `A_1=-kappa/12`, `A_2=-kappa v/30` give `X_1=1,X_2=v/2`
and yield `v^3/16` per solid angle. These values check the algebra against
the exact test-scalar comparator. They are **not** presented as constructed
global registered Einstein–scalar perturbations.

Thus the correct next theorem would establish a characteristic solution
space, prove the needed remainder bounds and gauge regularity, then compute
or constrain (5). The current and radial closure alone cannot provide that
existence theorem. Conversely a "zero cone flux" phase space can be chosen,
but its radiation restriction must be stated independently of the radial
limit-point classification and may exclude smooth radiative extensions.

## 7. Narrow conclusions and remaining work

1. The checkpoint's radial self-adjoint domains, essential thresholds and
   positive scalar bounds are not invalidated by this audit.
2. Fixed-time strong-domain Green form vanishing is rigorous. A pointwise
   assertion for arbitrary form vectors is false without added regularity.
3. Whole-cone flux is a different limit. The exact smooth flat-cone example
   demonstrates nonzero flux with finite-energy data and radial closure.
4. The registered generic scalar current has the conditional null limit (5).
   Proving registered radiative solution existence, the admissible gauge
   quotient, shell matching, vertex behavior, and the relation of the radial
   norm to a chosen global energy remains necessary.
5. No statement about global physical completeness, cosmological origin,
   nonlinear stability, or dynamics between static branches follows here.

The companion `check_domain_flux.py` provides exact algebra, deliberate
wrong-formula controls, source hashes, and numerical values exhibiting the
two limits. Its successful run validates those identities; it is not a
replacement for the functional-analytic or hyperbolic existence arguments.
