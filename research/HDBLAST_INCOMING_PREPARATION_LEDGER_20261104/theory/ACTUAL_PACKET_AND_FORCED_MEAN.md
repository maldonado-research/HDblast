# Actual incoming-packet test and a forced five-mean ledger

4 November 2026. Conditional technical note. Existing published files are
unchanged. No new field/source oracle is evaluated and no physical preparation
law or full nonlinear solution is asserted.

## 1. The fixed first packet is not time-symmetric or purely incoming now

The normalized master data are exactly

```
w(0,z)=q_ell(z), w_t(0,z)=exp(-z)q_ell(z),
q_ell=exp(4-ell^2/[(z-ell)(2ell-z)]), ell<z<2ell,
q_ell=0 outside, ell=10^-6.
```

These are equation(13) of the retained PHYSICAL_SHELL_RESPONSE input and the
scaled data of COMMON_GAUGE_FORCING/inputs/HIGHER_JETS.md. In x=t/ell,r=z/ell,
`W0=Q(r)`, `W_x0=ell exp(-ell r)Q(r)`. At r=3/2, `Q=1`, `Q_r=0`, so both
principal characteristic derivatives `W_x+W_r` and `W_x-W_r` equal the strictly
positive `ell exp(-3ell/2)`. The first perturbation is therefore neither
time-symmetric nor purely one-way on the present slice. Background time
symmetry `K0=Pi0=0` is a different statement.

For propagation toward the shell r=0, the free one-way condition would be
`W_x-W_r=0`; it would require `Q_r=ell exp(-ell r)Q`. Flatness and `Q(1)=0`
force the zero solution of that first-order ODE, contrary to this Q. The
opposite one-way sign also fails. The actual equation has a potential and a
shell, so characteristic mixing is expected: this calculation does not rule
out a carefully selected incoming past history that backscatters/reflects into
the actual mixed present data. It rules out relabeling the present packet as a
pure incoming pulse without solving that history.

The actual metric/scalar reconstruction is also fixed, not a test-scalar
choice. With p=phi0_y, a=rho_y/rho,

```
alpha=(rho/rho_b)^(3/2)p/p_b,
C=alpha/rho^2, B=3alpha/(rho^3 p), kbg=rho(5a/2-Vbulk'/p),
u=Cq, v=B(q_z+kbg q),
u_t=exp(-z)u,
v_t=exp(-z)[v-3u/(rho p)]=B exp(-z)[q_z+(kbg-1)q].
```

For Y=n1^2-n2^2 the physical first Cauchy data, in radial/angular frame order,
are `gamma1=diag(-4u,2u,2u,2u)Y`,
`K1=diag(-2u_t,u_t,u_t,u_t)Y/rho`, `phi1=vY`, `Pi1=v_tY/rho`.
In the retained action convention `pi^{ij}=sqrt(gamma)(K^{ij}-K gamma^{ij})/2`,
the normalized canonical momentum tangent (orthonormal components) is
`pi1/sqrt(gamma0)=diag(-3u_t/(2rho),0,0,0)Y`,
`varpi1/sqrt(gamma0)=Pi1`. Physical coordinate velocities are
`dot gamma1=2rho K1`, `dot phi1=v_tY`; lapse/shift values are
`N1=rho uY`, `beta1=0`. Harmonic completion fixes rather than freely chooses
`dot N1=rho u_tY`, `dot beta1^y=-6rho^2(u_y+4au)Y`,
`dot beta1^i=0`. These are the retained EXPLICIT_INITIAL_HARMONIC_COMPLETION
formulas(1),(4),(8), and EXPLICIT_COMPACT_SECOND_INITIAL_JET formulas(3),(4).

## 2. Exact inhomogeneous characteristic and shell laws

For the full five mean second coefficient U=(K,L,N,J,w), retain the same
scaled PDE and increasing-r boundary convention:

```
U_xx-U_rr=A U_x+B U_r+Dmat U+f,
PD U=L=0, PN U_r=M u+g, u=PN U(x,0).
```

The actual quadratic adapters give
`f=P0(-2ell^2 rho^2 Qred, ell^2 rho^2 Qphi)`. In mean order the physical
outward-y source rows convert with `partial_r=-ell rho partial_y`:
`g_K=-ell rho GN`, `g_N=2ell rho QR00`,
`g_J=(2/3)ell rho sum_a QRaa`, `g_w=ell rho QP`.
These are one-copy shell Hessian sources with actual-second-variation
normalization. The Dirichlet row has no free extra forcing. No physical
forcing factor1/2 is inserted into these PDE source adapters.

Set `Jplus=U_x+U_r`, `Jminus=U_x-U_r`. Exactly,

```
(partial_x-partial_r)Jplus = A U_x+B U_r+Dmat U+f,
(partial_x+partial_r)Jminus = A U_x+B U_r+Dmat U+f,
PN Jminus=PN Jplus-2(Mu+g),
u_x=PN Jplus-Mu-g,
PD Jminus=-PD Jplus.                                     (1)
```

Thus an incoming physical-shell characteristic combination is constrained by
M and g; it cannot be prescribed independently while keeping the shell law
unchanged. In affine null coordinates a characteristic-cut field trace and
its tangential derivative are linked, and null Einstein constraints impose
additional restrictions on metric/scalar traces. No arbitrary five-component
incoming trace is asserted to be admissible.

At a Cauchy-shell corner, initial data must satisfy
`PD U0=PD V0=0`, `PN (U0)_r=M(0)PN U0+g(0)`.
Differentiating the last law gives
`PN (V0)_r=M_x(0)u0+M(0)PN V0+g_x(0)` whenever required by regularity;
higher corner identities also use the PDE and differentiated sources. A
background-flat initial collar with f,g and their required jets initially zero
is sufficient for the corresponding zero corner rows. It is not automatically
available for a newly chosen past first state.

## 3. A continuous forced energy inequality

Use the same positive analytic energy as the H1 response theorem, with
k=alpha=1/2, V=U_x, Q=U_r, on `0<=r<=T-x`, T=5/4:

```
e=(|V|_H^2+|Q|_H^2+2k<V,Q>_H+|U|_H^2)/2,
E(x)=integral e dr + alpha |u(x)|_H^2/2,
H=diag(1,2,1,3,1).
```

All continuous operator bounds remain the same retained conditional assertions
`||A||,||B||<=10^-4`, `||Dmat||<=10^-6`, `||M||<=1/2000`, not new all-point
coefficient proofs. They cover the existing future cone. No such bounds on a
new past domain or for a different first state are inferred.

The homogeneous bulk rate is at most2.000702; the homogeneous boundary proof
leaves the negative term `-|PN Jplus|^2/8` and boundary rate2.003001. The added
bulk forcing is bounded sharply because

```
2e-|V+kQ|_H^2=(1-k^2)|Q|_H^2+|U|_H^2>=0,
<V+kQ,f>_H <= e+|f|_H^2/2.                                (2)
```

For the shell source, expansion about the homogeneous boundary expression
adds exactly

```
-(1/2)<PN Jplus,g>+<Mu,g>+(1/2)|g|^2-(1/2)<u,g>.
```

Using the reserved negative term and m=1/2000,

```
-|PN Jplus|^2/8+(1/2)|PN Jplus||g| <= |g|^2/2,
(m+1/2)|u||g| <= |u|^2/4+(m+1/2)^2|g|^2.
```

The `|u|^2/4` term is precisely the boundary energy summand. Therefore

```
E_x <= (301/100)E +(1/2)||f(x,.)||_L2(H)^2 +Cg |g(x)|_H^2,
Cg=1+(1/2+1/2000)^2=5002001/4000000.                       (3)
```

No selfadjointness, positivity or time independence of M is assumed. The
moving outer characteristic flux and the physical Dirichlet contribution
remain nonpositive, exactly as in the homogeneous proof. Integrating (3)
gives the actual continuous ledger

```
E(x)<=exp(Cx)E(0)+integral_0^x exp(C(x-s))
          [(1/2)||f(s,.)||_L2(H)^2+Cg|g(s)|_H^2]ds,
C=301/100, exp(CT)=exp(301/80)<44.                          (4)
```

The exponential bound is amenable to the same rational Taylor-tail method;
no new source evaluation is involved in deriving it. The boundary scalar is
bounded by `|w(x,0)|^2<=4E(x)`. All inequalities include the zero-energy case.
This is an analytic component norm, not an identified physical energy law.

For the unchanged future first packet, existing source certificates give
`FH=17387513156573`, `GH=3851400`. On the full cone,
`integral||f||_L2^2 dx <=25 FH^2/32`,
`integral|g|^2 dx <=5 GH^2/4`. Its exact zero wedge improves these to
`9 FH^2/32` and `GH^2/4`, respectively. Those future budgets cannot be used
for an uncomputed incoming past state or past source history. Moreover a
bound on the coordinate mean w alone does not omit the known first-quadratic
proper-clock, measure and Gauss terms in an absolute shell observable.

## 4. Why this does not yet define the retarded physical preparation

Zero second data for the forced mean PDE define an auxiliary solution, not
automatically the second derivative of an Einstein--scalar solution. The
second constraints are `P X2=-D^2Constraints[v,v]`. The retained actual packet
has a previously checked Hamiltonian Hessian witness `<-10^19` at its central
point, so setting all second Cauchy data to zero fails. A past prescription
must supply constraint-compatible free data and forced null/Cauchy fields,
the physical corner state, the shell law and an independently specified past
clock/domain. With zero first perturbation data, no incoming state and no new
physical source, uniqueness cannot create the nonzero fixed packet.

A genuine past problem could select a present second state. It need not lie
in the current compact h-only sector with the other second canonical sectors
fixed. Only after checking that sector condition may one reconstruct the
scalar-level area profile D and apply the H1 class bound. The source-to-shell
ledger(4) is a concrete prerequisite, not completion of that physical history,
constraint compatibility or cosmological-origin hypothesis.


## Retained source links

- [Actual first packet and initial master](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/HIGHER_JETS.md).
- [Physical first Cauchy data and nonzero second constraint witness](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md).
- [Exact harmonic lapse/shift velocities](../../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md).
- [Canonical action/momentum normalization](../../HDBLAST_SHELL_NOETHER_20261013/manuscript/SHELL_NOETHER.tex).
- [Five mean physical boundary rows](../../HDBLAST_PREPARED_SHELL_VOLUME_20261101/gauge/MEAN_INITIAL_GAUGE.md).
- [Conditional continuous operator premises](../../HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/EXPLICIT_PREMISES.json).
- [Retained coupled preparation comparison and proper scalar meaning](../../HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/DERIVATION.md).
- [Earlier direct H1 energy theorem](../../HDBLAST_H1_PREPARATION_RESPONSE_20261104/energy/PROOF.md).
