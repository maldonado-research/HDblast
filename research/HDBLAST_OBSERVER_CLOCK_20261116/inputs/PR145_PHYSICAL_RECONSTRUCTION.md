# Exact physical reconstruction and finite-regularity shell join

16 November 2026. New calculation at base `6ee27c36dac67497a2454f7c1283fd7177a91888`. This proves an exact conjugacy and boundary-trace join for the already obtained conditional reduced solution. It does not assume that an earlier generic evolution theorem constructs the current solution. It proves the physical **reduced** equations and the physical primal shell laws; vanishing of the relative harmonic vector is the separate initial-data/subsidiary-energy join.

## 1. Pinned inputs and statement

The source identifiers P1--P9 below refer to exact paths, byte counts and SHA256 values in `PROOF_SOURCE_PINS.json`. All source files are read-only.

| Pin | Source and role |
|---|---|
| P1 | `HDBLAST_NONLINEAR_PICARD_LIMIT_20261114/REDUCED_EVOLUTION_AND_CAUSAL_SUPPORT.md`: the actual conditional sixteen-field solution, original data and causal support |
| P2 | Same capsule, `TOP_ORDER_TIME_CONTINUITY.md`: individual strong continuity of every X14 jet |
| P3 | `HDBLAST_EVOLUTION_BOUNDS_20261106/COEFFICIENTS.md`: exact physical point map and all sixteen coordinate rows, including the varying reference connection |
| P4 | `HDBLAST_BOUNDARY_FLOW_PATCH_20261106/BOUNDARY_FLOW_PROOF.md`: actual proper-y, background-subtracted state flow and its inverse |
| P5 | Same capsule, `COORDINATE_DICTIONARY.md`: proper y versus conformal z, normal orientation and fixed shell normalization |
| P6 | `HDBLAST_TEMPORAL_JETS_SEED_20261114/MIXED_INITIAL_JET_PRODUCTS.md`: exact coordinate-first-jet to fixed-product normalization |
| P7 | `HDBLAST_UNIFORM_CAUCHY_DATA_20261113/UNIFORM_BOUNDARY_CAUCHY_TRANSPORT.md`: the nonzero-state initial chart and its exact time velocity |
| P8 | Same capsule, `UNIFORM_CAUCHY_DATA.md`: current half-action canonical normalization and completed raw physical Cauchy pair |
| P9 | `HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md`: definitions of the physical residuals and relative harmonic vector; its older canonical density convention is not substituted for P8 |

Let `ell=10^-6`, `T=13ell/8`, and `Sigma=[0,6ell] x S3`. P1--P2 provide the actual elapsed-past solution `U` with

```
partial_s^j U in C([0,T];H^(14-j)),       0<=j<=14,
sup_s X14(U;s)<=M<R/4,                   R=10^-12,
U_D=0, D_z U_N=0 at z=0;                U=0 at z=6ell.
```

There are four Dirichlet fields and twelve Neumann fields. The derivative tensors and their slots use the fixed positive product metric and full Frobenius norms. Define the original-time straightened field by

```
V(t)=J_E U(-t),       -T<=t<=0.                           (1)
```

The construction below yields a Lorentzian metric `g` and scalar `phi` on this original past slab such that

```
partial_t^j(g,phi) in C([-T,0];H^(14-j)),    0<=j<=14,  (2)
partial_t^j mathcal C in C([-T,0];H^(13-j)), 0<=j<=13,  (3)
mathcal C^A=g^BC(Gamma[g]^A_BC-Gamma[g0]^A_BC),
mathscr E_AB=Ric_AB-phi_A phi_B-(2/3)mathcal U(phi)g_AB
            =nabla_(A mathcal C_B),
mathscr S=Box_g phi-mathcal U'(phi)=0.                   (4)
```

Parentheses include the factor one half. `mathcal U` denotes the fixed scalar potential, avoiding confusion with the elapsed unknown `U`; the coefficient named C in P1/P6 is denoted `Ccoef` below. At the shell, with the physical outward normal,

```
g_ya=0,  K_ab=sigma(phi)h_ab/6,
n_out(phi)=-sigma'(phi)/2,  mathcal C^y=0.               (5)
```

These are ten Israel rows, one scalar row, one normal harmonic row and four coordinate rows. They are asserted as the actual finite-regularity traces of the reconstructed solution. No four additional tangential harmonic rows are imposed.

## 2. The literal point chart, its image and its tensor reflection

In a round orthonormal angular coframe `theta`, use the tangential de Sitter coframe `(dt,cosh(t)theta)` and put `J4=diag(-1,1,1,1)`. The background is

```
g0=rho(z)^2[dz^2-dt^2+cosh(t)^2 kappa],
phi0=1+eta0(z),         dy=-rho(z) dz.
```

This y is the fixed background coordinate, not perturbed proper distance. Write

```
u=F_z(V)=(v,Htilde,w,D),       N=1+v.
```

P4 supplies the flow of its displayed shell-frozen generator, on the whole audited collar where `theta(z)=z`, with `F_-z F_z=Id` on the actual image and `F_z(0)=0`. Reconstruction in the frame `(dz,dt,cosh(t)theta)` is exactly

```
G_zz=rho(z)^2 N^2,
G_za=-rho(z)rho_b D_a,
G_ab=rho(z)^2 J4_ab+rho_b^2 Htilde_ab,
phi=phi0(z)+w.                                         (6)
```

Here every `rho_b=rho(0)` factor is constant; no moving normalization replaces it. `N=sqrt(g_yy)` is the radial coordinate component root. It is neither the temporal ADM lapse nor, away from `D=0`, a radial Schur-complement lapse.

In a stereographic chart, `sph(xi)=2/(1+|xi|^2)` and

```
E=diag(1,1,sph cosh t,sph cosh t,sph cosh t),
g=E^T G E,                 X=(g,phi)=Psi(x,V).           (7)
```

Every base derivative of Psi differentiates `rho,phi0,E,F_z`, as prescribed in P3. The maps are equivariant under the orthogonal changes of round angular frame: the generator uses inverse, multiplication and trace against J4, and the reconstruction is tensorial. Thus the local descriptions glue on S3, including orientation-reversing angular reflections.

The current state lies strictly inside the certified image patch. Indeed the accepted tensor embedding gives `|V|star<=2*10^7 R/4=5*10^-6<10^-3`; P4 then bounds `|F_z(V)|star<.00101` and proves inverse-chart membership. P3/P4 give `||rho^-2 G-J5||op<1/40` and the corresponding Lorentzian/nondegenerate metric margins. The Jacobian `J=Psi_V` is invertible: its factors are the invertible flow derivative, the radial square map with derivative `2N`, the nonzero fixed scales and the coframe conversion. This is invertibility on the actual chart image, not an assertion that a ball maps onto itself.

For time reflection put `S_t=diag(-1,1,1,1)` on tangential slots. Explicitly

```
J_E(v,Htilde,w,D)=(v,S_t Htilde S_t,w,S_t D).
```

Because `S_t J4 S_t=J4`, the generator satisfies `B(J_E u)=J_E B(u)`. Uniqueness of its finite-dimensional flow gives `F_z J_E=J_E F_z`. Also `cosh(-s)=cosh(s)`, and the background is invariant under the coordinate reflection `r(s,z,omega)=(-s,z,omega)`. Consequently the reconstruction from the elapsed field is exactly the tensor pullback of the original reconstruction:

```
Psi_s(s,U(s))=r^*[Psi_t(t,V(t))],        t=-s.           (8)
```

In particular `U_s=-J_E V_t`, `U_ss=J_E V_tt`; the inverse metric has `q_s^ss=q_t^tt`, `q_s^si=-q_t^ti`, and unchanged spatial entries. Its normalized mixed principal coefficient therefore changes sign as required. The lower source transforms with its output, time derivative and tensor slots. The round product connection commutes with J_E, so this includes its first-jet correction. Equations (6)--(8) establish the same-state pullback, not a replacement by another parameter state or a future solution.

## 3. Undo the full normalized fixed-product equation

It suffices by (8) to work in original t. Use the exact normalized maps of P6, with `a=-q^tt>0`, where `q=g^-1` is the metric reconstructed in (6)--(7). They are

```
A^i=2q^ti/a,       B^ij=q^ij/a,
Ccoef=-cbar/a,     Ecoef=-bbar/a,     Gcoef=-dbar/a.
```

Thus the actual equation for V, including the minus sign on the **complete** lower source, is

```
V_tt-(2q^ti/a)D_i V_t-(q^ij/a)Hsp_ij V
     =-(cbar+bbar DV+dbar(DV,DV))/a.                    (9)
```

Multiplication by `-a`, retaining both symmetric time/spatial mixed entries and `H_ti V=D_i V_t`, gives

```
q:Hprod V=cbar+bbar Q+dbar(Q,Q),          Q=DV.          (10)
```

This operation is legitimate everywhere in the certified patch; P6 uses the actual all-state lower bound `a>=1/18000`. If (9) is written in the fixed spatial orthonormal frame, its contractions are converted tensorially to the coordinates in (7) before the next step.

In such coordinates let `D_i=partial_i+Omega_i`, and let `gamma^k_ij` be the product base connection. Define `K(V)_i=Omega_i V`, so the ordinary first jet is

```
P_i=partial_i V=Q_i-K(V)_i.
```

The exact P6 coefficient identities are

```
cbar=c_c-b_c K+d_c(K,K)+T V,
bbar Q=b_c Q-2d_c(K,Q)+S Q,        dbar=d_c,
S Q=q^ij(Omega_i Q_j+Omega_j Q_i-gamma^k_ij Q_k),
T=q^ij(partial_i Omega_j-Omega_j Omega_i).              (11)
```

They include the full substitution in the old coordinate source, not only a change in its first-order term. The product Hessian expands as

```
Hprod_ij V = partial_ij V+Omega_i P_j+Omega_j P_i
             -gamma^k_ij P_k
             +(partial_i Omega_j+Omega_i Omega_j
                              -gamma^k_ij Omega_k)V.
```

Substituting `P=Q-Omega V` cancels the `gamma Omega` terms and leaves exactly the `S Q+T V` in (11), including the ordered product `-Omega_j Omega_i`. Therefore (10)--(11) are equivalent to

```
q^ij partial_ij V=f_c(x,V,P),
f_c=c_c+b_c P+d_c(P,P).                                 (12)
```

No product-frame correction, direction-slot action or state-connection term is omitted in this reversal. The horizontal derivative convention of P6 is precisely the convention that makes these identities covariant along an arbitrary nonzero state.

## 4. Undo the point map and identify the literal physical residuals

P3 defines the actual coordinate source in (12) by the identity

```
f_c = J^-1 { Qphys(x,Psi,Psi_x+J P)
       -q^ij[Psi_ij+2Psi_iV P_j+Psi_VV(P_i,P_j)] },
J=Psi_V.                                               (13)
```

Here `Qphys` is the physical sixteen-row source below; it is distinct from the first-jet array Q in (10). With `X=Psi(x,V)`, the full point-map chain rule is

```
partial_ij X = J partial_ij V+Psi_ij
                 +Psi_iV P_j+Psi_jV P_i+Psi_VV(P_i,P_j).
```

Contracting with symmetric q, substituting (12)--(13), and cancelling the three displayed chain terms gives the exact coordinate equations

```
q^ij partial_ij X=Qphys(x,X,partial X).                 (14)
```

This cancellation differentiates the moving coframe, flow and entire background in Psi; it is more than matching scalar principal symbols.

To specify the physical source and its reference terms explicitly, write `k_iab=partial_i g_ab`, `p_i=partial_i phi`, `Kbar^a_bc=Gamma[g0]^a_bc`, and define

```
Gamma^a_bc=(1/2)q^ar(k_brc+k_crb-k_rbc),
q_i^ab=-q^ar k_irs q^sb,
T_i^a_bc=-q^ar k_irs Gamma^s_bc,
mathcal C^a=q^bc(Gamma^a_bc-Kbar^a_bc),
C0_i^a=q_i^bc(Gamma^a_bc-Kbar^a_bc)
                    +q^bc(T_i^a_bc-partial_i Kbar^a_bc),
mathcal C_b=g_ba mathcal C^a,
N0_ib=k_iba mathcal C^a+g_ba C0_i^a-Gamma^a_ib mathcal C_a,
R0_ab=T_c^c_ab-T_b^c_ac
                +Gamma^c_cr Gamma^r_ab-Gamma^c_br Gamma^r_ac,
L_ab=R0_ab-(N0_ab+N0_ba)/2-p_a p_b
                                  -(2/3)mathcal U(phi)g_ab,
Qphys_ab=2L_ab,
Qphys_phi=q^ij Gamma^a_ij p_a+mathcal U'(phi).            (15)
```

`T,C0_i,N0` in this display mean the expressions obtained by setting only the metric's second jet to zero; the reference derivatives `partial_i Kbar` remain. The scalar potential is exactly the retained one, `mathcal U=(W')^2/2-2W^2/3`, with `W=1/3+eta^2+eta^3/3`, `eta=phi-1`. The actual reference `g0=rho(z)^2(dz^2-dt^2+cosh(t)^2 kappa)` is position dependent. Neither its Christoffel symbols nor their derivatives are frozen in (15). Only the auxiliary boundary generator freezes shell coefficients.

Here is the general second-jet cancellation behind (15). The second-derivative part of Ric_ab is

```
(1/2)q^cd(partial_ac g_bd+partial_bc g_ad
                         -partial_cd g_ab-partial_ab g_cd),
```

whereas the second-derivative part of `nabla_(a mathcal C_b)` is

```
(1/2)q^cd(partial_ac g_bd+partial_bc g_ad-partial_ab g_cd).
```

Their difference is `-(1/2)q^cd partial_cd g_ab`. Every remaining term is exactly L_ab from (15). Thus the metric part of (14) implies

```
mathscr E_ab-nabla_(a mathcal C_b)
              =-(1/2)q^ij partial_ij g_ab+L_ab=0.
```

The scalar part gives `q^ij(partial_ij phi-Gamma^a_ij partial_a phi)-mathcal U'=0`. This proves (4) as an identity for the actual reconstructed fields. It is an arbitrary-jet algebraic identity, not a conclusion from a finite sample test or a proposed generic theorem.

## 5. Reconstruction and harmonic-vector regularity

The regularity needed in the preceding identities follows directly from the already proved strong X14 table. P4 supplies mixed flow derivatives through eighteen; P3 bounds Psi derivatives through eighteen and reference Christoffel derivatives through seventeen. On the compact actual state patch the required maps, including metric inversion, are smooth with finite uniformly bounded derivatives. The inverse metric has the independently established nondegeneracy margin; no new smallness condition is imposed here.

For each `j+k<=14`, differentiate `Psi(x,V(t))` j times in t and k times covariantly in space, keeping every slot. The result is a finite sum of an explicit-base term and coefficient arrays applied to mixed derivatives of V. The actual field factors have total derivative order at most fourteen. Select a largest-order factor for L2; every other factor has order at most seven and three spare spatial derivatives in the same actual jet's X14 norm. It is consequently in L-infinity, with maximum embedded rank `2+7+3=12`. The high factor uses its own spatial Sobolev regularity; unrelated temporal jets are not interpolated. Coefficient arrays are bounded directly and are not Sobolev-embedded. Finite telescoping of these products proves strong time continuity as well as boundedness. The explicit-base term is smooth. This proves (2), with all components in fixed background coordinate/bundle identifications.

The same argument applies to `q=g^-1` on the nondegenerate image and to the coordinate change y(z). The relative harmonic vector and covector have the form of a smooth coefficient in g multiplying one actual first derivative of g, plus a smooth coefficient multiplying the full fixed reference Kbar. A mixed derivative of total order at most thirteen therefore uses field derivatives of total order at most fourteen. The same allocation proves (3) for both upper and lower components. No derivative of the background reference is missing from this estimate.

In particular each spacetime derivative of `(g,phi)` of total order at most eleven is continuous: its spatial Sobolev remainder is at least H3 in dimension four. Thus `(g,phi)` is piecewise C11 up to the boundary, and `mathcal C` is piecewise C10. This is finite regularity, not spatial C-infinity. It is more than sufficient for the classical one-sided Bianchi and Codazzi identities and all pointwise products used below. The reduced residuals in (4) belong to `C H12`; in particular their initial values are strong H12 identities. Ordinary traces of first derivatives belong to `C H^(25/2)` on the shell, while values have `C H^(27/2)` traces. No ordinary trace is assigned to the L2 fourteenth time jet.

## 6. Exact shell traces from the fixed D/N rows

At z=0 the flow has `F_0=Id`, `(F_0)_V=Id`, and `partial_z F_z|0=B`. The chain rule and the actual field traces therefore give

```
u_D=0,                    u_N,z=B_N(u_N).               (16)
```

Conversely (16) gives the original fixed D/N rows by the inverse chain rule. The chain rule holds for these Sobolev traces, or pointwise using the C11 regularity just proved. It does not require that D vanish off the shell.

Put `H=J4+Htilde`, `sigma_b=sigma(phi0(0))`, and `a_b=rho_y(0)/rho_b`. The retained background relations are

```
sigma_b=6a_b,          p_b=phi0_y(0)=-sigma'_b/2,
rho_z(0)=-a_b rho_b^2.
```

The actual generator P4, with constant shell scalings, is

```
B_v=-rho_b{(2/3)[sigma(phi)N^2-sigma_b]
                       -a_b[N^3 tr(H^-1 J4)-4]},
B_H=-rho_b[sigma(phi)N H-sigma_b J4]/3,
B_w= rho_b[N sigma'(phi)-sigma'_b]/2,
B_D=0.                                                  (17)
```

Its tension is the retained exact polynomial in `eta=phi-1`, with the stated fixed delta and c; no expansion in the perturbation is used. The time-dependent tangential coframe is independent of z, so normal differentiation of its frame entries is the same as in (6).

First, (6) and `D|0=0` give `g_ya=rho_b D_a=0` and `g_za=0`, for all four tangential directions including time. Hence also `g^ya=g^za=0` on the shell, and

```
n_out=N^-1 partial_y=-(rho_b N)^-1 partial_z.             (18)
```

Second, the complete normal derivative of the tangential physical metric is

```
h_z=-2a_b rho_b^3 J4+rho_b^2 B_H
    =-2a_b rho_b^3 J4-(rho_b^3/3)
                                 [sigma N H-sigma_b J4]
    =-(rho_b sigma N/3)h.                               (19)
```

The cancellation uses `sigma_b=6a_b`; discarding the differentiated background term would give the wrong row. Since `partial_y=-rho_b^-1 partial_z` at the shell,

```
h_y=(sigma N/3)h,       K_ab=(1/(2N))h_ab,y=sigma h_ab/6. (20)
```

The formula for K follows from `K=(1/2)Lie_n h`; tangential derivatives of `N^-1` multiply the vanishing mixed trace and contribute zero. It does not assume a Gaussian gauge away from the boundary.

Third, the scalar subtraction gives

```
phi_z=-rho_b p_b+B_w=rho_b N sigma'(phi)/2,
n_out(phi)=N^-1 phi_y=-sigma'(phi)/2.                    (21)
```

No division by sigma' occurs, so zeros of sigma' create no missing branch.

Finally the v row, with its exact constant cancellation, yields

```
N_y=(2/3)sigma N^2-a_b N^3 tr(H^-1 J4)
    =(2/3)sigma N^2-a_b N^3 tr_h h0,b.                  (22)
```

Compute the relative harmonic normal component in proper-y coordinates. Tangential derivatives of the trace `g_ya=0` vanish, so

```
g^BC Gamma[g]^y_BC=N_y/N^3-(1/(2N^2))h^ab h_ab,y
                  =N_y/N^3-2sigma/(3N),
g^BC Gamma[g0]^y_BC=-a_b tr_h h0,b.
```

Here the reference has `GammaBar^y_yy=0` and `GammaBar^y_ab=-a_b h0,b,ab`; its cross contraction vanishes. Substituting (22) gives `mathcal C^y=0` exactly. Since the mixed trace is zero,

```
mathcal C_y=N^2 mathcal C^y=0,
mathcal C(n_out)=mathcal C_y/N=0,
mathcal C_z=-rho_b mathcal C_y=0.                        (23)
```

This proper-y calculation is consistent with, but does not erase, the nonzero conformal reference entry `GammaBar^z_zz=-a_b rho_b`. Assigning the vanishing proper-y radial Christoffel to z would alter (22). The full bulk reduction in sections 3--4 remains in its actual conformal coordinate, with its actual reference derivatives.

Equations (19)--(23) also prove the converse: given `g_ya=0`, the Israel and scalar laws and `mathcal C^y=0`, the raw perturbations satisfy exactly (17), hence (16), hence the fixed D/N rows. This closes the nonlinear primal boundary equivalence at the finite regularity of the current solution.

## 7. Exact return to the current physical Cauchy pair

This subsection identifies the data, not merely their tangent at zero. Let `(X_raw,Xdot_raw)` be P8's completed background-subtracted physical Cauchy pair. P8 uses the current half-action density conversion

```
K_initial=2lambda[gamma tau gamma-(tr_gamma tau)gamma/3],
Pi=lambda omega,               lambda=(det gamma/det gamma0)^-1/2,
alpha=rho(1+epsilon f),        beta=0.
```

Its full completion includes `gamma_t=2alpha K_initial`, `phi_t=alpha Pi`, the actual harmonic lapse/shift time derivatives, and generally nonzero `g_ti,t`. These current data are not replaced by the doubled-density formula in the older P9 note. The initial-data companion proves their exact physical constraints and initial harmonic identities.

Let I and L be precisely the fixed fiber identification and constant-shell normalization of P7. Set

```
P=L I X_raw,              Qraw=L I Xdot_raw,
r(P)=(sqrt(1+x)-1,Htilde,w,D),
T(z,P)=F_-z(r(P)).
```

Explicitly, in the raw coframe `(dt,dy,rho theta)` at t=0, this normalization is

| State entry | Raw completed perturbation entry |
|---|---|
| `x` | `delta g_yy` |
| `Htilde_tt` | `delta g_tt/rho_b^2` |
| `Htilde_tA` | `rho delta g_(t,phys A)/rho_b^2` |
| `Htilde_AB` | `rho^2 delta g_(phys A,phys B)/rho_b^2` |
| `w` | `delta phi` |
| `D_t` | `delta g_yt/rho_b` |
| `D_A` | `rho delta g_(y,phys A)/rho_b` |

I is a fixed fiber identification, not the coordinate substitution `dy=-rho dz`. The latter is applied separately in (6), producing `g_za=-rho rho_b D_a` and `g_zz=rho^2 N^2`. Thus the radial sign/factor and the angular rho powers in this table are both retained when checking the inverse.

P7 gives the exact nonzero-state Cauchy transformation

```
V(0)=T(z,P),                V_t(0)=T_P(z,P)[Qraw].       (24)
```

The composition in (6)--(7) is its literal inverse: `F_z F_-z=Id` on the image, the square map reverses the positive radial root, and the component scales/coframe reverse L I. Hence

```
Psi(0,V(0))=(g0,phi0)(0)+X_raw.                          (25)
```

For the velocity, the time-dependent raw normalization has angular factors sech(t) and sech(t)^2; its explicit derivative at t=0 is zero for every raw state. The background is even in t, the radial root is time independent, and the flow has no additional explicit t dependence in the de Sitter frame. Differentiating the exact inverse identity at fixed z therefore gives

```
Psi_t(0,V(0))+Psi_V(0,V(0))V_t(0)
                  =partial_t(g0,phi0)(0)+Xdot_raw.       (26)
```

In the coordinates used in (7), `Psi_t(0,V)=0` for every state, but (26) keeps that explicit derivative visible rather than suppressing it by a convention. The elapsed initial pair in P1 is `(J_E V(0),-J_E V_t(0))`; (1) returns exactly (24), including its velocity sign.

The second physical time jet is supplied by the actual solution, not independently prescribed:

```
X_tt=Psi_tt+2Psi_tV V_t+Psi_VV(V_t,V_t)+Psi_V V_tt.        (27)
```

Every term in (27) is retained; in particular evenness of the frame's first derivative does not remove Psi_tt. The current `J2` and strong X14 regularity make (27) a strong H12 identity at zero. Sections 3--4 then give (4) there in H12 for exactly the physical pair (25)--(26). This is the required premise for deriving the initial normal derivative of the harmonic vector from the current canonical constraints.

## 8. Outer collar, angular symmetry and identical copies

P1 proves actual support control for the nonlinear solution: `U=0` on the fixed open outer collar `80ell/19<z<=6ell` for the whole elapsed slab. By (1), `F_z(0)=0` and (6)--(7), the physical reconstruction there is exactly

```
g=g0,             phi=phi0,             mathcal C=0.     (28)
```

The last equality uses equality on an open spacetime collar, hence equality of first jets. An isolated Dirichlet wall value alone would not establish it. All corresponding time derivatives vanish there as well.

The generator and all point maps are natural under the retained angular reflections, so Gamma invariance passes to the reconstructed tensors. No finite angular truncation or symmetry restriction on arbitrary differentiated tensor words is used. Apply the same reconstruction to each identical interior copy, with its own positive inward z and outward normal (18). Copy identification gives identical h and phi; each copy obeys the same outward laws (20)--(21). This is the symmetric branch, not independent nonsymmetric data.

For orientation bookkeeping alone, if one chooses a single transverse normal from the minus copy to the plus copy, it equals the minus outward normal and the negative of the plus outward normal. Thus the corresponding oriented values are `K_minus=sigma h/6`, `K_plus=-sigma h/6`, and the scalar normal values are `-sigma'/2,+sigma'/2`. With plus-minus jumps this gives

```
[K_ab]=-sigma h_ab/3,    [K]=-4sigma/3,
[K_ab-Kh_ab]=sigma h_ab, [n phi]=sigma'.                 (29)
```

This is consistent with `S_ab=-sigma h_ab` in the usual Israel sign convention. Equation (29) records orientations; it is not a new global distributional-shell existence theorem or an additional incoming boundary row.

## 9. Handoff and exact scope

The current conditional reduced solution therefore reconstructs, by the actual pinned point chart, to fields with (2)--(3), the literal reduced physical equations (4), and the complete primal shell laws (5). It has exactly the prescribed current completed physical Cauchy pair, retains the same-state past reflection, and is the exact background on the established outer collar. All these claims follow by explicit conjugacies and trace identities above, rather than importing the conclusion of the older generic gauge/evolution note.

The initial-data companion must supply `mathcal C=0` and its normal time derivative at the initial section for the current canonical data. The subsidiary companion must derive and uniquely solve the homogeneous system with (23) and its tangential Robin rows. Only after those joins may one replace (4) by the unreduced equation `mathscr E=0`. No evolved parameter differentiability, top-data stability, proper observer tube or empirical conclusion is supplied here. No new amplitude restriction or enlarged admitted interval has been introduced.
