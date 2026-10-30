# Composition of the natural-field implementation joins

30 October 2026. The exact coefficient calculations here bind the current
first reconstruction and source maps to analytic identities already used in
the retained construction. They do not replace its PDE existence and
regularity hypotheses. All statements concern the physical collar with
nonzero rho and p, using outward proper y and inward scaled r with
`partial_y=-(ell*rho)^(-1) partial_r`.

## First natural bulk equations

Let `vartheta=cosh(t)^(-3/2)`, `H=tanh(t)`, and let Y be the original
degree-two harmonic. Put

```
P=alpha*vartheta*W*Y,  Psi=P/rho^2,
h_yy=-4Psi, h_y,mu=0, h_mu,nu=2rho^2 Psi gamma_mu,nu,
chi=-3 P_y/(rho^2 p),
Q=2(V'/p)^2-10a(V'/p)+(35/4)a^2-V''-p^2/6,
Rmaster=W_xx-W_rr+ell^2[rho^2 Q-(35/4)H^2]W.
```

For arbitrary smooth longitudinal fields the independently reconstructed
natural residuals factor through

```
M[P]=P_yy-(2a+2p_y/p)P_y+rho^-2(Box_gamma+4)P-(2/3)p^2 P,
M[P]=-alpha*vartheta*Rmaster*Y/(ell^2 rho^2),
E1_yy=2M/rho^2, E1_y,mu=0, E1_mu,nu=-gamma_mu,nu M,
S1=-3(D_y+2a)M/(rho^2p).
```

Thus the actual master equation on a neighbourhood supplies E1=0 and S1=0
on that neighbourhood, with the needed derivatives of E1. Setting Rmaster
to zero only at a point does not set its radial derivative to zero, and
does not suffice for the scalar equation. The scalar variation includes
the nonharmonic contraction `-DC[h].dphi0`; the longitudinal representative
cannot be put into a harmonic reduced operator as if DC[h] vanished.

The actual alpha coefficient here is `(rho/rho_b)^(3/2)*p/p_b`; it is
distinct from the Sturm–Liouville weight called alpha in the earlier
ambient proof. The separated calculation retains the radial +5/2 and
temporal -5/2 terms until they cancel. Ordinary physical derivatives act
on rho as well as on W. No interval endpoints or unknown Taylor remainders
are differentiated as though they were exact functions.

The previous exact first-wave audit binds the implemented harmonic field
to the plus-sign Lie pullback of this longitudinal pair. Natural covariance
on the exact background transfers the natural first equations to that
pair. This step retains the standard natural covariance identity. The
gauge-wave PDE supplies harmonicity, not the natural equation by itself.

## Physical first shell conditions

At the shell define the actual dynamic trace residual

```
Dtrace=W_xx-(ell/mu)(W_r-ell*gamma*W)
                +ell^2[5/2-(35/4)H^2]W.
```

The longitudinal metric junction vanishes by the same mixed
reconstruction constraint. The outward one-copy physical scalar junction
amplitude obeys

```
P1=3alpha*vartheta*(Rmaster-Dtrace)/(ell^2 rho_b^4 p_b).
```

These are identities of residuals before the equations are imposed. They
do not construct a new solution of the trace equation. The shell-fixing
gauge has f=0 and b_r=d_r=0 there, with tangential values b,d left free.
The physical shell residuals transform naturally, with the unit-normal
and lapse terms retained. A change of tangential labels is not a change
of the physical boundary law. Shell-only conditions are differentiated
only tangentially, not normally.

## Quadratic natural bulk source identity

For `E=Ric-dphi*dphi-(2/3)Vg`, `S=Box(phi)-V'` and
`B_g(T)=div_g T-(1/2)d tr_g T`, the background identity is
`B_g(E)=-S dphi`. Its affine second variation on an exact background is

```
B0(QE)+p Qphi_general e^y+2 deltaB_q(E1)+2 S1 dchi=0,
deltaB_q(T)_b=-q^{ac}nabla_a T_cb+(1/2)q^{ac}nabla_b T_ac-C1^d T_db.
```

Within the registered seven-profile quadrupole ansatz, the new map check
retains arbitrary coherent metric jets through order three and scalar jets
through order two. It independently constructs E1
and S1, includes their residual corrections, and checks all eight modal
covector rows with complete angular reconstruction. The background
Hamiltonian relation and de Sitter quotient are explicit.

The natural homogeneous specialization follows when E1 and its first
prolongation vanish and S1=0. The preceding local first-field factorization
provides the reason those equalities apply to the retained smooth actual
master, conditional on that actual solution. Harmonicity on a neighbourhood
additionally permits Qphi_general to specialize to Qphi. Harmonicity alone
does not remove the other first residuals.

The tensor input is natural QE, whose angular J/H coefficients are already
in the registered trace-free basis. It is not QEred or metric forcing and
must not be trace-reversed twice. The two transverse-traceless source
channels lie in B0's kernel; these eight rows cannot validate those
components. Their separate full-angular reconstruction checks remain
relevant.

## Remaining scope

The exact source/proof joins above do not change any numerical candidate,
source map, physical law or old error certificate. They give no improved
absolute or pointwise error bound. Actual background/master/gauge solution
existence, smoothness, corner compatibility and the subsidiary uniqueness
arguments remain the identified analytic inputs. The generic Bianchi
coefficient check does not supply a new numerical third-jet adapter for
the actual first fields or a numerical enclosure of all its radial
derivatives.

The next distinct implementation target is the reduced-source identity
with full-metric-lowered QC_actual, its covector wave operator and the
general PC1 terms retained before harmonic specialization. It must keep
the same first field and all needed prolongations; equality of QC values
at a point is insufficient inside a wave operator. This is separate from
the retained analytic subsidiary proof, which already exists.

No exact finite-amplitude nonlinear family, effective nonlinear remainder,
observational agreement or higher-dimensional Big Bang origin follows.
Internal AI reviews are implementation reviews, not external specialist
assessment. Frozen manuscripts and private editorial history are outside
this capsule.
