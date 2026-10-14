# Next target: complete one second-order state through shell contact

14 October 2026, America/Los_Angeles. Planning note outside the frozen tensor
packet. This specifies the remaining calculation; it asserts no solved
coupled or exact nonlinear contact theorem.

Choose the **copy-even aligned quadrupole**
\(A=B=\operatorname{diag}(1,-1,0,0)\), with
\(Y(n)=n_1^2-n_2^2\), the same smooth compact radial profile in displacement
and velocity, and the inherited fixed positive time scale. Put the profile
in a thin annulus just inside the actual shell. This choice passes the ten
initial consistency tests and has \(\zeta_1=0\) in the inherited longitudinal
first gauge. It is invariant under the coordinate reflections used by the
exact localized initial-data construction. Both copies have identical first
fields; choose identical second fields and initial corrections.

This is a tractable contact target because first-order master packets and
their coupled linear junctions are already known, the quadratic angular
products occupy a finite harmonic space, and the tensor coefficient is now
solved. It retains a nontrivial scalar boundary source. Do not replace it
by the copy-odd graph case or a different background during the calculation.

## Known inputs that can be reused

1. The actual registered Einstein–scalar background, its fixed potential and
   tension, and its controlled contact domain. A compact packet can be chosen
   with \(t_c=\log(b/R_{\max})\) and
   \(t_c<t_*<\operatorname{arsinh}(\tan\tau_*)\), where the inherited
   first-order solution is known on its conformal slab \(|\tau|\le\tau_*\).
2. Exact nonlinear initial data and sufficiently short **precontact**
   evolution for this reflection class, plus a localized solution of the
   second-order initial constraints. Tensor metric/momentum initial data
   are homogeneous constraint directions and can be set to zero separately.
3. The complete second-order *material* matching and Israel/scalar junction
   algebra. On identical copies, there is still one shell tension; dividing
   the two-copy equations gives
   \[
   K_{2,ab}-\frac{\sigma_0}{6}h_{2,ab}
      -\frac{\sigma'_0}{6}\psi_2h_{0,ab}
     =\frac{\sigma'_0}{3}\psi_1h_{1,ab}
      +\frac{\sigma''_0}{6}\psi_1^2h_{0,ab},
   \qquad
   p_2+\frac{\sigma''_0}{2}\psi_2=-2\psi_1^2.
   \]
   Here \(h_1=2\Psi_bh_0\), \(\psi_1=\chi_b\),
   \(\sigma''_0=4\varphi_0\), \(\sigma'''_0=4\). The metric source is
   angular trace; the scalar source cannot be dropped. More precisely,
   with \(\chi_b=v_b(t)Y\),
   \(Y^2=F_4+\tfrac12F_2+\tfrac16\), so its \(l=4\) coefficient is
   \(-2v_b(t)^2F_4\), nonzero wherever \(v_b\ne0\). Nonvanishing for a
   particular packet or at the causal support edge remains to be checked.
4. The complete longitudinal-gauge quadratic TT bulk source
   \(C=16u^2+8v^2/3\) for \(\Psi=u(y,t)Y\), \(\chi=v(y,t)Y\), and the
   forced tensor solution with
   \(\mathcal L_TH=2C/(\rho^2\cosh^2t)\), \(H_y|_b=0\).
   Its auxiliary energy, source normalization, corner compatibility and
   controlled finite-time domain are supplied by the accepted tensor note.
5. The exact residual identities
   \(D_bI^b{}_a=\sum_iE_{n_i a}+P D_a\varphi\) and the normal-force
   identity. These check boundary consistency. They do not yet provide
   a hyperbolic boundary theorem for the remaining fields.

## Minimum subsystem still required

**The complete forced angular scalar equations.** Compute the quadratic
Einstein residual components \(E_{yy},E_{yt},E_{tt},E_{yi},E_{ti}\), the
angular trace and scalar-derived trace-free part of \(E_{ij}\), and the
scalar-field equation. Keep the most general second-order scalar-harmonic
metric fields on the two-dimensional \((y,t)\) base, including lapse,
shift, angular trace and angular Hessian amplitudes, plus the scalar-field
coefficient. Use the linear constraints and a stated gauge to identify a
closed forced hyperbolic subsystem. The first-order scalar-master
reconstruction is a sufficient linear ansatz; an arbitrary quadratic source
is not automatically in its image. In particular, do not assume a single
forced copy of that master supplies every scalar-derived gravitational
polarization.

For this single \(Y\), scalar products such as \(Y^2\) and \(|DY|^2\)
have degrees \(l=0,2,4\). Retain the monopole and treat its constraint/gauge
degeneracies directly, without division by a harmonic eigenvalue that
vanishes. The solved \(l=2\) angular TT coefficient does not settle these
scalar-derived sectors.

**The vector sector, or a proved zero-vector reduction.** In the selected
first gauge, quadratic mixed angular one-forms are expected to be
proportional to \(Y D_iY=\tfrac12D_i(Y^2)\), suggesting zero coexact-vector
forcing. This is a calculation to verify for the complete bulk equations,
boundary conversion and initial constraint correction. Reflection symmetry
alone does not prove it. If all such projections and the chosen initial
vector data vanish, the identically zero coexact sector is a valid selected
solution. Otherwise its forced equations and shell conditions must also
be solved.

**Scalar boundary geometry and gauge conversion.** Keep the second normal
embedding \(\zeta_2\), normal lapse/shift, induced scalar trace and the
pullback/evaluation terms until a complete boundary gauge is fixed. Unlike
the angular TT projection, a scalar Hessian of \(\zeta_2\) survives here.
Convert the displayed material equations to that gauge. In a bulk wave
gauge, transform the complete quadratic bulk source together with the
material boundary source; the longitudinal source is not gauge-universal.
Determine the boundary displacement/normalization rather than inserting
\(\zeta_2=0\) without checking all scalar Israel components.

**Gauge and constraint propagation at the shell.** If using a fixed-shell
generalized harmonic reduction, there are fifteen metric wave unknowns and
one scalar wave unknown before angular decomposition. The eleven physical
Israel/scalar laws must be completed by a compatible boundary coordinate
and constraint prescription. That arithmetic is only a starting count:
compute the actual principal rank and incoming characteristic subspace.
Do not append five harmonic-vector conditions blindly to an already fixed
set of sixteen boundary rows. Derive the second-order subsidiary gauge
equation from the differentiated Bianchi/scalar identities and prove that
the selected boundary prescription supplies zero incoming gauge/constraint
data when the initial constraints vanish. The four tangential Codazzi
relations and normal constraint must be retained in this check.

**A boundary estimate and construction for the resulting subsystem.**
Show that its physical boundary flux and its incoming gauge/constraint flux
are controlled on the same finite contact slab. A dynamic scalar boundary
law needs its actual boundary mass, potential and corner hierarchy. A
negative bounded potential can often be handled by a bulk shift and trace
estimate; a wrong principal boundary sign or an uncontrolled incoming mode
cannot. The tensor estimate and a vacuum timelike Einstein-boundary theorem
do not provide this missing coupled estimate.

## The next falsifiable calculation

Start with the **even \(l=4\) scalar sector** generated by \(Y^2\). Extract
its complete bulk residual source and its scalar junction source,
then substitute a general scalar-harmonic second-order metric/scalar ansatz.
Eliminate only justified constraints and gauge redundancies. Produce the
closed evolution equations and the converted shell boundary principal
matrix, and freeze that matrix at a regular shell point for its characteristic
rank and energy-flux calculation.

The calculation succeeds if it identifies every independent scalar
polarization needed by this source, supplies exactly the required incoming
conditions, and closes the constraint/gauge flux with the correct scalar
boundary source and all corner identities. A nonzero residual outside a
proposed one-master ansatz, a boundary rank defect, or an incoming
constraint mode is a concrete failure of that proposed reduction. Such a
failure requires revising the ansatz or boundary prescription; it is not
evidence against the already proved tensor result and not a cosmological
observation. After this sector closes, repeat for the \(l=2\) and exceptional
\(l=0\) scalar sectors and verify the candidate zero coexact-vector reduction.

Only after these coefficients jointly satisfy every bulk equation,
configuration match, Israel/scalar junction, initial constraint and boundary
constraint/gauge condition may the result be called a **complete
second-order contact solution**. Even that would remain a perturbative
coefficient construction; exact nonlinear continuation through contact
would still require a separate nonlinear existence argument.
