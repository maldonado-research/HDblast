# Next target: an angular-uniform estimate for nonlinear contact

14 October 2026. Planning note. The complete selected second Taylor
coefficient is constructed in the new checkpoint. Exact nonlinear
continuation through shell contact remains unproved.

The immediate target is a boundary estimate for the actual coupled
Einstein–scalar system whose constants remain controlled as the angular
cutoff grows. The finite-mode characteristic proof supplies the selected
second coefficient in a fixed finite angular space; it does not supply this
uniform estimate. Exact nonlinear
products generate additional angular modes, so a fixed truncation cannot
be treated as an exact solution of the full equations.

Keep the registered background, canonical scalar potential, field-dependent
tension and identical-copy shell. Reuse the exact reference harmonic
boundary completion: four mixed metric coordinate conditions, ten Israel
conditions, the scalar Robin law and one normal harmonic condition. The
tangential subsidiary conditions must still follow from Codazzi and the
scalar law. These tangential conditions are not additional primal rows.

The scalar RW boundary calculation gives a useful concrete starting point:
here \(N\) is the scalar-harmonic amplitude of \(h_{2,yy}\), and \(w\)
is the scalar-field amplitude. This \(N\) is distinct from the normal lapse
used in the exact wave-chart boundary formulas.

\[
 \begin{pmatrix}N_y\\w_y\end{pmatrix}
 +\begin{pmatrix}2a_b&-4p_b/3\\-p_b/2&\sigma''_0/2\end{pmatrix}
 \begin{pmatrix}N\\w\end{pmatrix}
 =\begin{pmatrix}F_N\\F_w\end{pmatrix}.
\]

The exact algebraic weight \(\operatorname{diag}(3,8)\) makes this boundary
matrix symmetric. That fact alone does not establish that \(N,w\) are two
closed hyperbolic physical channels or provide a bulk energy. Derive their
full principal operator and constraint reconstruction before using the
weight in an estimate. Treat the exceptional monopole separately, without
division by a zero harmonic eigenvalue.

The next calculation should:

1. Derive the full coupled principal bulk and boundary symbols in the same
   gauge, including the scalar, the metric trace/lapse block and all gauge
   fields. Track angular-frequency dependence and the constraint subspace.
2. Examine the glancing regime with the correct energy weights. A vanishing
   unweighted Neumann symbol by itself proves neither ill-posedness nor an
   estimate. Look for uncontrolled incoming constraint modes, derivative
   loss or growth that increases without bound with angular frequency.
3. Construct a justified positive control norm and source estimate, with
   constants independent of the cutoff. Include boundary terms only when
   their flux identity or a trace theorem supports them. The nonsymmetric
   full Robin matrix cannot be assigned the scalar energy cancellation
   without proof.
4. Verify that inhomogeneous boundary data, corners and the subsidiary
   system obey the same regularity estimates. The actual Taylor sources
   must be used; changing first gauge changes both bulk and boundary sources.
5. Only after a suitable linear estimate closes, prove a quasilinear local
   construction with compatible initial data, Lorentzian signature, a fixed
   timelike shell and constraint propagation. Specify the Sobolev or other
   regularity class and any derivative loss explicitly.

A successful uniform estimate would justify attempting exact nonlinear
contact existence. An uncontrolled incoming constraint mode or genuine
uncontrolled high-frequency growth would falsify that proposed formulation
and require a revised gauge or boundary reduction. Neither outcome changes
the established fixed-finite-mode second-coefficient theorem or supplies a
cosmological observation.

Fournodavlos–Smulevici's vacuum umbilic method and the coupled scalar/master
literature are useful comparisons, with their original dimensional and
background restrictions retained. No new quantum, string or matter sector
is needed to test this mathematical target. Such additions would require
their own physical motivation and new field and shell equations.

External mathematical review should examine the frozen coefficient proof
and its source formulas before broader novelty or submission claims. This
planning note asserts no new uniform or nonlinear theorem and does not
request or create a Zenodo version or APS submission.
