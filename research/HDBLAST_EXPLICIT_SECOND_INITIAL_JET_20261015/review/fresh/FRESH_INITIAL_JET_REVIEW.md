# Fresh independent review: explicit compact second initial constraint jet

15 October 2026. **Accepted within the stated conditional initial-constraint
scope. No mathematical defect found.** The reviewed final source pins are in
`FRESH_INITIAL_JET_ACCEPTANCE.json`. The accepted derivation is SHA-256
`2cca9c9b87dfa5c2c21712ae81ea0c91d01b939b902e82d54dd62a6eae6bf0ed`,
the producer is `954eeb402c3f06c72951bdfb6051d3f9e346b43ae2487324d67370021262425e`,
and the capsule README is
`c61051888de1b73f9958802f4ec41575347ebadec42cf5062e1961b2b739a06f`.

The registered physical constraints are correct for
`K=+Lie_n(gamma)/2`, `Pi=n phi`, and the retained Einstein–scalar normalization.
At the background `K=Pi=0`, the linear Hamiltonian and momentum corrections
decouple. The construction explicitly specifies all four physical Cauchy
second fields; the second scalar and scalar normal velocity are allowed to
vary. No canonical momentum-density interpretation is used.

I independently rederived the full Hamiltonian source from the radial
foliation geometry and the scalar/K/Pi expansions. With
`j=a+eps*(u_y+2au)Y+6eps^2*a*u^2*Y^2`, the radial curvature Hessian is
`[-24u_y^2-288auu_y-24uu_yy-(192a_y+384a^2)u^2]Y^2`.
The angular conformal-curvature contribution is
`rho^-2*(-208u^2*Y^2+12u^2*G_Y)`; the angular lapse term contributes
`rho^-2*(-128u^2*Y^2+24u^2*G_Y)`. Together they give the stated `-336,+36`
coefficients. The scalar gradient gives
`-(32p^2u^2+16puv_y+2v_y^2)Y^2-2v^2*G_Y/rho^2`, the potential gives
`-2U''v^2Y^2`, and the K/Pi terms give
`-(12u_t^2+2v_t^2)Y^2/rho^2`. Every row of (6) is retained.

For the momentum source, independently raising the affine physical covariant
K gives `P^y_y=-3eps*sY+6eps^2*usY^2` and
`P^i_j=12eps^2*usY^2 delta^i_j`, where `s=u_t/rho` and
`P=K-tr(K)Id`. Its full divergence, including the perturbed connection, yields
the stated radial `G` and angular `F*d(Y^2)`. In particular the
`-6u_yu_t/rho` radial term and the angular Hessian factor are correct.
`K2=(F/3)Y^2 gamma0`, `Pi2=(G-F_y)Y^2/p` cancels all four momentum rows.

The angular identities are exact: `Delta N=4-8N` and
`Delta(Y^2)=8N-24Y^2` give degrees 0, 2, 4 with eigenvalues 0, 8, 24.
The three Hamiltonian sources project with exactly the displayed coefficients.
For `delta gamma_ij=f W gamma_ij`,
`delta R=-3f_yy-12af_y+2(lambda-3)f/rho^2`.
The compensator `zeta=-3f_y/(2p)` cancels both radial derivative terms using
the **spacetime background** scalar identity `p_y+4ap=U'`.
The denominators `-6,10,42` and source sign are correct; the monopole is
included. The exceptional dipole eigenvalue 3 does not occur in this source.
The actual retained box proves `-7/50000 <= p <= -3/25000 < 0`, so all
divisions on the packet collar are justified.

The physical-Cauchy-affine and spacetime-affine paths differ exactly as stated.
The latter induces `K2=-2uY K1`, `Pi2=-2uY Pi1`, whose linear momentum
contribution is `+6u_yu_tY^2/rho` by the first momentum identity.
The coordinate conversions follow from `dot gamma=2NK`, `dot phi=NPi`:
`dot gamma2=2rho(K2+2uY K1)` and
`dot phi2=rho(Pi2+2uY Pi1)`. Choosing `N2=0` gives
`g_tt,2=-2rho^2u^2Y^2`. There is no missing angular conversion term.

I separately verified the actual packet-maximum obstruction with the exact
authenticated coefficient endpoints, deriving rather than assuming the needed
`u_y,u_yy,v,v_t` caps. It proves `2.3e9<|v_y|<3e9` and
`Q_H<-10579999999995000000<-1e19`. Thus the obstruction is about the actual
specified packet, rather than a generic rational jet. The vanishing of K0 and
Pi0 explains why lapse-induced second velocity changes cannot remove this
Hamiltonian obstruction.

All profiles use finite derivatives and smooth coefficients on the collar.
The bump is flat at both annulus endpoints; the profiles extend by zero
smoothly with the same closed support `[ell,2ell]`. They vanish on open initial
shell and vertex collars, hence all their physical corner jets vanish there.
The postulated local compatible gauge completion remains a separate step.

The fresh standard-library acceptance checker proves 18 exact all-point
Laurent-polynomial identities supplementing these analytic derivations, then
executes the source checker in normal and optimized Python with explicit
`--background`. Each run passes all 48 finite coordinate/sector checks and
rejects all seven specified errors for their intended reasons without output
reports. The fresh outputs match the retained source receipts byte for byte.
The fresh acceptance checker itself also gives byte-identical reports under
normal and optimized Python. All source hashes are rechecked after replay;
scratch and selected output files live outside the source capsule.

The README accurately states the scientific scope. This is an explicit compact
second differentiated initial constraint jet under the retained actual-branch
assumptions. It gives no exact nonlinear initial family, second-field norm
enclosure, completed gauge bridge, coupled second evolution, U2 amplitude,
invariant force or nonlinear remainder. It does not alter the earlier frozen
checkpoint.

For a fresh capsule layout with the producer at its root:

```sh
python -I -S -B review/fresh/check_initial_jet_acceptance.py --root . --output /tmp/fresh-initial-acceptance.json
```

For the local preassembly source tree, also pass `--readme` pointing to the
separately staged packaging README. The checker writes only `--output` and
uses the system temporary directory for all subprocess scratch.
