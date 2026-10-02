# Independent static bridge and source-specification cross-review

2026-10-02 UTC. Reviewed staged
`junction-next/ACTION_AND_JUNCTION_SPECIFICATION.md` and
`static-bridge/STATIC_DESITTER_SOURCE_BRIDGE.md` against the actual inherited
action and source code. All calculations in this review are algebra or matrix
reanalysis; no new background BVP, sensitivity integration, quantum sum, state
expectation, mode evolution or physical source has been executed.

The displayed action/junction specification needs no equation or sign
correction. The static bridge's expansion, Euclidean scaling formulas and
archived matrix forcing also agree with an independent derivation. The tiny
archived scalar-current radius response is not an accurately resolved physical
susceptibility. The final staged static note now states this limitation and
requires prospective derivative checks rather than treating it as a result.

## Static expansion checked independently

At phi_b=1+eta the original registered d=0 tension and bulk potential obey

    sigma = 2/3+delta(1+c)+delta c eta+2eta^2+O(eta^3),
    sigma_phi = delta c+4eta+2eta^2,
    U = -2/27+(14/9)eta^2+O(eta^3).

Let source-weighted shell values be s=kappa5^2 rho and t=kappa5^2 J, and
define d=delta c, nu=d+t, tau=delta(1+c)+s. In this paragraph rho,J are
frozen supplied shell values, and delta,s,t have a common formal small order.
The normal points toward the shell, so the scalar junction is
w=-(sigma_phi+t)/2. The regular linear cone mode has asymptotic normal
logarithmic derivative 14/9. Matching it gives

    eta = -9nu/64+O(epsilon^2).

The independently checked Gegenbauer polynomial C_14^(2)(cosh(y/9)) solves
the linear radial equation. Its large-radius limit supports this slope; the
finite-curvature response requires the actual radial solution.

Before substituting eta, the second-order static constraint reduces to

    H^2 = tau/27+tau^2/36+d eta/27-nu eta/6-nu^2/48+O(epsilon^3).

The eta^2 coefficient cancels between the metric, scalar and potential terms.
An unknown second-order coefficient in eta also cancels at this order. Thus

    H^2 = tau/27+tau^2/36-d nu/192+nu^2/384+O(epsilon^3).

This equals tau/27+tau^2/36+(t^2-d^2)/384. It is a response expression
for supplied shell values, not a general prescription for differentiating a
quantum action. Generic source variations with phi and H must be included
according to their assigned perturbative order.

For the special action Gamma=-int sqrt(-h)[V0+J0(phi-1)], rho=V0+J0 eta and
J=J0. Writing s0=kappa5^2 V0 and t=kappa5^2 J0 adds t eta/27 to the
constraint. The result becomes

    H^2 = tau0/27+tau0^2/36-nu^2/384+O(epsilon^3),
    tau0=delta(1+c)+s0.

This potential-only result must not be generalized to curvature or general
state-dependent sources. For Gamma_local=int sqrt(-h)[-V+FR], at constant
phi and exact de Sitter,

    rho=V-6FH^2,  J=V_phi-12F_phi H^2,
    (partial rho/partial phi)_H=V_phi-6F_phi H^2.

J is not the fixed-H derivative of rho when F_phi is nonzero. All three
expressions and the distinction were independently verified.

## Euclidean-state route checked, not executed

For S4 radius 1/H and x>0, lambda_l=H^2 l(l+3)+x and
d_l=(l+1)(l+2)(2l+3)/6. The Euclidean-state construction has a strictly
positive mass gap. The minimally coupled x=0 zero mode prevents using these
formulas as the same ordinary gapped vacuum prescription. The prospective
branch also requires H^2>0.

The formal renormalized Gamma_E=Vol W(H,x), with
Vol=8pi^2/(3H^4), gives

    Q=2 W_x,
    rho=W-H W_H/4,
    rho_x=Q/2-H Q_H/8.

Holding x and the reference/renormalization data fixed during metric variation
is essential. Direct insertion of the Euclidean image
W_local=V-12FH^2-144alphaH^4-24gammaH^4 produces the Lorentzian local
rho and Q above. Constant alpha R^2 and Euler terms produce no static local
source on exact de Sitter; that statement does not remove the full anomaly.

The inherited matched trace expression specializes correctly to

    -4rho=-xQ+A_r,
    A_r=[(x-r)^2/2-2(x-r)H^2+(29/15)H^4]/(16pi^2).

It is conditional on establishing the same renormalization prescription for
the prospective Euclidean calculation. It is not a evaluated quantum source
or a valid identity for an arbitrary finite harmonic cutoff. The staged note
properly requires an independent metric variation, matching, ordered limits,
resolved tails and source derivatives. Euclidean/BD state selection differs
from the smooth control's exactly static incoming plane-wave preparation.

## Archived Jacobian provenance and limits

The archived producer is
`research/HDBLAST_CHECKPOINT_20260927/stability_gauge_invariant/spectrum.py`,
lines 100–107. It uses centered finite differences of step 1e-6 in
(log10(abs(eta_h)), y_b), holding the sign of eta_h fixed. The residuals in
`gi_core.py:88` are unscaled

    E1=q-sigma/6,  E2=w+sigma_phi/2.

The separate root polish divides residuals by max(delta,.001); that scaling
does not appear in this archived Jacobian. Adding sources gives
E1-s/6 and E2+t/2, hence the forcing matrix diag(1/6,-1/2) is correct.
Inversion of the saved matrix is valid matrix algebra, not an independent
reproduction of its derivatives.

The continuum endpoint column derives from q'=-w^2/4-U/6-q^2 and
w'=U_phi-4qw, with vacuum junctions imposed:

    (E1)_y=-H^2,
    (E2)_y=(phi''/phi'+sigma_phiphi/2)w=Bw.

At delta=.001 the archived endpoint column differs from those exact formulas
by approximately (7.37e-12,-4.84e-13). The matrix condition number is about
16.94 in these specific coordinates; it is not a physical response condition
number in phi_b and H^2.

The first-column metric derivative is -6.938893903907228e-12. With h=1e-6
this represents a residual difference of approximately -1.388e-17, at rounding
scale for the terms being differenced. The resulting tiny dy_b/dt value,
about -8.5e-5, depends entirely on that entry. A moderate overall condition
number does not establish the accuracy of an individual near-zero component.
It must remain unresolved until independently refined sensitivities are
available. Endpoint output derivatives phi_ell and (H^2)_ell are also needed
to turn coordinate changes into physical branch changes; the archive does
not certify them.

The static note correctly avoids a local existence certificate, a quantum
stability claim, a new static solution, or a physical first-order response
with resolved error bars. Its proposed loop expansion is a useful route once
the source and sensitivity inputs are prospectively fixed.

## Action/junction specification cross-review

Checked the moved local metric terms and the retained -kappa5^2 F_phi R scalar
term; proper-time/conformal-time/bulk-chart mode transformation and Wronskian;
physical versus kappa5^2-weighted code source normalization; and signed local
bulk, Codazzi and Weyl relations. All displayed signs and factors agree with
the independent action derivation. Existing B1 helper limitations are narrowly
scoped. The difference between the original delta=.001,d=0 branch and B1's
tuned delta=.1,d=-3.106933495673783 branch remains explicit. The source
specification also correctly requires a larger constrained system or justified
order reduction rather than a direct higher-derivative adapter insertion.

## Evidence

`independent_static_bridge_checks.py` derives and verifies the scalar cone
equation, expansions with an unknown second-order scalar response, local
Euclidean variations, general source integrability, anomaly, endpoint column
and archived inversion. It passes 16 checks and rejects two deliberate
mutations under both ordinary Python and `python -O`. Detailed evidence and
the matrix caveat are recorded in `INDEPENDENT_STATIC_BRIDGE_CHECKS.json`.
This review does not rerun the historical radial solver or spectrum.
