# Canonical normalization clarification for the physical response handoff

15 October 2026 UTC. Separate audit clarification. This does not alter the
frozen accepted response manuscript, its proof notes or existing manifest.
The response manuscript reviewed here remains SHA-256
`a78f2239893fdd0e5ff653d4d5ce5cd9d3244c57c794ed2d3bb43bd4f5e1484b`.

The older action convention is

`S_s=integral [R/2-|dphi|^2/2-U] + GHY - common_shell_tension`.

Let `mu=sqrt(det gamma)`, let `kappa=+Lie_u gamma/2` on the spacelike
initial section, and write `v=u(phi)`. Its canonical densities are

`pi_s^{ij}=mu(kappa^{ij}-tr(kappa)gamma^{ij})/2`,
`varpi_s=mu v`.

The later gauge companion displays the uniformly doubled bulk convention,
with densities

`Pi^{ij}=mu(kappa^{ij}-tr(kappa)gamma^{ij})=2pi_s^{ij}`,
`Upsilon=2mu v=2varpi_s`.

This must be an overall action rescaling, including GHY and shell action;
doubling only the bulk while retaining a singly normalized shell action
would change the physical junction law and is not the model used here.
Uniformly multiplying the full action by two leaves every physical field
and shell equation unchanged.

In four spatial dimensions, put
`Q(pi)=pi^{ij}pi_{ij}-(gamma_{ij}pi^{ij})^2/3`. The half-action constraint
Hamiltonian density can be written

`H_s=2Q(pi_s)/mu + varpi_s^2/(2mu) - mu R_gamma/2`
`    + mu |Dphi|^2/2 + mu U`.

The doubled-action density is

`H_d=Q(Pi)/mu + Upsilon^2/(4mu) - mu R_gamma`
`    + mu |Dphi|^2 + 2mu U`.

Consequently, for `T2(gamma,pi_s,phi,varpi_s)=(gamma,2pi_s,phi,2varpi_s)`,

`H_d = 2 H_s composed with T2^{-1}`,
`M_d = 2 M_s composed with T2^{-1}`,

where `M_i=-2 D_j pi^j_i + varpi D_i phi` in each convention, using that
convention's consistently scaled momenta. Thus an exact canonical curve
`x_s(epsilon)` gives the exact doubled-variable curve
`x_d(epsilon)=T2 x_s(epsilon)`. This transports its tangent, all higher
parameter derivatives, support and reflection symmetry. The finite Sobolev
spaces and weighted inverse estimates are unchanged up to constant norm
equivalence. No new existence or invertibility hypothesis is required.

The geometric reconstruction explicitly agrees:

In the following formula `pi_s` without indices means its gamma-trace.

`kappa_ij = 2(pi_s,ij - pi_s gamma_ij/3)/mu`
`         = (Pi_ij - Pi_trace gamma_ij/3)/mu`,
`u(phi) = varpi_s/mu = Upsilon/(2mu)`.

The response theorem specifies its tangent through physical metric/scalar
master fields. Their ADM sectors must therefore be computed with one fixed
convention, or transported by `T2` as an entire curve before using the other
convention. Its physical first jet, chosen lapse/shift completion and shell
scalar response are then identical. The earlier gauge companion expressly
allows retaining the checkpoint's canonical-to-geometric convention
throughout and requires that canonical normalization and constraints change
together. That instruction is sufficient for this application.

Blindly inserting the numerical half-action momenta into the doubled-action
reconstruction or Hamiltonian would be incorrect. It would not represent
the same physical initial jet or exact constraint curve. No such insertion
is necessary or asserted in the response proof. The physical observable
uses the one common shell metric/scalar and receives no factor of two from
this action convention or from the two interior copies.

Verdict: consistent after the explicit constant identification above; no
blocker or new physical assumption. This is a clarification of how the
retained general initial-data theorem is instantiated from a physical first
jet, not a modification of the frozen mathematical claim.
