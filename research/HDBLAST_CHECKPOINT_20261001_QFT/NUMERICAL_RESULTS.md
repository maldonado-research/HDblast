# Executed finite-band numerical results

Registration commit: 99dd818755c493127ae4b4e0cad7bbf21fd141de.
First executed source: c71054d10c4ecd6fd993251946a0e9d02119903b.
[Successful first run 36938696104](https://github.com/maldonado-research/HDblast-archive/actions/runs/36938696104).

## Setup and controls

The original source-free Y=0 A1 histories have 305 records each and later terminate non-finite. The selected endpoint s=6.9 is earlier than their last finite near-shell constraint records. Constraint checks are intermittent. The background comparison is spatial-resolution sensitivity of these archived data, not a new continuum certificate.

The finer Hermite crossing is s=1.4861216315, a_star=4.0424956912, q=128.31912836 and k_scale=a_star sqrt(q)=45.79259617. The fixed primary band is kappa=[0.05,3], with 32 Gauss-Legendre nodes; one case uses 64.

Both exact state families use the same operator. The in-state has the specified WKB-amplitude derivative at s=0. The out-reference is an exact solution evolved backward from physical-Hamiltonian data at s=6.9. It is not a reference reset at every intermediate time. Its Bogoliubov coefficients should be constant along the history.

All ten cases completed with no errors and passed their coded controls. The maximum paired path residual across them is 1.7309965e-7, also below the registered 1e-5 criterion even though the pass flags directly gate endpoints. The maximum raw path residual is 1.1923858e-9.

For the primary case, occupation agreement after mapping the identical state is 7.98e-9 in canonical proper time X and 1.40e-9 in conformal time u, across the four specified modes. The conformal null control has occupation <=1.05e-16 and relative mode error 1.77e-8. Exact-pair alpha/beta drift is below 8.4e-11. Endpoint relative energy equals the Hamiltonian particle-energy integral to 9.0e-11 relative.

Work integrals are auxiliary quadratures evolved alongside the modes with the same DOP853 controller. Energy is independently reconstructed from mode amplitudes/velocities. This is stronger than an algebraically substituted differential identity; it is not a separate post-processing quadrature. Independent RK4 source and comparisons are provided separately.

## Endpoint table

Number is in H0³; energy, pressure and current are in H0⁴. Delta quantities compare the specified exact states. Particle energy equals endpoint Delta rho for the Hamiltonian out basis, within numerical error.

| Case | Number density | Particle energy | Delta current | Delta pressure |
|---|---:|---:|---:|---:|
| primary | 0.76474290 | 41.116871 | 84.612524 | -4.3843727 |
| tight | 0.76474290 | 41.116871 | 84.612524 | -4.3843727 |
| quadrature64 | 0.76468790 | 41.113889 | 84.648351 | -4.4073141 |
| band2p5 | 0.76466353 | 41.112532 | 84.716561 | -4.4459162 |
| band3p5 | 0.76458513 | 41.108434 | 84.638665 | -4.4059507 |
| coarse | 0.76472138 | 41.115707 | 84.324170 | -4.2307664 |
| cubic | 0.76472720 | 41.116022 | 84.602885 | -4.3799234 |
| hamiltonian_in | 0.75751543 | 40.727366 | 78.626975 | -1.4555028 |
| start0p2 | 0.76468463 | 41.113724 | 84.626320 | -4.3945401 |
| out_end5 | 1.1406513 | 61.382453 | 118.19201 | -2.3465554 |

Each variation changes a specific numerical or physical assumption. In particular, a new initial state, band or out slice changes the state comparison.

The tighter integrator changes the energy by 4.3e-12 relative. Doubling quadrature changes it by −0.00725%; the alternate spline by −0.00207%; the coarser background by −0.00283%. Pressure/current are more phase-sensitive. For example, the coarser background changes the current by −0.341%.

The Hamiltonian initial preparation changes the count by −0.945% and current by −7.074%. Moving the WKB-amplitude initial slice to s=0.2 changes the count by −0.00762%, but still defines a different state.

Comparing out time 5 with 6.9 requires dilution to be accounted for. The physical number density differs by about +49%, whereas the comoving number differs by **−0.1510%**. Do not describe the density change as a 49% instability in particle production.

Band-edge variations change the state family and also move the 32 quadrature nodes. Their small nonmonotonic count differences do not establish a continuum cutoff limit. The spectra contain oscillatory interference, so retain the 64-node comparison and avoid claiming all displayed digits are physical precision.

## The older estimate and the actual force

The matched-band Gaussian number is 0.7799871968; the actual count is 0.7647428950, a ratio 0.9804557026. The earlier local correction factor is 0.9809491460, giving actual/corrected estimate=0.9994969735. This ~0.05% agreement concerns the integrated number for this declared case, not pointwise spectral agreement or the full current.

The particle-only current is 76.25370585 and its coherent addition is 8.35881776, giving 84.61252360. Coherence supplies 10.96% of the gas estimate, or 9.879% of the total current. The pressure correction is −4.46340337 and reverses the +0.07903063 gas pressure. Occupation alone therefore does not determine the force and pressure required by the shell equations.

The relative endpoint equation of state is −0.106632; the particle-gas value is 0.0019221. Mean endpoint energy per particle is about 53.7656 H0. This test contains no decay or thermalization calculation and does not produce a radiation-fluid source.

## Independent numerical implementation and endpoint-basis test

A separately written JavaScript RK4 implementation uses the same archived Hermite data, with 100,000 and 200,000 uniform steps for four registered momenta. Maximum occupation change is 2.1e-11. The independent 32-node calculation at 200,000 steps agrees with DOP853 occupations within 1.9e-11. Direct endpoint stress reconstruction agrees within 3.5e-8 for Delta rho, 9.4e-10 for Delta p, and 6.1e-8 for Delta J. These comparisons also test the squeezing phase, which occupation agreement alone cannot check.

After inspecting the first result, an explicitly exploratory projection onto endpoint WKB-amplitude data was added. At kappa=2.995964, occupation changes from 5.93812e-7 in the Hamiltonian basis to 3.77826e-10, a factor about 1,572. The whole-band comoving number changes by −0.06354%.

A slowly varying WKB-amplitude mode has real derivative offset gamma=(3H+omegadot/omega)/2, giving Hamiltonian occupation gamma²/(4omega²) even in that local adiabatic description. The observed tail is therefore strongly sensitive to the finite-time basis. The alternative tail is not declared the unique physical spectrum or a UV-renormalized result.

## Scale screen

The represented primary band reaches physical momentum 137.3778 H0 at its upper edge, with peak mass about 53.6102 H0. The sampled peak total frequency at the band edge is 146.0922 H0.

If every represented mode energy is required to remain below M5, this supplies the necessary screen b=(H0/M5)³ <=3.2072e-7. This is not a UV certificate, and it is not an affordability result for an absolute source. In particular, the earlier characteristic mass-only screen cannot automatically be reused for this full band.

## Limits that remain

No absolute vacuum stress/current or covariant finite matching was computed. No matter source was inserted into a new bulk/shell evolution. The prescribed Y=0 geometry would have to respond in a physical coupled experiment. Same-state-operator differences cancel local ambiguities but omit the reference vacuum contribution to absolute backreaction.

No new hot Big Bang mechanism, decay bath, observational prediction, global solution, or external mathematical novelty follows from these results.
