# D1: finite-band quantum modes and stress/current exchange

Registered 1 October 2026, before any quantum-mode integration in this round. Input archives and classical crossing summaries have been inspected; earlier Gaussian estimates are known. This disclosure makes no novelty claim and does not preregister discovery of a hot Big Bang.

## Fixed model and data

Use the source-free (Y=0) tuned A1 background: delta=0.1, c=0.5975949350280132, d=-3.106933495673783, rb=7.838285538073243. This remains a model change from the original delta=0.001 hypothesis. Pinned archive input commit: 8f67197b730d4e1c43554b86f224c29cc72629eb. Input SHA-256 values are in INPUTS.json.

Use main_dstar_Y0_dc1e-2_dzf5e-4 as primary and dzf1e-3 as a background-resolution comparison. Both have 305 saved shell records; they later stop non-finite. Our common final time s=H0*tau=6.9 precedes their last finite near-shell constraint checks (~6.9603 and ~6.9596). Constraints are sampled intermittently, not certified between records.

Candidate spectator action on the induced spatially flat FRW shell: one real minimal scalar, xi=0, m_hat²=G²(phi-phi_star)², G=100, phi_star=0.5, m0=0. No decay, new matter backreaction or new five-dimensional evolution is included. s=H0*tau, G=gbar/H0, a=exp(ln_a), a(0)=1. Stress and current are in H0^4 units.

## Geometry, state and momentum

Primary reconstruction: CubicHermiteSpline for ln_a with archived H/H0 derivatives, and for phi with archived v/H0 derivatives. Derive all evaluated H and phi_dot from these same interpolants. A CubicSpline reconstruction supplies a separate control. Their finite smoothness is explicitly insufficient for a blanket Hadamard/high-order adiabatic claim.

Locate the primary phi=0.5 crossing using its interpolant. Set q=G*abs(phi_dot_at_crossing), k=a_cross*sqrt(q)*kappa. Primary fixed comoving Gauss-Legendre band kappa in [0.05,3], 32 nodes. Freeze this physical band from the primary history for comparisons.

At initial s=0 specify
chi=1/sqrt(2*a^3*omega),
chi_dot=(-3H/2-omega_dot/(2omega)-i*omega)*chi,
omega²=k²/a²+m_hat²,
omega_dot=(-H*k²/a²+G²*(phi-phi_star)*phi_dot)/omega.
This finite-time WKB-amplitude state is explicit finite-band data, not a proved infinite-order adiabatic vacuum.

Specify the exact out-reference at s=6.9 by physical Hamiltonian data with the same amplitude and chi_dot=-i*omega*chi; evolve it backward under the same operator. The in-minus-out stress/current is a finite relative diagnostic. At the endpoint its excess energy equals the finite-band integral omega*|beta|²/a³. Local state-independent counterterms cancel for a legitimate same-theory state pair; this does not determine the absolute reference stress/current. Common Hadamard completion is a conditional construction on sufficiently smooth backgrounds, not a numerical claim about our spline.

## Registered numerical controls

Primary DOP853 rtol=1e-9, atol=1e-11; tighter rtol=1e-11, atol=1e-13. Use a conservative phase step and explicit treatment of interpolation knots where canonical frequencies contain derivative jumps.

Required controls:
1. Maximum Wronskian drift <=1e-6.
2. Independently integrated finite-band comoving energy ledger normalized residual <=1e-5, for raw states and state difference. Normalize by a sum of absolute energy/work terms with a stated floor, never a cancelling net work alone.
3. Mapped-clock occupation absolute difference <=1e-5 on the declared checked modes. Transform both initial amplitudes and derivatives; do not choose separate clock vacua.
4. Exact free/conformal geometry control, with occupation/Wronskian absolute errors <=1e-6. A conformal null occupation is not a zero renormalized stress claim.
5. Report exact-pair Bogoliubov coefficient drift, endpoint energy/occupation agreement, and explicit current/pressure coherence terms. Their computation is a diagnostic, not a substitute for the independent integrated ledger.

Prespecified sensitivity variants: 64 momentum nodes; kappa_max=2.5 and 3.5; alternate C2 spline; coarse background; initial physical-Hamiltonian state; initial slice s=0.2; out endpoint s=5.0; tighter mode tolerance. Changed endpoints/initializations/bands change the declared state comparison. Report all variants rather than enforcing a discovery threshold. A new infrared lower-edge variant, if added, must be labeled exploratory.

Compare the endpoint spectrum with the earlier exp(-pi*kappa²) single-crossing ansatz descriptively. Quantify mode integration, quadrature, reconstruction, background and state sensitivity separately. Report the actual peak physical momentum and mass on the time interval; the earlier characteristic cutoff screen is not automatic UV certification for this new band.

## Exact ledger and scope

For each minimal scalar mode:
e=(|chi_dot|²+omega²|chi|²)/2,
p=(|chi_dot|²-(k²/(3a²)+m²)|chi|²)/2,
j=G²*(phi-phi_star)|chi|².
For fixed comoving k weights:
rho_dot+3H(rho+p)=phi_dot*j.
Integrate a³rho with scalar work and pressure work, retaining both. The identical law holds for the difference of the two exact states. A moving physical cutoff needs an explicit mode-boundary flux.

Never replace the absolute five-dimensional source by this difference without treating the reference polarization/current. Occupation does not establish a thermal bath. Passing a finite-band Ward identity cannot certify renormalization or UV validity. Preserve failed controls, registration deviations and the original September/first-October verdicts.
