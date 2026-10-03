# Primary prospective active-source method

Fixed interval: eta=-4.5 to -3.5, saved unconstrained mode anchors u_1,w_1 and u_3,w_3. All two sources, both inherited grids, and K=64,128,256 are included. The density/pressure/action/contact convention, exact represented epsilon and pi, momentum nodes/weights, physical mass x=r=2, H=1, xi=0 and incoming state remain inherited.

No physical source, mode, stress, or real integral has been evaluated during this design.

The primary route uses native binary80 source jets and complex160 phases with exact phase variation-of-constants on each inherited time cell. Fixed GL24 and GL32 controls both evaluate the actual analytic bump jets directly; there is no source interpolant. Complex mode propagation is empirical binary80 arithmetic. Exact-ratio integer weighted accumulation uses MP80/100 coefficient seeds and p=266/333 fixed fractional bits, with no conversion of a physical input to binary64. MP precision comparisons diagnose accumulation/contact arithmetic; they do not enclose phase/source roundoff or quadrature error.

Write U=u/epsilon, W=w/epsilon and mu=weight*k^2/(2*pi^2). With C=sum(mu*Re(U)/(2k)), B=sum(mu*Re(W)/(2k)), E=sum(mu*k*Re(U)), J=sum(mu*Im(W)/2), separately expanded direct bare mode stresses are Rm=E+3L^2 C-L B-J and Pm=E/3-L^2 C-L B-J. G=3L^2 C-L B. Real forcing implies the canonical constant E-J is invariant and G'=L*(Rm-3Pm)+M*L*g, M=sum(mu/(2k)).

Canonical I is computed from the independent phase-aware Green endpoint as delta G minus M times a separately integrated L*g plus the directly integrated full metric/subtraction contact work. Measured stored endpoint density is never used to define I. Full R/P profiles and F=L*(R-3P)-3 h' L^4*(rho0+p0) are reconstructed independently and compared with saved history.

Signed flow projection uses the separately reconstructed final mode stress evaluated on measured u_3,w_3 minus the same operator evaluated on the predicted Green endpoint. Initial profile mismatch, direct contact integration, the weighted canonical invariant drift, and the mode-order control are reported separately. This distinguishes the inherited mode trajectory discrepancy from quadrature discrepancy. A finite-band analytic contact reduction is acceptable only after a symbolic action-matching proof; its original discrete momentum quadrature comparison is empirical.

Root must publish and verify the complete registration before enabling the physical path. Entire route budget remains 900 seconds / 262144 KiB for all twelve cases, both GL controls and both declared accumulation/contact precision levels.

## Final pre-freeze scope corrections

The direct contact integrands used in quadrature are evaluated in vectorized native binary80. Their exact represented values and rule weights are accumulated in MP80/100; profile and knot contacts are evaluated directly in MP80/100. Fixed GL24/32 controls use the exact same analytic source function at their own nodes. The entire all-case timer covers controls, both arithmetic levels, discrete action knot audits, signed projections, triangle envelopes, augmentation, and serialization.

Analytic contacts use M_A=K^2/(8pi^2), whereas the frozen finite momentum rule has M_d=sum(weight*k/(4pi^2)). The canonical raw analytic integral is I_A=delta G_d-M_d J+Q_contact, J=int Lg. The corrected primitive check includes (M_A-M_d)J. The matched physical target is I_matched=I_A+E_momentum, E_momentum=delta(C_R^d-C_R^A)+(M_d-M_A)J. The raw analytic I_A,D_cont_A,E_Q_A and all analytic profiles remain explicit. Principal I_ab,D_cont,E_Q refer to the matched target.

E_operator is the remaining stored-minus-direct discrete action residual at the two endpoints; all three knot residuals are separately reported so endpoint cancellation cannot hide a large operator mismatch. E_reconstruction=delta(E-J)+delta C_R^A+M_A J-Q_contact includes both invariant drift and contact quadrature/roundoff. Raw closure is D_cont_A=E_flow+E_momentum+E_operator+E_reconstruction, and matched closure is D_cont=E_flow+E_operator+E_reconstruction. Fixed midpoint flow, momentum-contact and operator pieces are separately reconstructed from retained u_2,w_2.

Numerical consistency, native GL controls, source/contact reconstruction, and cross-route comparisons use the fixed empirical 2e-7 gate. MP80/100 arithmetic gaps, serialization, and definitional closures use 1e-12 fatal arithmetic gates. The strict attribution trigger retains abs(D_S)>2e-6, abs(D_cont)<=0.1*abs(D_S), abs(E_Q)>=0.9*abs(D_S), applied to the matched target only after every fixed consistency check passes. The original metric results remain FAIL.
