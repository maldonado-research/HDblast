# B1 matter work, tension, and the signed five-dimensional interface flux

The B1 particle source has a consistent local energy bookkeeping interpretation, but positive matter work does not imply net depletion of a bulk reservoir. For the selected archived pair on x=H0 tau in [1,4], the inferred two-sided bulk-to-shell flux is **negative**, approximately -75.57 and -75.59 in common code energy units. Positive matter work is only about 0.01473. The tension contribution is opposite in sign and about 5131 times larger. An independent integration-by-parts calculation from the saved tension, expansion, and endpoint energies confirms this sign and magnitude.

This is a finite numerical comparison of a specified historical model and a local boundary identity. The flux is inferred using the imposed scalar junction. No separately measured decrease of a bulk energy integral, physical reheating success, or globally conserved positive energy is established.

## 1. Action, orientation, and the two-copy factor

Use the published 22 September action, signature (-++++), two identical interiors Z<0, and **one** shell action:

\[
S=\sum_{a=1}^2\left\{\frac1{2\kappa_5^2}\int_{M_a}\sqrt{-g}
[R-(\nabla\phi)^2-2U]+\frac1{\kappa_5^2}\int_\Sigma\sqrt{-h}K_a\right\}
+\int_\Sigma\sqrt{-h}[-\sigma(\phi)/\kappa_5^2+\mathcal L_m].
\]

The outward normal of each retained interior is n=exp(-B_b) partial_Z, directed toward the shell. The shell proper-time vector is u=exp(-B_b) partial_T. Put v=u(phi), w=n(phi), lambda=sigma/kappa5^2, and j=-delta S_m/[sqrt(-h) delta phi]. Variation gives

\[
K_{ab}-Kh_{ab}=\frac{\kappa_5^2}{2}[-\lambda h_{ab}+T_{ab}],
\qquad w=-\frac{\sigma_{,\phi}+\kappa_5^2j}{2}.
\]

The scalar coefficient in the variation is -2w/kappa5^2-sigma_phi/kappa5^2-j. The factor two comes from the two bulk pieces; it does not double the shell particles. Trace reversal for T^a_b=diag(-rho,p,p,p) gives

\[
k_s=n(A)=\frac{\sigma+\kappa_5^2\rho}{6},\qquad
k_0=n(B)=\frac{\sigma-\kappa_5^2(2\rho+3p)}6.
\]

Reversing the normal requires reversing K and w together. With the retained normal, the scalar stress has T_bulk(n,u)=wv/kappa5^2. The physical outward energy current of one interior is -T_bulk(n,u), so the two-copy flux **into the shell** is

\[
\mathcal F_\Sigma=-2T_{\rm bulk}(n,u)
=-\frac{2wv}{\kappa_5^2}
=\left(\frac{\sigma_{,\phi}}{\kappa_5^2}+j\right)v.
\tag{1}
\]

Matter diffeomorphism invariance supplies dot(rho)+3H(rho+p)=jv. Adding dot(lambda) yields the total shell budget

\[
\dot{(\lambda+\rho)}+3H(\rho+p)=\mathcal F_\Sigma.
\tag{2}
\]

Equivalently, Codazzi is dot(k_s)+H(k_s-k_0)=-wv/3. Direct substitution in the conformal PDE momentum constraint gives

\[
\frac{M_b}{e^{2B_b}}=-\frac{\kappa_5^2}{2}
[\dot\rho+3H(\rho+p)-jv].
\]

These identities independently fix the signs and factors. They also explain why evaluating the boundary momentum row with the same differentiated junction and matter RHS cannot constitute an independent measurement of the donor flux.

## 2. What the B1 source pays for

For a particle node, omega_i^2=p_i^2+m_chi^2, m_chi^2=m0^2+G^2(phi-phi_star)^2, and

\[
\rho_\chi=\sum_i n_i\omega_i,\quad
p_\chi=\sum_i\frac{n_i p_i^2}{3\omega_i},\quad
j_{\rm ad}=G^2(\phi-\phi_*)\sum_i\frac{n_i}{\omega_i}.
\]

Redshift gives dot(p_i)=-Hp_i. Time-dilated decay uses dot(N_i)|decay=-Gamma m_chi N_i/omega_i and Q=Gamma m_chi sum n_i. The resulting energy equations are

\[
\dot\rho_\chi+3H(\rho_\chi+p_\chi)
=j_{\rm ad}v+P_{\rm prod}-Q,
\qquad \dot\rho_r+4H\rho_r=Q.
\tag{3}
\]

B1 fills each cohort with a prescribed field-space ramp. Its bounded force j_prod satisfies j_prod v=P_prod, including the implemented stall switch. Thus j=j_ad+j_prod enters the scalar boundary condition. The adiabatic term pays for changing particle mass; the production term pays for changing particle number. They are distinct terms in the derivative of sum n_i omega_i. Q is an internal chi-to-radiation transfer and cancels when the two matter sectors are added.

The curved-crossing correction D is used once in the target population, N_fin proportional to (1+D/q). The same N_fin enters both the added particles and the ramp work. There is no separate D energy deposit to add to (3). Adding one would double count this production prescription. Likewise, tension is included once through sigma and its derivative; it is not a second copy of the particle energy.

This is the particle/ramp approximation actually implemented. The finite production force is not a derivation of the full renormalized, coherent quantum stress/current during a nonadiabatic event. In a calculation with exact quantum modes, production work is already in the matched Ward identity; adding this ramp force would double count it. A local quantum potential absorbed into sigma must also be removed from an explicit matter potential/source. Conservation alone cannot exclude a consistently doubled action term; the operator and matching specification matters.

There is a small finite-step exception to perfectly paid ramp initiation. `add_cohort` starts at the post-crossing progress u0 and inserts N_fin f(u0), rather than starting at exactly zero. The solver records this unpaid jump. For the selected coarse/fine runs the inserted code densities are 3.669506031898071e-13 and 7.250085330163294e-14. They are tiny, but not zero. The separately recorded instant-insertion control is a much different, deliberately unpaid prescription. No continuum exactness is inferred from these small numerical values.

## 3. One normalization for both energies

The selected pair is `main_ps0.5_G100_y1_bcut_dzf{1e-3,5e-4}`: delta=0.1, d=-3.106933495673783, G/H0=100, phi_star=0.5, y=1, b=kappa5^2 H0^3=6.49e-6, and rho_b=H0^-1=7.838285538073243. Both use the same normalization a=exp(A_b-A_b(0)), a(0)=1, and x=H0 tau. These are not B2's delta=0.001 runs.

The saved code variables have

\[
R=\kappa_5^2\rho_r,\quad X=\kappa_5^2\rho_\chi,
\quad P=\kappa_5^2(p_\chi+\rho_r/3),
\quad J=\kappa_5^2(j_{\rm ad}+j_{\rm prod}),
\quad \widehat v=d\phi/dx,\quad h=H/H_0.
\]

`p_hat` is chi pressure. **`P_hat` is production power**, defined in the producer as `jprod_hat*v_over_H0`; it is not pressure or an independent measurement. `R_T` and `rho_T` are direct RHS diagnostics. The independent endpoint comparisons below do not use those three columns.

The two comoving energies, with the same volume and common kappa5-squared multiplier, are E_lambda=a^3 sigma and E_m=a^3(X+R). Their work budgets are

\[
I_\lambda=\int_1^4a^3\sigma_{,\phi}\widehat v\,dx
=\Delta E_\lambda-\int_1^4 3h a^3\sigma\,dx,
\tag{4}
\]
\[
I_m=\int_1^4a^3J\widehat v\,dx
=\Delta E_m+\int_1^4 3h a^3P\,dx.
\tag{5}
\]

The minus sign in (4) is vacuum pressure p_lambda=-lambda. Consequently the tension's **comoving stored energy can grow** even when its pressure-corrected work is negative. Dropping the pressure term would compare different budgets.

## 4. Numerical result and the independent form

The seven-column tables retain all saved interior nodes on [1,4], with separate linear interpolation of each weighted column at the endpoints. Each contains 155 rows. Composite trapezoidal quadrature is used consistently. The explicit tension law sigma=2(1-phi+phi^3/3)+delta(1+c phi+d phi^2/2) and its derivative reproduce the saved sigma/sigma1 within rounding error, independently of the producer import.

| Quantity | Coarse | Fine |
|---|---:|---:|
| Delta E_m | 0.0117979693500 | 0.0118018861525 |
| Matter pressure work | 0.00293270021881 | 0.00293351239149 |
| I_m, scalar-source form | 0.0147310252609 | 0.0147357339706 |
| I_m, endpoint plus pressure | 0.0147306695688 | 0.0147353985440 |
| Delta E_lambda | 132.009086720 | 132.055558213 |
| Integral 3h a^3 sigma | 207.591139365 | 207.663410280 |
| I_lambda, scalar-force form | -75.5806502653 | -75.6065602970 |
| I_lambda, endpoint minus expansion | -75.5820526450 | -75.6078520668 |
| Inferred two-copy bulk-to-shell flux I_lambda+I_m | **-75.5659192400** | **-75.5918245630** |
| Total endpoint plus pressure form | -75.5673219754 | -75.5931166683 |

The relative differences between the independent tension forms are 1.86e-5 and 1.71e-5; the matter forms differ by 2.41e-5 and 2.28e-5. These are finite saved-data discrepancies, not error bounds. Retaining every other saved node still gives negative net flux, -75.56223 and -75.58811. Its magnitude is far larger than either numerical discrepancy or the positive matter work.

Omitting tension would report the positive numbers 0.014731/0.014736 as the bulk-to-shell flux and would reverse the inferred net direction. Missing one of the two interiors would instead halve the negative flux. With the correct orientation and complete action, the selected interval shows **net local energy transfer from the shell into the interiors**, while matter simultaneously gains positive scalar-coupling work. This does not identify a unique global partition of the energy into gravity, scalar kinetic energy, and potential energy.

## 5. What the saved trace cannot establish

The full-cadence records contain matter densities, pressures, scalar velocity and source, tension, and sparse normalized constraint diagnostics. They do not contain an independently evaluated w=n(phi), boundary scalar flux, bulk energy integral, or outer-boundary flux history. Computing w=-(sigma1+J)/2 from those same records simply reapplies the boundary condition.

The NPZs additionally contain A,B,phi,A_T snapshots at approximately computational T=0,1,2,3 on every fourth spatial node. They lack bulk phi_T and B_T and do not resolve the full production interval in time. One could form uncontrolled finite differences from these sparse snapshots, but they are insufficient for the requested independently controlled bulk energy balance. The recorded Weyl scalar is itself calculated algebraically from H, scalar, and source terms; it is not a bulk-energy reservoir observable. Sparse momentum/Hamiltonian residuals are useful solver diagnostics, not a replacement for that missing budget.

More fundamentally, a bulk reservoir requires a specified energy current. For a chosen time vector xi, define J_bulk^A=-T_bulk^A_B xi^B. On the scalar equation,

\[
\nabla_A J_{\rm bulk}^A=-T_{\rm bulk}^{AB}\nabla_{(A}\xi_{B)}.
\tag{6}
\]

For xi=u at the shell, the boundary flux is the one in (1). Unless xi extends as a Killing field, an integrated bulk budget also has a geometry-work term, as well as any other boundary flux. In the evolving conformal chart neither is automatically zero. The scalar potential need not provide a positive density either. Hence a local jv term cannot by itself be equated to the negative change of a positive global bulk reservoir.

## 6. Evidence, replay, and retained limits

The controlling source files match the public commit `e1d1285c65675899b9a1e3958171c11a8fbab961`; `SOURCES.json` gives immutable URLs and hashes. Unlike the former B2 source-version issue, all three B1 producer/helper SHA256 values recorded in the selected summaries match those public files exactly. The optional extractor uses the authenticated original NPZ files. The compact derived tables suffice for the public replay.

Run without third-party packages:

```sh
python replay_signed_budget.py --output SIGNED_BUDGET.json
```

`check_boundary_algebra.py` provides an optional independent SymPy check of the displayed finite identities and wrong-sign/factor controls. These finite checks are not analytic existence proofs or evidence that a microscopic reheating model is correct. The old B1 registered outcome remains INCONCLUSIVE; this note does not promote its descriptive results to a successful radiation-era claim.

One inherited code limitation remains explicit: `hat_rates` omits the time derivative of j_prod, so the returned differentiated scalar Neumann helper is incomplete during production. The current constraint calculation consumes the differentiated A boundary row, not that scalar helper; this is not evidence that the current evolution or its matter work is wrong. Any future use of the differentiated scalar row must include the missing derivative.

Primary retained sources: the [original action and matter extension](https://github.com/maldonado-research/HDblast/blob/e1d1285c65675899b9a1e3958171c11a8fbab961/hdblast/checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md), [B1 matter adapter](https://github.com/maldonado-research/HDblast/blob/e1d1285c65675899b9a1e3958171c11a8fbab961/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/shell_matter.py), [coupled solver](https://github.com/maldonado-research/HDblast/blob/e1d1285c65675899b9a1e3958171c11a8fbab961/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/evolve_b1.py), and the [later action/junction audit](https://github.com/maldonado-research/HDblast/blob/e1d1285c65675899b9a1e3958171c11a8fbab961/research/HDBLAST_CHECKPOINT_20261002_JUNCTIONS/ACTION_AND_JUNCTION_SPECIFICATION.md).
