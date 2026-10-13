# Independent review of the tensor-partner threshold

Reviewed 13 October 2026 UTC. **No critical correction is required in the current proof or checker.** The factorization, endpoint mappings, rational potential gate, and stated distinction between the dimensionless threshold and the physical gap are correct under the proof's assumptions. This review does not independently reconstruct the prior static-branch theorem or the full physical mode classification.

The exact reviewed-source hashes, checker commands and exits are pinned in `TENSOR_PARTNER_INDEPENDENT_RECEIPT.json`. The checker was copied to `tensor_partner_scratch/`; root sources and the earlier scalar replay artifacts were preserved.

## Factorization and boundary maps

With `d_z = rho d_y`, `psi = rho^(3/2) h`, and `S = (3/2)rho'`, the exact identities are

    D psi = rho^(5/2) h',
    H_+ = D*D = -d_z^2 + S^2 + S_z,
    H_- = DD* = -d_z^2 + S^2 - S_z.

The tensor equation maps to `H_+ psi = mu_T^2 psi`, and its norm becomes the ordinary `L^2(dz)` norm. For each fixed positive detuning, bounded `S` and `S_z` make the stated first-order Sobolev operators closed.

On `(-infinity,0]`, integration by parts pairs a free endpoint trace for `D` with a zero trace for `D*`, and vice versa. Consequently the even Neumann tensor sector maps to a Dirichlet partner; the odd Dirichlet tensor sector maps to the partner condition `chi'(0)+S_b chi(0)=0`. The odd partner form is exactly

    ||D* chi||^2 = int (|chi'|^2 + V_- |chi|^2) dz
                   + S_b |chi(0)|^2.

The endpoint term is positive when `S_b>0`. This sign is essential and is correct in the reviewed document. Continuity and the sum junction on identical glued interiors give the stated even Neumann and odd Dirichlet tensor sectors, conditional on the pure transverse-traceless reduction and absence of additional tensor shell sources.

The background substitutions give

    V_- = rho^2(3a^2/4 - 3a'/2)
        = 9/4 + rho^2(9p^2/16 - U/8).

Thus `U<=0` yields `H_->=9/4` in both partner domains. The polar decomposition compares the full nonzero spectra of the closed `D*D` and `DD*`, rather than only their square-integrable eigenfunctions. The even kernel is the single radial profile `psi proportional to rho^(3/2)`; the odd endpoint condition eliminates that kernel. The formal partner zero mode `rho^(-3/2)` fails square integrability at the cone.

## Threshold and fixed-detuning scope

At the cone, `S->3/2`, `S_z->0`, and `V_- - 9/4 = O(exp(2z))`. A bounded potential tending to zero at the infinite endpoint is a relatively compact perturbation of the half-line Laplacian. Endpoint changes have finite-rank resolvent differences. Hence the essential spectrum starts at `9/4`; the partner form excludes spectrum below that value, and polar mapping gives the same nonzero tensor threshold.

The statement is a dimensionless spectral gap from the even zero mode to `9/4`, with an odd lower threshold of `9/4`. At each fixed positive detuning, `H_b>0`, so it becomes the physical gap `m_T^2 >= (9/4)H_b^2` for nonzero tensor masses. **There is no uniform strictly positive physical mass bound as detuning tends to zero when `H_b` tends to zero.** The reviewed text states this distinction correctly.

Although `U<0` makes `V_->9/4` at every finite interior point, its excess tends to zero at the cone. This strict pointwise inequality cannot raise the spectral threshold above `9/4`. The reviewed text does not make that invalid upgrade. It also correctly limits the Higuchi comparison to this reduced tensor spectrum and does not claim complete absence of ghosts or full nonlinear stability.

## Exact polynomial gate and execution

Independent standard-library rational arithmetic, using the displayed coefficients and `r=1/100`, gives

    U(0) + sum_{j>=1} |U_j| r^j
      = -997874976899 / 13500000000000
      = approximately -0.07391666495548148 < 0.

This matches the stated rational value exactly and certifies negativity for every `eta in [-1/100,1/100]`. It is an interval upper bound, not a claim that the maximum equals that rational number. The extension to the registered detuning interval remains conditional on the prior theorem keeping the entire branch inside this scalar interval.

The unchanged copied checker passed under normal Python and `python -O`, each with exit 0, **8 identities/sign checks**, and **4 symbolic wrong-formula rejections**. The resulting JSON dictionaries agree exactly. Its failures use explicit exceptions rather than optimization-sensitive assertions. These executions validate the algebra and exact coefficient gate; operator domains, spectral completeness and physical mode classification remain mathematical/physical arguments, as the checker explicitly states.
