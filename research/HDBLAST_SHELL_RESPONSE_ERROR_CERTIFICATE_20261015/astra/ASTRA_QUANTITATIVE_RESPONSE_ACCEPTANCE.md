# Independent acceptance of the quantitative shell-response enclosure

15 October 2026 UTC. This acceptance applies only to the exact source hashes
listed below. It is an internal mathematical and executable review, not an
external review or a formal verification of the retained analytic theorems.

**Decision:** accept the conditional statement
\(10^{15}<\mathcal C_{3\ell}<4\cdot10^{17}\) for the newly registered
\(\ell=10^{-10}\) pulse on the specified near-\(\phi=+1\),
\(\delta=1/1000\) branch, with the exact-data, gauge-completion, evolution,
and full-correction causal-margin hypotheses retained in the manuscript.
No substantive mathematical blocker was found in this continuation.

## Reviewed claims

The positive bulk-and-boundary shifted energy has the required shell source
sign \(+Gv\). The moving right edge contributes
\(-[(e_t-e_z)^2+\lambda e^2]/2\), so no unmentioned artificial boundary law
enters the estimate. The factor \(a/2\) in the norm differential inequality
is correct. The incoming characteristic identity controls the missing normal
derivative trace in \(L^2_t\) using bulk residuals and initial data, rather
than assuming that the displacement energy already supplies this trace.
The shell velocity/value contribution has exactly the factor
\(\sqrt{1/m+\gamma_*^2/\nu}\).

The smooth odd-extension free wave has the exact registered initial data,
zero shell value and acceleration, and normal derivative
\(q'(t)+e^{-t}q(t)\). Its true error sources are \(F=-Wu\) and
\(G=\widehat D\). The displayed bound \(|q'|\le16/\ell\) is valid.
The manuscript's \(2/\ell\) lower bound for the squared trial norm is
conservative and valid. The independent proof in
`RESIDUAL_TO_PHYSICAL_RESPONSE.md` also gives the stronger \(4/\ell\)
bound, but that improvement is not needed for the accepted interval.

The weight cancellation and normalization
\(\mathcal C_T=3\pi^2(\rho_b^2p_b^2)^{-1}\|w_z-(C_b/2)w\|_2^2\)
are correct. Coefficient uncertainty is included in the general residual
formula; in this comparator the shell value and acceleration vanish, which
justifies their absence from the specialized residual bound.

The new shell/collar argument uses the pinned analytic branch theorem's
Eq. (29) with outward rational constants at an admissible detuning.
Its first-exit inequalities and regular ODE continuation justify the short
collar coefficient box. I checked the source formulas and replayed the new
rational producer without running a historical numerical trajectory.
This review retains the analytic branch theorem as an input; it does not
claim a new full proof of that historical theorem.

The manuscript distinguishes the energy length from the correction-annulus
endpoint and supplies the strict geometry
\(L_{\rm corr}+c_*T_{\rm evol}=59\ell/8<8\ell\), with
\(T_1<T_{\rm evol}\). Its dependence on localization/inverse hypotheses
and an existence-smallness parameter interval is correctly retained.

## Replays and arithmetic result

`ASTRA_REPLAY_RECEIPT.json` records 12 executed runs: two successful runs
of the 17-identity independent checker; six expected failures for its three
injected formula errors in normal and optimized Python; two successful
25-group root-certificate runs; and two successful fresh background
producer runs. Normal and optimized root-certificate outputs agree byte
for byte. The freshly produced collar JSON agrees byte for byte with the
reviewed archive JSON. Copies and outputs were written only below `astra/`.

The exact rational certificate implies \(\Delta<8\),
\(\mathcal C_{T_1}>10^{15}\), and
\(\mathcal C_{T_1}<4\cdot10^{17}\). Approximate displays are
\(\Delta=7.081473609780481\), lower endpoint
\(1.7603897802771098\cdot10^{15}\), and upper endpoint
\(3.163172650869521\cdot10^{17}\); the rational JSON, not these displays,
is the acceptance basis.

## Source pins

| File | SHA-256 |
|---|---|
| `SHELL_RESPONSE_ERROR_CERTIFICATE.md` | `5f01be536c77123e9f667f4cd5514b52067641b2b389ba4e6eb655cb24654022` |
| `certify_response.py` | `f634507663644e204d65d12c996fce2b160fee36331d066f654ed5d44f37f386` |
| `check_response_certificate.py` | `953c8a51b21ed7c2cc3c51517a6df6d96ae0166c6f6a3990a6790fd68b7d3d36` |
| `archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md` | `0832a16c4bdf18cc188861a2c2860e7e6ff736dec39426d6ca7c589dc0434d5c` |
| `archive/derive_shell_box.py` | `3a4d2a71847114b8dce9f910679a1f4ce14898212001b0355279ec633ebed110` |
| `archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json` | `8ccba65a1c7fd66e175703cd4dfb8e1adb029b52d455a2b9a440f346bf2e4bc8` |
| `astra/RESIDUAL_TO_PHYSICAL_RESPONSE.md` | `6fd15b36168287b9283a1ad292ff71f2a77a190fa11a09324fd7fe6feeb2eb7b` |
| `astra/check_residual_response.py` | `d2eef3e0682a93db7676151d8f914d84b579f47c71ed26b18e8362f78c32ddf3` |

The main limitations are substantive and correctly stated: no resolved
waveform, no effective nonlinear \(\varepsilon\)-radius or remainder
bound, and no observational unit conversion follows. The large coefficient
is attached to a narrowly supported normalized mathematical tangent; it is
not a finite-amplitude physical prediction.
