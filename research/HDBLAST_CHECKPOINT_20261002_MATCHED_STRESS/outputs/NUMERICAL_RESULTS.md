# Registered matched stress and current calibration

The prescribed-geometry calibration passed both normal and optimized validators. At the twelve registered source/observation points, independent forced canonical modes with direct minimal stress agree with the logarithmic-memory and matched closed-stress route. Independent raw-mode and complete direct-stress histories support the finite-band trace and Ward checks.

This is a linear homogeneous scalar-to-stress/current response on fixed de Sitter geometry with an unchanged incoming BD state. It does not compute metric response, coupled shell dynamics, stability, particle yield, radiation transfer, thermalization, or heating.

## Registration and scales

H=1, r=2, a=-1/eta, epsilon=1e-4, and b=1 for the current translation. The compact sources s=a²delta_x are epsilon*B(eta+4) and epsilon*(eta+4)B(eta+4), supported on -5<eta<-3. Only eta={-5.5,-4.5,-4,-3.5,-2.5,-1.5} and K={64,128,256} enter the displayed response plots.

**Every numerical table below uses the registered normalization:** q=a²deltaQ/epsilon; q_prime and q_second differentiate a²deltaQ before division by epsilon; R=a⁴delta_rho/epsilon; P=a⁴delta_p/epsilon; J=a²delta_j/epsilon. Q0K alone is a physical baseline without epsilon division. Thus delta_rho=epsilon*R/a⁴, delta_p=epsilon*P/a⁴, and delta_j=epsilon*J/a². The changing conformal factors matter when comparing physical amplitudes at different times. `SUMMARY.json` records both normalized and physical continuum values.

## Independent agreement and uncertainty

| Quantity | Max direct-mode / primary finite-K gap | Max refinement gap | Max K=256 / continuum gap | Max analytic K=256 tail |
|---|---:|---:|---:|---:|
| q | 2.19008839e-17 | 1.13299127e-16 | 1.14468763e-07 | 2.34077946e-06 |
| q_prime | 3.98986399e-17 | 2.31065167e-15 | 3.49416208e-07 | 7.06641613e-05 |
| q_second | 3.74414041e-13 | 2.62358746e-11 | 4.24875653e-07 | 3.91015160e-03 |
| rho | 4.51028104e-17 | 3.51715185e-16 | 4.60093054e-08 | 2.51152401e-05 |
| p | 6.23993045e-14 | 4.37268906e-12 | 1.09617987e-07 | 6.75766826e-04 |
| current | 2.29850861e-17 | 1.15359111e-16 | 1.14468759e-07 | 2.34077946e-06 |

The fourth column is an observed difference, whereas the final column is a constructive analytic enclosure of the omitted combined momentum band. The derivative/stress tail bounds include the higher source derivatives and finite-band local-contact differences. They are not the earlier variance bound reused as a stress bound. Quadrature estimates, finite-step refinement and floating-point allowances remain numerical evidence; neither tiny observed agreement nor the UV theorem certifies the complete numerical result.

The maximum fine Ward-ledger endpoint residual is 1.985067613e-10; the maximum coarse/fine Ward-ledger difference is 1.269287688e-07. The maximum finite-band trace residual is 6.071532166e-17. These quantities use the registered normalization. The ledger independently integrates the saved direct-stress history; the density was not defined by that ledger. The finite-band checks retain Q0K and the finite-band local trace remainder.

The maximum linear Wronskian residual/epsilon is 5.611593276e-17, and the independently reconstructed physical residual/epsilon is 5.585139177e-17. Observation modes and direct-stress histories are archived; global all-time Wronskian maxima are recorded producer evidence.

## All registered stress/current values

| Source | eta | f=s/epsilon | R continuum | P continuum | J continuum | K=256 rho UV bound | K=256 p UV bound |
|---|---:|---:|---:|---:|---:|---:|---:|
| positive_B | -5.5 | 0 | 0 | 0 | 0 | 0.00000e+00 | 0.00000e+00 |
| positive_B | -4.5 | 0.71653131 | 0.0020042480543 | -0.021124739832 | 0.014263685337 | 4.01940e-06 | 3.29737e-04 |
| positive_B | -4 | 1 | 0.0042842889586 | 0.0015190139588 | 0.0027350803214 | 4.53416e-06 | 3.30360e-04 |
| positive_B | -3.5 | 0.71653131 | 0.0038084165019 | 0.0021106410509 | -0.013719856771 | 5.26041e-06 | 3.31167e-04 |
| positive_B | -2.5 | 0 | -0.0043929987994 | -0.0033768492009 | -0.011039893523 | 1.46946e-05 | 6.66012e-04 |
| positive_B | -1.5 | 0 | -0.0050684511809 | 0.00012862698832 | -0.0062778286939 | 2.51152e-05 | 6.75767e-04 |
| signed_uB | -5.5 | 0 | 0 | 0 | 0 | 0.00000e+00 | 0.00000e+00 |
| signed_uB | -4.5 | -0.35826566 | -0.0034560770229 | 0.0082167670363 | -0.0052522175109 | 3.80561e-06 | 3.13871e-04 |
| signed_uB | -4 | 0 | -0.0012085173869 | -0.012201246729 | 0.0076427855755 | 4.29530e-06 | 3.14395e-04 |
| signed_uB | -3.5 | 0.35826566 | 0.0042093071881 | -0.01374726327 | 0.0065460176075 | 5.00783e-06 | 3.15322e-04 |
| signed_uB | -2.5 | 0 | -0.00071298241828 | -0.0011744335324 | -0.0012742691331 | 1.38850e-05 | 6.33854e-04 |
| signed_uB | -1.5 | 0 | -0.0003881376203 | -9.9384502027e-05 | -0.00040900058374 | 2.37571e-05 | 6.43070e-04 |

Before the pulse the perturbative responses vanish. At observations after support, the source and all source jets are zero, so the local source contacts vanish while retarded coherent variance/stress terms can remain. A signed first-order stress perturbation is not a positive particle-energy yield. The Ward source at this order is Q0*delta_x_prime/2; the product deltaQ*delta_x_prime first appears at second order.

## Control accounting

The validator records 20 controls: 15 nonzero algebraic or saved-operator sensitivity witnesses and 5 explicitly rejected synthetic metadata/raw-format mutations. Of the sensitivity witnesses, 15 exceed their corresponding operational cross-route tolerance somewhere on the registered grid. The sensitivity witnesses are not all described as operational gate rejections. They reuse saved modes or exact contact algebra; they are not additional altered-state physical experiments.

| Control | Classification | Recorded outcome |
|---|---|---|
| missing_density_mass_operator | Saved-operator/algebraic sensitivity | max residual 4.1501157e+02; operational gate exceeded: True |
| missing_pressure_mass_operator | Saved-operator/algebraic sensitivity | max residual 4.1501157e+02; operational gate exceeded: True |
| missing_pressure_W4 | Saved-operator/algebraic sensitivity | max residual 6.1739934e-02; operational gate exceeded: True |
| improved_stress_density | Saved-operator/algebraic sensitivity | max residual 1.2572394e-02; operational gate exceeded: True |
| improved_stress_pressure | Saved-operator/algebraic sensitivity | max residual 5.1277432e-02; operational gate exceeded: True |
| omit_density_local_contact | Saved-operator/algebraic sensitivity | max residual 5.6932737e-04; operational gate exceeded: True |
| omit_pressure_local_contact | Saved-operator/algebraic sensitivity | max residual 1.7303416e-03; operational gate exceeded: True |
| wrong_density_contact_sign | Saved-operator/algebraic sensitivity | max residual 1.1386547e-03; operational gate exceeded: True |
| wrong_pressure_contact_sign | Saved-operator/algebraic sensitivity | max residual 3.4606833e-03; operational gate exceeded: True |
| omit_current_quadratic_mass_law_contact | Saved-operator/algebraic sensitivity | max residual 2.1108580e-03; operational gate exceeded: True |
| wrong_current_contact_sign | Saved-operator/algebraic sensitivity | max residual 4.2217160e-03; operational gate exceeded: True |
| continuum_Q0_at_finite_K | Saved-operator/algebraic sensitivity | max residual 3.7262582e-10; operational gate exceeded: True |
| continuum_anomaly_at_finite_K | Saved-operator/algebraic sensitivity | max residual 2.9900875e-07; operational gate exceeded: True |
| wrong_variance_derivative_sign | Saved-operator/algebraic sensitivity | max residual 9.7146323e-03; operational gate exceeded: True |
| missing_stress_conformal_factor | Saved-operator/algebraic sensitivity | max residual 4.2805380e-03; operational gate exceeded: True |
| moving_reference_scale_declaration | Metadata/raw-format mutation | Synthetic guard rejected |
| nonfinite_raw_value | Metadata/raw-format mutation | Synthetic guard rejected |
| insufficient_raw_precision | Metadata/raw-format mutation | Synthetic guard rejected |
| altered_initial_data | Metadata/raw-format mutation | Synthetic guard rejected |
| altered_raw_Wronskian | Metadata/raw-format mutation | Synthetic guard rejected |

## Reproducibility and presentation

Public prospective freeze: `4a5dad8dda6d0a57a9cf88cc8b4a9b8feeaea45d`. Full-registration SHA256: `4a1533fa7bf9df41cc635d76fa61a22b3aa119d8fbf215349fe09cf3ee54c893`. Primary runtime 1.273779 s; independent runtime 274.054955 s; independent peak RSS 102016 KiB. Timings describe this execution and are not portable performance guarantees.

`stress_response.svg/.png/.pdf` plots normalized density, pressure and current at registered times. `stress_uv_bounds.svg/.png/.pdf` separates observed K=256 discrepancies from analytic removed-band bounds. `stress_checks.svg/.png/.pdf` summarizes registered finite-K agreement and conservation residuals. Lines only guide the eye; no denser response history was evaluated for presentation. Zero values are omitted on logarithmic axes.

This report and renderer are explicitly post-run. No source, quadrature, mode, or response evaluation was added by them, and no frozen scientific source was modified. Full-precision values, all raw archive hashes and all frozen source hashes are retained in `SUMMARY.json`; figures have their own rendering provenance.

## Exact raw and presentation hashes

| Artifact | SHA256 |
|---|---|
| FULL_REGISTRATION.json | `4a1533fa7bf9df41cc635d76fa61a22b3aa119d8fbf215349fe09cf3ee54c893` |
| outputs/primary/results.json | `1d1d67ae723b46383d24e125177fc8d30c0c138fbed1347a43a9691628cf79f5` |
| outputs/independent/results.json | `1575915b6c361ff3efc4e3881760865a268db96d88fd77fe6d4f3f191f7a4bcf` |
| outputs/validation/CHECKS.json | `bebeebbd788498cac527db67fc19d907851e076006946b4f62289e9eacdfe0e6` |
| outputs/validation_optimized/CHECKS.json | `590d8bbc16f66a3ce175d887e456f3fb23de1b745255f51ced7372260abc07cf` |
| postrun/summarize_stress.py | `412508e4506d5454e906806271c0563279ccc70f48e1d075c6396337bf9c6857` |
| outputs/independent/stress_modes_positive_B_coarse.npz | `1b728240a92d481f035dca97487289e7b68176fdf77c1aecd9c01c4626821fc7` |
| outputs/independent/stress_modes_positive_B_fine.npz | `229be7333ee18273c035eb4a46e9ed96bcc9b39376243ead64a94e7c0ee9a729` |
| outputs/independent/stress_modes_signed_uB_coarse.npz | `9b90460430ca6580976a66051d561e0bbccfe8858077fc83ebafdcf1218fd576` |
| outputs/independent/stress_modes_signed_uB_fine.npz | `a25acf489745ce0801197db270069b7588ce91e5dae4e3107bafcfa9e631c34c` |
| outputs/primary/started.json | `1c2920cf229122adac95d9f8128694acb72990b5ce82ec12873b53c1271aeefd` |
| outputs/primary/partial_results.json | `57d8980edb42ef695716b7b2f05897f219e7171f64846311e514cb9103ca7360` |
| outputs/primary/EXECUTION.json | `a35b997acc6edc849555c9faa1e1cb816e42655c8c540c72cdb51b062144b133` |
