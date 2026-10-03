# Independent source/operator formula review

Disposition: **PASS mathematical formula review and fabricated tests**, conditional on the new independently verified remote freeze and the root entry-point integrity guard. No registered source builder was called successfully. No physical bump values, retained scientific arrays or physical coefficients were evaluated.

Reviewed source pins:

- `source_models.py`: `475822f5b251dbcbdb72beddb85e49775f25520eb8d68ca77848272b79defdb8`.
- `generic_operator.py`: `9250165a7058441d79467ad6f16d56e3dabbbb559368d9600d2c2c155d87c6c1`.

## Source formulas

The registered builder uses `z=center+4`, hence the derivative scaling is exactly one. The formal denominator `(1-z²,-2z,-1)` is the Taylor denominator of `1-(z+x)²`. Its reciprocal recurrence is exact. The normalized exponential recurrence follows `(exp f)'=f' exp f`; its zero coefficient is one and the separate scalar factor is `exp(f0)`. For the allowed real centers, `f0` lies in `[-1/3,0]`, inside the scalar enclosure's required `[-1,0]` domain. An even degree-200 alternating sum is an upper bound, and adding the next negative term is a lower bound; the factorial terms decrease on this domain.

The positive profile needs normalized bump coefficients through 26 to obtain g through 24. The signed profile's coefficient `z beta[n]+beta[n-1]` is exactly the coefficient of `(z+x)B`. The derivative terms in gamma have the correct `(n-j+1)` and `(n+1)(n+2)` multipliers. The L series is the varying geometry `-1/(center+x)`, not a frozen midpoint value. Its square and `Lg` use correct coefficient convolutions.

Multiplying every normalized coefficient by the one scalar exponential enclosure is valid because all profile, differentiation and geometry operations are linear in that common scalar. Signed multipliers swap interval endpoints correctly. A chosen dyadic point may lie just outside its tiny scalar-product interval; `max(chosen-lo,hi-chosen)` still encloses every possible coefficient. Summing these coefficient radii times `half_width^n` is a uniform polynomial error bound, separate from the analytic tail.

The prior independent analytic proof supplies `|g_B|<=2737816576/62462907<44` and `|g_zB|<=2321190208/62462907<38` on each stated complex radius-1/8 disk. Therefore M=64 is conservative for g. Since `|L|<=8/27`, M=32 is conservative for Lg: even the larger bound gives `|Lg|<352/27<32`. At halfwidth 1/128, degree 24, their Cauchy tails are respectively `1/(15*2^90)` and `1/(15*2^91)`. These proofs concern the formulas and analytic domains, not sampled physical-source values.

## Generic operator and outward carry

The actual chosen polynomial is checked against `w'=i omega w-g`, even if a different solver frequency generated it. The code computes every coefficient of `W'-i omega W+g` and encloses its complex norm by the sum of the absolute real and imaginary parts. `U'=W` holds exactly, including divisions after coefficient quantization. For real omega the propagator has unit complex modulus. Thus a uniform source error rho and defect coefficients q_n contribute

`rw = incoming_rw + rho H + sum q_n H^(n+1)/(n+1)`,

`ru = incoming_ru + H incoming_rw + rho H²/2 + sum q_n H^(n+2)/((n+1)(n+2))`.

This formula includes all errors from the chosen coefficients; it has no division by omega and is stable at zero and tiny phase. The whole-interval fixture begins with exact zero states, carries both endpoint centers, adds their exact displacement under point quantization, and carries both proven norm errors. Its final negative states are therefore source/Duhamel operators over the whole prefix. For an arbitrary nonzero incoming state without a declared zero origin, the same outputs instead enclose negative endpoint states; the revised documentation and scope field correctly distinguish this.

`outward_radius(v,512)` is exactly `ceil(v*2^512)/2^512` for nonnegative rational v. It preserves zero and satisfies `v<=rounded<v+2^-512`. Its denominator is at most `2^512`, so repeatedly rounded prefix bounds avoid multiplication of unrelated exact source-error denominators. The rounding is applied after adding the exact endpoint center displacement; no error is dropped. Over 64 panels of width 1/64, pure radius-ceiling overhead is less than `64*2^-512` for w and `95.5*2^-512` for u. Including two-component point floors gives conservative additional overhead less than `192*2^-512` and `286.5*2^-512`, respectively. These terms are already included automatically in the carried radii; they are not substitutes for source or defect errors.

## Fabricated checks

- `test_fake_source_formal_series.py`: 38 checks, exact fractions, normal and optimized outputs byte-identical. Independently compares reciprocal recurrences with convolution, exponential recurrences with a finite power expansion, scalar alternating enclosures, dyadic points and Lg normalization. The actual geometry/g/Lg coefficient statements were extracted and evaluated only with artificial formal h and the nonregistered positive center 13/7. No authorization or registered source constructor was executed.
- Fresh `generic_operator.py --fabricated-only` receipts: 48 checks, four mutation controls, 1152 local cases and 18 whole-prefix cases, 2304 certificates, mode degree 96 and chosen coefficient precision 512 bits. Normal and optimized stable fields agree. Stable serialized moment SHA-256: `d281d6f0f57621ebc2cb2f09ae659c17d6223a1e8988da8c9b624313fcdc29f9`.
- `test_generic_rounding_carry_review.py`: 38 checks, normal and optimized outputs byte-identical. Intentionally coarse 16-bit chosen coefficients produce nonzero actual defects. Independently bounded entire kernels enclose the results with inherited initial-state and source uncertainty, at zero, tiny and high phase. Exact upward rounding, denominator bounds and invalid-radius rejection are also checked.

## Limits of this disposition

The source module's authorization routine checks the shape of a trusted root-provided attestation; it does not itself fetch or authenticate a remote freeze. That responsibility remains with the root integrity guard, and this review is conditional on its success. Nine rational frequency probes do not by themselves certify a continuous frequency interval. The generic mathematical argument is valid for each exact real frequency; the source approximation error is uniform in the real cell coordinate. Full direct pressure/metric-response certification and observational or nonlinear conclusions remain outside this checkpoint.

No defect was found in the reviewed formulas under these conditions. This review does not claim a registered numerical result or external mathematical novelty.
