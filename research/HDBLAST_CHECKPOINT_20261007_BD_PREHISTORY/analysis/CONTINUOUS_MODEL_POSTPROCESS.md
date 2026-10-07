# Continuous source representation from the registered coefficient models

This extension uses only the independent route's registered JSON output.
It requires no additional physical source, phase, mode, or saved-array
evaluation. The synthetic preparation does not evaluate those outputs.

Before its momentum loop, the independent route constructs one exact real
dyadic polynomial of degree at most 112 for each of the eleven parents and
each of the two sources. Its 113 coefficients come from source derivatives
through degree 114, and every child polynomial is an exact translation of
that same parent polynomial. Thus it defines one coherent source function
on the whole interior interval, independent of which momentum is later
probed. The frozen algorithm determines these polynomials; the JSON identifies
each complete ordered coefficient list by SHA256. These hashes are not
standalone exported coefficient vectors and do not by themselves prove
coefficient completeness. The reviewed source and highest-degree fabricated
tests establish that inventory.

For parent `j`, the builder computes an exact coefficient error norm `e_j`,
the already proved rational analytic tail `T_j=M_j*2^-112`, and exports

    coefficient_error_j = ceil_dyadic512(e_j),
    uniform_source_error_j = ceil_dyadic512(T_j+e_j).

Each upward rounding adds less than `q=2^-512`. Consequently

    T_j+coefficient_error_j−q < uniform_source_error_j
                              < T_j+coefficient_error_j+q.

Either exported expression bounds the true source error; the postprocessor
uses the conservative, transparent choice

    E_j = max(uniform_source_error_j, T_j+coefficient_error_j).

The postprocessor checks every parent endpoint, center, analytic radius,
halfwidth, child count, source majorant, analytic tail, degree, dyadic
precision, nonnegative radius, complete-polynomial hash encoding and source
identity. All 22 parent records must occur exactly once. It binds the budget
bytes to an externally pinned successful registered independent entry receipt
and checks the complete 22-event source journal. It does not fetch GitHub;
the supplied receipt pins retain the root authentication trust boundary.

For each source define `G_poly=sum_j width_j*E_j`. The exact Duhamel target
of this real polynomial family differs from the analytic BD target by the
flat cap plus the enclosed real source residual. Hence, for each
`K in {64,128,256}`, the existing proof gives

    C_K=C_cap(K)+G_poly,
    |delta W(a,k)|<=C_K,
    c=0,
    |delta A(k)|<=C_K/(2k),                 0<k<=K,
    |delta R_K(t)|<=C_K*(K²/252+K/294),
    |delta P_K(t)|<=C_K*(K³/162+5K/1764),  a<=t<=b.

The new contribution is that `G_poly` includes the actual registered real
coefficient choice, in addition to the analytic Taylor tail. The reference
on the polynomial side is its **exact** Duhamel integral at every momentum.
Its numerical phase/ODE evaluation is outside this certificate and remains
separately enclosed only at the nine registered probes. No stored binary80
mode is read, fitted, interpolated, normalized, or claimed accurate. The
imaginary constant invisible to these linear stresses and nonlinear
normalization beyond first order retain their existing limitations.

`certify_continuous_model.py --self-test --output FRESH.json` checks zero
coefficient-error metadata and independently upward-rounded coefficient
errors, monotonicity of the resulting bounds, a highest-degree nonzero
moment, and complete-polynomial hash changes under coefficient omission.
Its 53 checks reject 14 mutants, including incomplete parents, changed
degree/precision, a halved tail, changed analytic disks, negative or omitted
errors, incompatible rounding, malformed hashes, noncanonical rational
encodings, status promotion and a fabricated/physical scope mismatch.
Normal and optimized Python receipts agree exactly.

After the registered independent run, the CLI requires `--budget`,
`--budget-sha256`, `--entry-receipt`, `--entry-receipt-sha256`, and a fresh
`--output CONTINUOUS_MODEL_CERTIFICATE.json`. It computes exact rational
bounds and upward-rounded decimal displays. There is no new fitted
acceptance threshold. The inherited 2e-8 full-integral gate is unchanged.

The resulting finite-band representation certificate does not enclose
all-k numerical mode error, original binary80 state error, momentum
quadrature, full pressure/contact work, a UV tail, or a nonlinear physical
state. The full twelve-case certificate remains UNRESOLVED, metric
calibration FAIL, and a higher-dimensional Big Bang cause NOT_ESTABLISHED.
