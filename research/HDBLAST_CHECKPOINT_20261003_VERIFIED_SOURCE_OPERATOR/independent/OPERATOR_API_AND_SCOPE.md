# Independent source/operator API and scope

This is a narrower prerequisite to the unexecuted direct-pressure design. It certifies a source polynomial enclosure and Duhamel operators at a fixed rational probe universe, not any retained trajectory or full twelve-case pressure ledger. No registered source callback has executed during preparation.

The settled candidate settings are 64panels, source Taylor degree24, exact alternating scalar-exponential degree200, dyadic chosen coefficients512bits, ODE candidate degree96, and the nine momenta `0,1/2^40,1/2^12,1/4,1,16,64,128,256`. The interval is `[-9/2,-7/2]`. Panel width is1/64 and half-width1/128. The resource budget is120seconds and131072KiB for a full route; complete fabricated construction/operator/work benchmarks justify engineering feasibility, not real-source success.

The implementation uses Python exact `fractions.Fraction` and integer dyadic point/radius operations. There are no third-party numeric dependencies, primary helper imports, registered source callbacks at import, retained arrays, momentum quadrature nodes, stored modes or endpoint-density inputs.

## Source construction

`source_models.registered_source_model(center, source, authorization, degree=24, bits=512)` returns:

    {
      "forcing": {
        "center_coefficients": tuple[Fraction],
        "coefficient_radii": tuple[Fraction],
        "uniform_source_error": Fraction,
        "analytic_tail": Fraction
      },
      "Lg": { same fields },
      "center": Fraction,
      "source": "positive_B" or "signed_uB",
      "half_width": Fraction(1,128),
      "degree": 24,
      "coefficient_bits": 512,
      ...
    }

Both chosen polynomials use the local **centered** coordinate, not the left coordinate. `uniform_source_error` includes analytic truncation, the exact rational scalar-exp enclosure, and every chosen coefficient error weighted by the half-width power. It is not just the Cauchy tail.

The independent model derives coefficients from `d=1-z²`, an exact formal reciprocal, a separately implemented normalized exponential recurrence, actual `h′,h″`, and the variable geometric series `L=-1/(center+x)`. Scalar exp is enclosed with exact alternating rational sums. The source forcing has independently reviewed disk boundM64 and tail `1/(15*2^90)`. The Lg model has boundM32 and tail `1/(15*2^91)`. `Lg` is real `L(t)*g(t)`, not an arbitrary homogeneous phase/work substitution.

`authorization` must be provided only after the production wrapper has independently authenticated and read back the entire new registration:

    {
      "status": "PASS_REMOTE_REGISTERED_SOURCE_GO",
      "freeze_commit": "<40 lowercase hex>",
      "registration_sha256": "<64 lowercase hex>",
      "input_scope": "NO_RETAINED_ARRAYS_SOURCE_OPERATOR_ONLY"
    }

This metadata check inside the builder is an additional execution guard; it does not fetch or authenticate remote bytes by itself. The root production entry-point must supply that independently verified result. The present fixtures deliberately exercise only rejected authorization, artificial formal expressions, and fake exponential-polynomial inputs.

## Local operator

Import the independent modules only:

    import generic_operator as op
    import source_models as sm
    from fractions import Fraction as Q

After the new registration passes the external guard, for each declared center:

    model = sm.registered_source_model(center, source, authenticated_receipt)
    force = model["forcing"]
    p_left = op.centered_to_left(force["center_coefficients"], Q(1,128))
    cert = op.defect_certificate(
        p_left, force["uniform_source_error"], Q(1,64), 2*k, 96,
        coefficient_bits=512)

The source target and operators are:

    M0(a,b)=integral_a^b g(s)ds,
    Mexp(a,b;k)=integral_a^b exp(2ik(b-s))g(s)ds,
    Mu(a,b;k)=integral_a^b phi1(2ik(b-s))*(b-s)*g(s)ds.

Mu is the stable entire kernel `(exp(2ik(b-s))-1)/(2ik)` including its limit `(b-s)` at k=0. No division by k occurs in the independent numerical implementation.

The certificate fields are `M0:(center,radius)`, `Mexp:((real,imag),radius)`, and `Duhamel_u:((real,imag),radius)`. Map `Duhamel_u` to the registered output name **Mu**. These radii are absolute source/state-defect norm bounds. `actual_residual_coefficients`, `chosen_w_coefficients` and `chosen_u_coefficients` are also returned (field names are `residual_coefficients`, `chosen_w_coefficients`, `chosen_u_coefficients`). The exact chosen polynomial defect is reconstructed against the true declared omega, rather than trusting the recurrence. The source polynomial and all residual coefficients are exact rational objects.

`op.output_moments(cert)` emits exact rational endpoints with integer-string `numerator`/`denominator` fields. A common wrapper may format those into its common real/imag `{lo,hi}` schema, keeping every endpoint exact. Imaginary parts of M0 and Lg are exactlyzero. At k=0 the real-source Mexp and Mu imaginary parts are exactlyzero by the real ODE identity, rather than a numerical cancellation.

The actual production `route.py` applies one additional proved output operation: it floors every lower endpoint and ceils every upper endpoint to the512bit dyadic grid. This preserves inclusion and avoids large exact denominator strings. Each coordinate half-width increases by less than2^-512. The registered `total_absolute_radii` values are the **L1 sum of the actual exported rectangle half-widths**, computed from those endpoint bytes. Thus a two-component rectangle derived from a complex norm radiusR normally reports about2R, while a proved real target reports aboutR. Core disk radii and exported rectangle radii are distinct quantities. The core disk radius remains the one propagated in the ODE proof; the larger exported L1 radius is used for the common width gate.

The separate method evidence represents nonnegative source/residual/point-displacement/rounding components as independently upward-rounded512bit bounds. It explicitly identifies those values as upper bounds rather than exact signed components. No negative signed scientific residual is rounded with that nonnegative helper.

## Whole interval

Start each source/probe prefix with exact zero complex incoming modes and zero norm radii. Carry modes in the same direction as the original source operator:

    whole = op.defect_certificate(
        p_left, source_error, Q(1,64), 2*k, 96, coefficient_bits=512,
        incoming_w=w, incoming_u=u,
        incoming_w_error=rw, incoming_u_error=ru)
    w_exact_endpoint = op.scale(whole["Mexp"][0], Q(-1))
    u_exact_endpoint = op.scale(whole["Duhamel_u"][0], Q(-1))
    w_next = op.chosen_complex(w_exact_endpoint, 512)
    u_next = op.chosen_complex(u_exact_endpoint, 512)
    rw_next = op.outward_radius(
        whole["Mexp"][1] + op.norm_upper(
            op.add(w_exact_endpoint, op.scale(w_next, Q(-1)))))
    ru_next = op.outward_radius(
        whole["Duhamel_u"][1] + op.norm_upper(
            op.add(u_exact_endpoint, op.scale(u_next, Q(-1)))))

At the final panel, the whole-interval Mexp/Mu centers are `-w_next,-u_next`, with radii `rw_next,ru_next`. This is a zero-origin source prefix. For arbitrary nonzero incoming modes, the same core instead encloses a general negative endpoint state; it must not be mislabeled a source-only operator. That limitation is also marked in `incoming_state_scope`.

For M0, add each exact local polynomial integral center. For Lg, translate its own model and call `op.exact_integral()`; local integral radius is `H*uniform_source_error`. Add cumulative radii using `op.outward_radius(old+new)` after every panel. Exact rational radius sums across varying exponential enclosures otherwise grow into large denominator products. The chosen512bit upward radius operation adds at most2^-512 each step, with an explicit exact test; it preserves all prior errors and avoids that denominator growth.

The state center floor chooses a precise dyadic point. Its exact displacement is always added to the norm radius before the upward radius operation. Do not assign512bit accuracy to a state without that displacement/error term. All source/model, polynomial-defect and cumulative rounding components should be retained in separate error-budget/proof receipts.

## Output universe

The full wrapper must output1152panel moment rows,18whole moment rows,128panel source-work rows and2whole source-work rows. Preserve root's common row/schema order and metadata. Labels should use canonical reduced exact momenta, fixed panel membership and fixed interval endpoints. It must never substitute a favorable subset or compare only midpoints.

`fabricated_complete_route.py` benchmarks own fake model generation, all1152local and18prefix rows, plus128+2fake geometry-work integrals and serialization. Its sources are explicitly finite artificial polynomials times exp(q) with declared uniform source uncertainty. They are not registered B/zB callbacks, hidden decoded arrays or physical evaluations. The source-formula AST audit separately checks the actual derivative/geometry code using fabricated jets.

The first complete fabricated fixture exposed a Python4300decimal-digit serialization limit caused by exact prefix-radius denominator products. That pre-freeze execution failure is preserved separately; the correction was explicit upward dyadic radius rounding. No physical output motivated the change. Neither raising the string limit nor dropping error terms was used.
