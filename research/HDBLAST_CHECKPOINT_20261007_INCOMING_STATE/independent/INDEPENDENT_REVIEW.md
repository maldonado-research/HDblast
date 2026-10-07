# Independent exact incoming-state review

The algebra below was derived independently from the displayed ODE and direct
stress operators before the primary verifier or theorem was read. Afterwards,
the primary Laurent-polynomial implementation, theorem and result statement
were inspected. No mathematical discrepancy was found. This is internal
AI-assisted independent verification, not external peer review or a
proof-assistant formalization.

## Exact decomposition and normalization

For the same forcing and contacts, write the incoming errors as `u=deltaU(a)`
and `w=deltaW(a)`. For `k>0`, set

```
A = w/(2ik),  c = Re(u)-Im(w)/(2k),  E = exp(2ik(t-a)).
deltaW = 2ik A E,
deltaU = u-A+A E,
Re(deltaU) = c+Re(A E).
```

If `X=Re(AE)` and `Y=Im(AE)`, direct substitution, before any conservation
check, gives

```
deltaR = (k+3L²/(2k))c + (3L²/(2k))X + LY,
deltaP = (k/3-L²/(2k))c - (2k/3+L²/(2k))X + LY.
```

Thus the density's leading phase contribution cancels, while the pressure's
leading phase contribution persists. The imaginary constant in `u-A` is
invisible to these linear stress operators, rather than being an inferred
zero error.

In the fixed basis `q0=exp(-ikt)/sqrt(2k)`, take `q=q0(1+eps U)` with
`U'=W`. The exact Wronskian norm is

```
i(q* q'-q q*') = |1+eps U|² - Im[(1+eps U)* eps W]/k
              = 1 + eps(2 Re U-Im W/k)
                + eps²(|U|²-Im(U*W)/k).
```

Its first-order incoming difference is `2 eps c`. Exact first-order
normalization establishes `c=0`; arbitrary initial error, finite-precision
normalization, and nonlinear normalization alone do not establish that
linear condition without further analysis. The homogeneous ODE conserves
`2 Re(deltaU)-Im(deltaW)/k=2c` exactly.

## Stress and spectral hypotheses

The phase coefficients have exact squared norms

```
R_A² = L² + 9L⁴/(4k²),
P_A² = 4k²/9 + 5L²/3 + L⁴/(4k²).
```

The positive pressure majorant `2k/3+5L²/(4k)` has square exceeding `P_A²`
by `21L⁴/(16k²)`. For density, putting `q²=9L²/4` gives
`k sqrt(k²+q²) <= k²+q²/2`. Both inequalities preserve the singular
small-momentum terms rather than silently discarding them.

Uniform independently established `|c|<=chi`, `|A|<=sigma`, the measure
`k² dk/(2 Pi²)`, `Pi>=3`, and `0<L<=2/7` imply the following rational
all-time finite-band coefficient bounds:

| Contribution | Coefficient multiplying its uniform envelope |
|---|---|
| Density, `c` | `K⁴/72+K²/294` |
| Density, `A` | `K³/189+K/686` |
| Pressure, `c` | `K⁴/216+K²/882` |
| Pressure, `A` | `K⁴/108+5K²/1764` |

The density amplitude integral itself is exactly
`L[(K²+9L²/4)^(3/2)-(9L²/4)^(3/2)]/(6 Pi²)`. The rational density
coefficient is a valid upper bound on that integral. The pressure `c`
coefficient uses a triangle inequality, and all coefficients are
sensitivities, not achieved physical state-error bounds.

For nonuniform envelopes, a sufficient finite-band domination premise is
`integral_0^K (k³+k)(|c|+|A|) dk < infinity`. Near zero, independent
power envelopes `k^-p` satisfy that premise for `p<2`; `p=2` makes the
individual envelope integral diverge logarithmically. This is a sufficient
condition for the bounding method, rather than a necessary condition on
every correlated exact state. Bounded raw `w` generally gives `A=O(1/k)`,
which remains integrable here, but does not give a uniform `sigma`. At
infinity, appropriate separate moments or decay would need to be proved;
this finite-band result supplies no ultraviolet tail.

## Direct endpoint work and pressure-time transfer

Differentiating the original operators using `L'=L²` gives, with a general
complex forcing,

```
R' - L(R-3P) = (k Im g+L Re g)/(2k).
```

For incoming-state differences with `delta g=0`, this proves
`deltaR'=L(deltaR-3deltaP)`. With `L=-1/t`, it also gives
`(t deltaR)'=3 deltaP`. Under sufficient absolute domination,

```
integral_a^b L(deltaR_K-3deltaP_K) dt = deltaR_K(b)-deltaR_K(a),
integral_a^b deltaP_K dt = [b deltaR_K(b)-a deltaR_K(a)]/3.
```

Here `a=-9/2`, `b=-7/2`. The direct-work canonical coefficient is exactly
`64c/(1323k)` per mode: the `kc` sector cancels between endpoints. The
oscillatory coefficient is

```
[(6/(49k)-2i/7) exp(2ik) - (2/(27k)-2i/9)] A.
```

The rational integrated envelope coefficients are

| Quantity | `c` coefficient | `A` coefficient |
|---|---|---|
| Direct continuity work | `16K²/11907` | `16K³/1701+536K/250047` |
| Signed time-integrated pressure | `K⁴/216+K²/1134` | `K³/81+65K/23814` |

The smaller integrated-pressure growth is proved from its primitive. It is
not a bound on `integral |deltaP| dt`. Both endpoint amplitude columns use
valid triangle bounds and need not be attainable by one phase function.

Conservation is insufficient to identify the correct operators: removing
`3L²c/(2k)` from density and adding `L²c/(2k)` to pressure still preserves
the homogeneous Ward relation. The independent verifier demonstrates this
blind spot and rejects both altered expressions against the specified
direct operators.

## Verification and review receipts

`verify_exact_state.py` uses SymPy 1.14.0 and explicit exceptions throughout.
The pinned Python 3.12.14 interpreter ran it in ordinary and optimized modes:

```
/workspace/hdblast-research-work/environment/bin/python -B \
  independent/verify_exact_state.py --output /tmp/fresh-state-normal.json
/workspace/hdblast-research-work/environment/bin/python -B -O \
  independent/verify_exact_state.py --output /tmp/fresh-state-optimized.json
```

Run these commands from the checkpoint directory and choose fresh output
paths. Each mode passes **54 exact checks** and rejects **17 altered
formulas or invalid-premise controls**. The complete receipts agree after
excluding only `python_optimization`; both have scientific SHA-256
`f37114569eab847c17ac2afb6b2b888151763d81296eccefe30f2841f546051a`.
Their expressions are verified symbolically; exact rational witnesses are
used only to prove that mutations are incorrect.

`REVIEW.json` pins the reviewed files and records **24 exact rational
coefficient matches** against the separate standard-library implementation:
all eight density, pressure, continuity-work and pressure-time coefficients
at each of `K=64,128,256`. Its ordinary and optimized receipts also agree
except for the optimization flag. The reviewed theorem correctly states
the normalization, integrability and time-integrated-pressure limitations.

This work made **zero physical source evaluations, saved-array decodes or
physical trajectory evaluations**. It supplies no measured incoming-state
enclosure. The full twelve-case numerical pressure/contact certificate
remains **UNRESOLVED**; the historical metric result remains **FAIL**; a
higher-dimensional origin of the Big Bang remains **NOT_ESTABLISHED**.
