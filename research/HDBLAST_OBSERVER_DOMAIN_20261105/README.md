# Conditional finite-cylinder observer-domain criteria

5 November 2026, America/Los_Angeles. Working notes, continuous arguments and
standard-library exact algebra checks. These are sufficient geometric and ODE
criteria, not an evaluated nonlinear amplitude range or a new field-equation
existence theorem.

The [previous checkpoint][budgets] bounds the actual-area shell quadrupole
moment only on a common admitted observer domain. This checkpoint explains
how to establish important parts of that domain from more primitive inputs,
while keeping the unevaluated model-specific assumptions explicit.

## Three complementary results

1. **A quantitative embedding and cut-metric test.** On the entire finite
   cylinder `[0,T_q] x S3`, embed the unit label sphere in `R4` and write
   `Psi_e=(-t_e,omega_e)`, `Psi_0=(q,n)`, with `q=-tau/rho_b`.
   A uniform bound `||D(Psi_e-Psi_0)||<=kappa<2/pi` gives global injectivity,
   immersion and compact embedding, including points on different proper-time
   cuts. The proof uses the product geodesic/chord comparison; it does not
   integrate along a chord outside the sphere. Projected velocity and
   label-Jacobi parameter bounds give an explicit sufficient test.
   [Proof](embedding/OBSERVER_DOMAIN.md).

2. **A stronger structural route for the actual normal observers.** Suppose
   the regular normal proper geodesic flow already exists, its cuts are
   spacelike, and the fixed reference coordinate vector `partial_t` is future
   timelike. Gauss orthogonality and a global angular covering argument make
   each cut a graph whose height satisfies

   ```
   Theta_q = -rho_b^2 / h_e(partial_t,X_q) < 0.
   ```

   Strict graph monotonicity proves injectivity across the whole finite
   cylinder. No Schur smallness condition is required for this physical route.
   Physical injectivity after composition also requires the retained shell
   identification to be an embedding. The independent Schur/degree estimate
   is retained as an alternative for maps without Gauss structure, alongside
   an explicitly nonphysical crossing example and inverse-derivative costs.
   [Proof](normal/OBSERVER_CLOCK_DOMAIN.md).

3. **Noncircular existence of the candidate flow.** A first-exit and compact
   continuation argument starts with coefficient bounds on a fixed regular
   neighborhood of the background trajectory. Its nonexit guard is

   ```
   a [exp(L S) v_disp + phi_L(S) P1] < r,
   phi_L(S) = (exp(L S)-1)/L,   phi_0(S)=S.
   ```

   It obtains the whole candidate interval before using parameter-variation
   bounds. Initial displacement and initial derivative bounds are kept
   distinct. Proper duration is `S=rho_b*T_q`; coefficient-only extensions,
   state/parameter segments, chart changes and the continuous-clock `C2`
   argument are charged explicitly. It does not infer injectivity from ODE
   uniqueness alone. [Proof](bootstrap/FLOW_BOOTSTRAP.md).

## Exact sufficient geometric margins

For the retained coarse central target, one convenient horizon is

```
ell=10^-6,   T_coord=13ell/8,   T_q=25ell/16,
T_coord-T_q=ell/16,
T_q-6145ell/4096=255ell/4096.
```

The metric estimate gives `1/2 I < G_e < 2 I` under its full-pullback
smallness inputs `b<=1/8`, `eta<=1/8` and the stated time bounds. These are
conditions to validate, not measured properties of the HDBLAST family.
A uniform velocity parameter bound giving `c=rho_b*a*U1<1/50` implies
`|t_e+q|<ell/32`, leaving a perturbed reference-slab margin greater than
`ell/32`; the background margin `ell/16` must not be confused with it.
The bootstrap, regular-field and common parameter hypotheses must all apply
to the same chosen preparation and domain.

This horizon covers the coarse central interval used by the earlier moment
zero criterion. It does not cover the later linear zero or the tail near
`|t|=2ell..3ell`. No new nonlinear moment zero or numerical amplitude window
is obtained merely from these symbolic inequalities.

## Reproduction and remaining work

See [REPRODUCE.md](REPRODUCE.md), [REFERENCES.md](REFERENCES.md),
[SOURCE_DEPENDENCIES.md](SOURCE_DEPENDENCIES.md) and
[NEXT_STEPS.md](NEXT_STEPS.md). The embedding checker authenticates its own
nested source capsule; its canonical reports stay outside that inner capsule.
The outer replay authenticates the complete checkpoint and runs the three
complementary components in normal and optimized Python modes.

The methods are established differential topology, Lorentzian Gauss calculus,
ODE continuation and exact inequality bookkeeping. Internal review and finite
algebra checks do not certify the retained PDE theory or external originality.
The actual regular coefficient neighborhood, initial-data/evolution constants,
parameter radius, pullback/modulus bounds and numerical remainder test remain
unevaluated. Physical preparation, radiation production, physical scales and
observational evidence for a higher-dimensional Big Bang remain unresolved.
Earlier frozen work is preserved. This incremental GitHub checkpoint includes
no new Zenodo edition or journal submission.

[budgets]: https://github.com/maldonado-research/HDblast/tree/d62bdda754507a29c198213946906cea6d9a9dbf/research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105
