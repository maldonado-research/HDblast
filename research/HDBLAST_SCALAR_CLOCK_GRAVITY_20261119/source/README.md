The 32 tables contain the fixed-trajectory explicit order-y² Yukawa forcing on
the retained native quiet/plus/minus/far histories, on both meshes, for comoving
cutoffs Λ=60 and 120 and bare/static prescriptions. They include causal parent,
daughter, interaction and static-completion responses. They do not evolve the
causal moments on the perturbed wall.

Each `tables/{grid}_{arm}_{cutoff}_{scheme}.npz` exports `s,rho,p,j,drho,dp,dj,
rad,drad`, with derivatives taken with respect to proper time `s=H0*tau`.
All densities/currents are hatted physical coefficients per `lambda=y²`, with
native q weights already applied once. Add `lambda*(b/r_b)` times the triple at
the native matter interface. `rad` is only the daughter density; its pressure
is `rad/3`. `energy_comoving` and `daughter_comoving` are the respective physical
densities times `a³`. The suffixes `_parent`, `_daughter`, `_interaction` and
`_counterterm` give the corresponding components of all six source fields.

`SourceTable(path).evaluate(s)` returns `(rho,p,j,rad)` using cubic Hermite
interpolation with the exported analytic slopes. `evaluate(s,1)` returns their
proper-time derivatives. Convert derivatives to native coordinate time using
`ds/dT=exp(B_b)/r_b`. The source is exactly zero before turn-on and before the
fork. Evaluation after the full archived endpoint raises an error. The tables
continue after the readout `s=2.1` to each native archive endpoint, approximately
`2.1678–2.1685`.

One scalar function `F(phi)` is calibrated on **fine quiet** and used unchanged
on every grid/arm. Its inverse clock is a C2 cubic spline through 8193 quiet
points over `s=1.35–2.05`, with exact inverse endpoint slopes. The clipped
`sin⁴` gate and its first two field derivatives vanish outside the field
support. `Fphi` and `Fphiphi` are exported; the scalar current includes the
`Fphi` interaction and static-completion terms. Interaction pressure is zero
and static-completion pressure is minus one third of its density. Every
native background Hermite knot is included in the source integration/table
grid; switch events and the readout are also exact knots. Nominal RK4 steps
are 4096 and bath Gauss order is 64.

`SOURCES.json` records exact switch bounds (`gate_on,gate_off`), endpoint,
readout components and SHA256 for each case. `VERIFICATION.json` checks the
Ward identity, component sums, interpolation and independent complex-mode
spectral daughter energies. The fine-quiet spectral values at readout are
21.5778929151979 (Λ60) and 22.5105693261827 (Λ120); the causal tables differ by
less than 5.4e-10. The shared gate's sampled fine-quiet calibration error is
below 8.2e-13.

For an untouched archive, extract the exact tables with:

```sh
python source/archive.py extract --archive source/TABLES.zip --output source
```

To regenerate or verify, point `NATIVE` at the preceding covariance-feedback
capsule's `native/` directory containing `FORK.json` and `WALLS.zip`:

```sh
python source/source.py --native NATIVE --output TEMP/tables --jobs 4
python source/verify.py --native NATIVE --tables TEMP/tables --output TEMP/VERIFICATION.json
```

The scripts require Python, NumPy and SciPy. The native input hashes are
recorded in the manifest. No network or installed copy of the previous
solver is required. `source.py` packages its native archive reader;
`verify.py` independently solves complex free modes and integrates spectral
amplitudes as in the preceding clock study.
