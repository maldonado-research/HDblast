Matching the zero-frequency value and slope makes the held-out bath response
converge with pair cutoff. It does not make the finite-cutoff responses equal.
For all 24 native momenta and the predeclared imaginary-frequency probes
`ν/q = 0.25, 0.5, 1, 2`, the largest remaining cutoff effect is 6.70% at
Λ=60, 1.73% at Λ=120, 0.436% at Λ=240 and 0.109% at Λ=480.

| Comoving pair cutoff Λ/H₀ | Largest missing continuum response | Native-weighted relative error, largest across ratios |
|---:|---:|---:|
| 60 | 6.6989% | 4.8530% |
| 120 | 1.7295% | 1.2443% |
| 240 | 0.43595% | 0.31308% |
| 480 | 0.10921% | 0.078397% |

The largest cutoff error occurs at the largest native q and `ν/q=2`.
At that probe (`q/H₀=7.98075`), raising Λ/H₀ from 60 to 480 changes the raw
response from 34.7304 to 2890.683 H₀² and the static-subtracted response from
−6.80283 to −19.95093 H₀². The value-and-slope matched response changes from
3.04606 to 3.26119 H₀², approaching the analytic continuum value 3.26476 H₀².
The weighted error is `Σw_q |R_Λ−R_∞| / Σw_q R_∞` evaluated separately at
each frequency ratio. It is a diagnostic of the response kernel across the
band, not an energy density or a response at one common frequency.

The experiment freezes `a=1` and `F=1`; retarded oscillator response is
independent of the prepared parent Gaussian state. The two spatial meshes
have exactly the same 24 q nodes and native weights. `NATIVE_MODES.json`
records those values and the originating fork hash. `PROTOCOL.json`
predeclares the probes, cutoffs, quadrature refinements and zero-frequency
matching conditions. No nonzero probe sets a matching parameter.

With retarded frequency `Ω=iν`, define `z=ν²>0`. The finite oscillator
Hamiltonian has precisely the previous spectral weights

```math
g_j^2=\frac{\omega_j(\omega_j^2-q^2)}{4\pi^2}\Delta\omega_j,
\qquad q<\omega_j<\Lambda,
\qquad \Sigma_{N,\Lambda}(z)=\sum_j\frac{g_j^2}{\omega_j^2+z}.
```

The matching coefficients are `C=Σg²/ω²` and `D=Σg²/ω⁴`, so
`Σ(0)=C` and `Σ′(0)=−D`. The twice-subtracted response is evaluated
under the sum to avoid cancellation of large terms:

```math
R_{N,\Lambda}(z)=\Sigma(z)-C+zD
=z^2\sum_j\frac{g_j^2}{\omega_j^4(\omega_j^2+z)}>0.
```

Its value and first derivative are exactly zero at `z=0` for every finite
bath. This common target is fixed before evaluating the four positive-frequency
probes. Adding any common constant and linear term to every matched response
would leave the reported cutoff differences unchanged. For a parent inverse
response at `a=1`, the associated bare parameter changes per `y²F²` are
`δm²=C` and `δZ=−D`: they give
`D_E=q²+m²+ν²−y²F²R(ν²)` after matching. The kinetic condition is an
additional matching choice beyond the prior static completion.

The continuous spectral integral gives

```math
\begin{aligned}
C_\Lambda&=\frac{\Lambda^2-q^2-q^2\ln(\Lambda^2/q^2)}{8\pi^2},\\
D_\Lambda&=\frac{\ln(\Lambda^2/q^2)+q^2/\Lambda^2-1}{8\pi^2},\\
R_\Lambda(z)&=\frac{(q^2+z)[\ln(1+z/q^2)-\ln(1+z/\Lambda^2)]-z+zq^2/\Lambda^2}{8\pi^2},\\
R_\infty(z)&=\frac{(q^2+z)\ln(1+z/q^2)-z}{8\pi^2}.
\end{aligned}
```

Only the first two polynomial terms diverge with cutoff. The held-out missing
tail is positive and has the dimensionally correct large-cutoff limit

```math
R_\infty-R_\Lambda
=\frac{(q^2+\nu^2)\ln(1+\nu^2/\Lambda^2)-\nu^2q^2/\Lambda^2}{8\pi^2}
=\frac{\nu^4}{8\pi^2\Lambda^2}
\left[1+O\!\left(\frac{q^2+\nu^2}{\Lambda^2}\right)\right].
```

Here `q,Λ,ν` are conformal energies, with the tables expressed in H₀ units.
`g²` has energy dimension four, `Σ,C,R` have dimension two, and `D` is
dimensionless. At frozen scale factor `a`, physical momenta and frequencies
are `p=H₀q/a`, `L=H₀Λ/a`, `ν_p=H₀ν/a`; the proper-time response is
`R_p=H₀²R/a²`. Thus its missing tail is `ν_p⁴/(8π²L²)`, an energy squared.
At `ν_p≪p`, `R_{p,∞}=ν_p⁴/(16π²p²)+O(ν_p⁶/p⁴)`.

Finite oscillator resolution is a separate error. The original Gauss rule
with N=64 reproduces the finite-cutoff response poorly at the smallest native
`q/H₀=0.0192511`: its worst relative discrepancies are 15.8%, 32.8%, 70.3%
and 89.8% as Λ increases from 60 to 480. Exact matching of two anchors in a
discrete bath does not resolve the modes near its threshold. Those errors
have very small native weight: the corresponding largest weighted relative
errors are 7.75e−11, 2.39e−9, 7.62e−8 and 2.43e−6. Refinement through N=2048
reduces the largest individual discrepancy to 8.51e−10. This test samples a
different kernel from the prior switched-yield calculation; its N=64 error
is specific to these low-frequency probes. An independent smooth-variable
integral and the separate Hamiltonian audit verify the continuous reference.

This establishes convergence of a finite, modewise matched response on a
fixed background. Matching the whole value and slope independently at each q
also fixes q-dependent finite parts. It does not construct a common local,
covariant counterterm action or its stress tensor, and does not establish
cutoff independence of the dynamical wall or of an assigned reheating budget.
All native nodes have q>0; a zero-momentum slope condition would need its own
infrared treatment.

Reproduce with NumPy and SciPy:

```sh
python response/response.py --output TEMP/response
```

`RESPONSE.json` contains summary errors and input hashes; `TABLES.npz`
contains every response and matching coefficient with explicit axes;
`MODES.csv` gives all 384 held-out mode/cutoff/frequency cases. Generate
the standalone figure with Matplotlib:

```sh
python response/plot.py --data TEMP/response/TABLES.npz --output TEMP/response.png
```
