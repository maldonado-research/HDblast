# Two field patterns at the same point on a spherical shell

Charged daughters born at different times experience different portions of the incident waves. A pair of fields can therefore change the distribution of sideways momenta. This experiment follows the two charges of each decay along the unit three-sphere with the measured electric and magnetic histories on the expanding shell. The question is whether both directions perpendicular to the dominant particle stream develop the pressure of a drifting thermal gas.

## What the comparisons show

* Two quadratures of the same early pulse at angular wavenumber 7 make the pressure averaged over the shell look thermal in both lateral directions. On independent birth samples the two pooled principal pressures are 0.995 and 1.005 times their targets for either interaction history. Rotating the measurement axes at each observation point reveals a weaker pressure below 0.000016 of target in the dedicated case and 0.000023 in the shared case. A single rotated mode produces these two spatially varying quadratures. An observer who pooled the two fixed axes would miss the imbalance. The mean energy ratios are 1.022 and 1.007, respectively, and the discrepancies in the cumulative distribution of lateral momentum are 36.80 and 36.94 percentage points.
* Two genuinely different time profiles can fill both local pressure directions. At wavenumber 2 the field components are spatially uniform in the invariant frame, so the calculation is the same at every point. An early wave and a wave timed to reach the shell during production were adjusted on a coarse birth sample to match the two-by-two target tensor. On an independent 4,096-pair sample, each sign's two principal pressures agree with their targets to within 1.3%. The momentum-distribution gaps are still 31.26 and 27.04 percentage points. Moreover, the two selected launches would cost 1,356 and 492 times the emitted energy in the original pulse. These large-amplitude controls use a prescribed original history. They illustrate the limited information in a pressure measurement. They do not give a self-consistent source for a hot shell.
* Faster angular variation can populate the weaker lateral direction even with a much smaller incident-energy budget. We paired an early wavenumber-2 field with a prepared wavenumber-515 pattern distributed across the sphere. At 128 common observation positions and 64 birth samples per position and sign, the mean local weak and strong pressures were **0.132 and 2.182** for dedicated interactions, and **0.192 and 2.249** for shared interactions. The transverse distribution gaps from those same samples were 14.20 and 15.36 percentage points. The respective launch costs are 0.2729 and 0.0950 times the original emitted energy, *in addition* to the original axial wave. At wavenumber 258 the corresponding weaker pressures were 0.031 and 0.050. Each wavenumber uses a different launch amplitude, delay and propagated response. These examples sample the joint effects of wavenumber, launch strength and timing. A new endpoint scramble gives weak pressures 0.137 and 0.198. At matching dedicated sample settings a finer radial field mesh changes the weak pressure from 0.139 to 0.161 and the momentum gap from 14.97 to 18.38 points.
* The evolving-mean source cases replay the electric, magnetic and particle-birth histories from the companion source calculation. Diverting 0.64% of the original emitted power into two early wavenumbers, 5 and 7, gives a 29.21-point distribution gap and an aligned weak pressure at most 0.000098 of its target. Refining the source radial mesh gives 29.25 points and 0.000099. The mean early-time geometry and birth history respond to the extra radiation in these cases. The orbits here provide no transverse current back to the fields.

The target is the parent model's moving, classical Maxwell–Jüttner population, evaluated at the marker for 1% of the original full-run absolute electrical work. A discrepancy of 31 percentage points means that for some cutoff in sideways momentum the cumulative fractions below that cutoff differ by 0.31. The comparison is a population diagnostic of the simulated particles. It has no interpretation as a p-value for cosmological data. [Numerical output](PAIRS.json) · [Timed propagation and source normalization](../fields/README.md) · [Evolving spherical source](../source/README.md)

## Dynamics and field geometry

Write the shell's angular position as two complex numbers with $|z_1|^2+|z_2|^2=1$. The orthonormal tangent frame is $e_1=(-\bar z_2,\bar z_1)$, $e_2=i e_1$, $e_3=-i(z_1,z_2)$. The original axial fields multiply $e_3$. The two added fields use horizontal, divergence-free curl eigenforms of eigenvalue $-\kappa$ in the orientation of this frame. The first shape has horizontal components $Y_1+iY_2=\sqrt{\kappa-1}\,\bar z_1^{\kappa-2}$. The second can use the same family, with an additional angular phase, or the spatially distributed member
$$
Y_1+iY_2=\overline{\sum_{j=0}^{n}
  e^{i\pi\sqrt{2}\,j^2}\sqrt{\binom{n}{j}}\,z_1^jz_2^{n-j}},
\qquad n=\kappa-2 .
$$

The two independent phase angles on the sphere make distinct monomials orthogonal. Their individual norms squared integrate to $j!(n-j)!/(n+1)!$ using normalized Haar volume. The distributed member consequently has unit mean squared magnitude. For every degree-$n$ homogeneous holomorphic polynomial $f$, the frame derivatives satisfy $e_2f=i e_1f$ and $e_3 f=-i n f$. Together with $[e_1,e_2]=2e_3$ and its cyclic versions, these identities give zero divergence and curl eigenvalue $-(n+2)$ for the horizontal vector with complex components $\bar f$. A conjugate homogeneous polynomial gains phase $+n\psi$ along the axial flow $z\mapsto e^{-i\psi}z$. A particle can therefore see a rotating lateral field while traveling a small distance on the shell. The distribution over monomials also lets many different portions of the surface feel the high-wavenumber wave.

The radial transfer is degenerate over the angular eigenspace on the prescribed round background, so each new angular shape can multiply the same computed radial waveform at its wavenumber. For unlike wavenumbers, or orthogonal quadratures of an equal wavenumber, the squared emitted amplitude of a pair is the sum of the squared launch coefficients. The large uniform control uses nonperpendicular vectors with **disjoint time supports** at the emitting tube, which also gives this sum. Live-source traces already include their launched amplitudes.

The prescription takes birth times from cumulative daughter production and uses unit Haar birth positions. The rest-frame birth kick has magnitude $\sqrt{R^2-1}$ in $mc$ units and uniform transverse azimuth, with the inherited polar acceptance proportional to $\min(1,|c-d\mu|/\mathrm{target})$. A scrambled Sobol sequence selects these variables. The two charges receive opposite kicks. Each time step drifts by a great-circle geodesic for half a step, applies a centered relativistic electric and magnetic Boris update to comoving momentum $\mathbf u=a\mathbf p/(mc)$, and drifts again. A partial first step starts at the actual birth time. The round metric has a bi-invariant frame, in which its free geodesics have constant momentum components. [Orbit integrator](sphere.cpp) · [Sampling, reference distributions and controls](pairs.py)

## Measuring local pressure

At one position, each species contributes a pressure tensor proportional to $\langle p_i p_j/\gamma\rangle$. Momenta here are dimensionless. In the two lateral directions the target pressure per lab-frame particle, in units of $mc^2$, is $\theta/\Gamma$ for comoving temperature $\theta$ and drift Lorentz factor $\Gamma$. The target lab-frame mean energy is $\Gamma(e_{\rm rest}+\theta)-\theta/\Gamma$, where $e_{\rm rest}=K_1(1/\theta)/K_2(1/\theta)+3\theta$. The two charge species use opposite reference drifts. The normalized full three-dimensional comparison transforms the momentum numerators to that species' drift frame while retaining the invariant $d^3p/\gamma$ tensor measure.

There are three spatial comparisons in the output.

1. The pooled tensor averages over births and positions in the invariant frame and then takes eigenvalues. In a varying field its eigenvalues give little direct information on the sidewise pressure at one point.
2. The aligned tensor rotates each particle's lateral axes according to its sign and observed position. The local principal direction is estimated from the impulse covariance along axial ultrarelativistic test characteristics with the recorded birth distribution. The heat diagnostic itself uses full tracked particle momenta. The smallest eigenvalue of a sum is at least the sum of the smallest eigenvalues. Therefore the aligned pooled weak eigenvalue *upper-bounds* the mean of the true local weakest pressure, for the underlying population and in the same number weighting. The empirical estimate of that bound has finite sampling error. A small value demonstrates local failure even when these approximate axes differ from the exact eigenvectors.
3. The conditioned calculation specifies common *final* positions, inverts the birth-position map for every charge and birth label, and computes the pressure at each position before taking eigenvalues. A finite-difference Jacobian in orthonormal geodesic coordinates gives the factor $1/|\det J|$ that converts initially uniform volume to final local number density. We average the local eigenvalues by number density, and separately record the 10th, 50th and 90th percentiles over uniformly sampled positions and both signs. The maximum endpoint residual, determinant range, empirical maximum Frobenius norm of the differential minus the identity, and total density normalization are provided. This inverse construction assumes the sampled branch covers the particles reaching that position. Near-identity differentials and agreement with forward sampling test the assumption along the sampled trajectories. They do not prove global invertibility.

The late uniform pressure control fits the three independent components of the mean of the two species' tensors by nonlinear least squares over the early amplitude and the late vector's two Cartesian components. The equal-early control varies the common amplitude to set the trace to two targets. Default runs replay the frozen coefficients on independently scrambled samples. The high-wavenumber survey chose amplitudes of the order needed to produce thermal-scale lateral momenta and used the same fixed pattern of angular phases. The densities, source timing, particle masses and original dominant drift remain benchmark inputs. In particular, collisionless anisotropic distributions with two selected moments close to thermal remain compatible with a substantial discrepancy in the full distribution.

## Reproduction

Use the parent study's compact dedicated and shared histories. After the primary and storage-scan fields and live-source traces have been generated, run

```sh
OPENBLAS_NUM_THREADS=1 python orbits/pairs.py --history history \
    --fields output/primary/traces --scan output/scan/traces \
    --source output/source/source_traces --output output/PAIRS.json --jobs 2
```

The optional `--refit` repeats the two pilot solves. Standard runs include a second scramble for the locally uniform control, endpoint and birth-sampling changes for wavenumber 515, a smaller orbit time step with identical initial conditions, a coarser source time grid, a finer radial field mesh, and the two source-model interaction histories. The outputs record the seed, sample count and step count for each row. The C++ standard library, a C++17 compiler, Python, NumPy and SciPy suffice.

## How much fine-grained phase space can all the births occupy?

The angular and pressure comparisons use a few observables. A separate phase-space calculation tests the complete joint distribution of a particle's position and momentum. Assume the injected daughters have a definite momentum magnitude at each birth time, as in the inherited decay law. Take one charge species on the prescribed round wall and let the fields depend on position and time. The following bound includes arbitrarily complicated subsequent collisionless magnetic deflection.

Use normalized Haar volume $dx$ on the unit sphere and comoving momentum $\mathbf u=a\mathbf p/(mc)$. The Lorentz flow preserves $dx\,d^3u$. The angular flow in the invariant frame has zero spatial divergence, the electric acceleration is independent of momentum, and the relativistic magnetic rotation has zero momentum divergence. Let $\Phi_{T,s}$ be the phase flow from birth time $s$ to the marker $T$. Allow the full sphere of initial directions $\mathbf n$ and all birth positions. All the particles of the selected sign at the marker lie in the image
$$
S_T=\{\Phi_{T,s}(x,U(s)\mathbf n)\},\qquad
U(s)=a(s)\sqrt{\max(0,R(s)^2-1)}.
$$

At each $s$, derivatives in the three initial spatial coordinates and the two direction angles span the tangent space to this five-dimensional injection shell. When the birth time varies, pulling the derivative back through $\Phi_{T,s}$ gives the velocity of the injection shell minus the phase-flow velocity at birth. The component normal to the momentum sphere is $U'-q a\mathbf E\cdot\mathbf n$, with $|q|=Q$. The magnetic term is tangent to this sphere. The absolute swept-volume Jacobian therefore equals $U^2 |U'-q a\mathbf E\cdot\mathbf n|$ relative to $dx\,d\Omega\,ds$. Volume preservation and the area formula imply
$$
\mathrm{Vol}(S_T)\ \le\ 4\pi\left[
  \tfrac13\int |d(U^3)|+
  \tfrac{Q}{2}\int U^2 a\,\|\mathbf E\|_{\rm RMS}\,ds
\right].
$$
Overlapping images from different birth times only increase the volume counted on the right. The factor $1/2$ uses the spherical mean of $|\cos\vartheta|$. The RMS refers to Haar volume at the time of injection. The factor bounds the mean of $|\mathbf E|$ even with angular concentration. The empirical birth probability and angular acceptance can select a subset of this support.

The target has uniform position and the moving Maxwell–Jüttner momentum density, whose maximum per $d^3[p/(mc)]$ equals $1/(4\pi\Gamma\theta e^{1/\theta}K_2(1/\theta))$. Divide the square bracket above by $a(T)^3\Gamma\theta e^{1/\theta}K_2(1/\theta)$ to obtain an upper bound $b$ on the target's probability of lying in the particles' occupied set. The source assigns probability one to that set. Hence the total-variation distance of the full same-sign position–momentum distributions is at least $1-b$.

| Source or illustrative example | Dedicated upper bound $b$ | Shared upper bound $b$ |
| --- | ---: | ---: |
| Original axial model | $2.282\times10^{-7}$ | $5.099\times10^{-7}$ |
| Two time profiles fitted to both lateral pressures | $2.332\times10^{-7}$ | $5.144\times10^{-7}$ |
| Any unit-energy combination in the primary bank | $9.291\times10^{-6}$ | $2.078\times10^{-5}$ |
| Any unit-energy combination in the storage bank | $2.121\times10^{-5}$ | $4.729\times10^{-5}$ |

The table uses the half-step histories. The two bank bounds allow arbitrary members of each angular eigenspace with the tabulated electric transfer, with a total additional emitted energy up to one original pulse. The round background gives the same electric transfer for the two helicities in the inherited Maxwell equations. In an orthonormal angular basis, the source energy for the time coefficients $c$ of any one component is $c^\mathsf{T}Wc$. If the unit launches give the vector $e_\kappa(s)$ of electric transfer coefficients, the additional field RMS is bounded by $\epsilon\max_\kappa\sqrt{e_\kappa^\mathsf{T}W_\kappa^{-1}e_\kappa}$ at total emitted energy $\epsilon^2$. This follows by Cauchy–Schwarz in each component and summing the squared norms of orthogonal angular components. The bound adds the axial norm by the triangle inequality. The original force and the two trial transverse forces are combined with their actual overlaps for the individual cases.

The calculation encloses the piecewise-linear input histories. Within each cell it takes maxima of the nonnegative factors in the electric integral. The total variation of $a^3(R^2-1)_+^{3/2}$ is bounded using the endpoint maximum of each factor times the other factor's absolute change. The electric response norm, including its maximum over wavenumbers, is convex along a linear segment and is bounded by the interpolated endpoint envelope. The code also reproduces the exact swept volume of a force-free, growing momentum sphere.

This is a restriction on **fine-grained**, collisionless phase space in a specified injection model. A folded distribution can agree with a smooth thermal distribution after binning, projection onto momentum alone, or observational averaging. A spread of parent momenta, spatially varying birth energies, stochastic scattering or a collision operator changes the assumed source or flow. This calculation supplies a precise reason to test those ingredients.

The same prepared inputs reproduce the enclosure in a separate step.

```sh
OPENBLAS_NUM_THREADS=1 python orbits/liouville.py --history history \
    --fields output/primary/traces --scan output/scan/traces \
    --source output/source/source_traces --pairs output/PAIRS.json \
    --output output/LIOUVILLE.json
```

[Phase-space calculation](liouville.py) · [Numerical bounds](LIOUVILLE.json)

## What survives coarse momentum measurements?

The support-volume theorem concerns joint phase space. Counting in just 20 broad momentum cells gives a different, directly checkable comparison. We use five ranges of reference-rest-frame kinetic energy $x=(\gamma'-1)/\theta$ with boundaries at 0, 1, 2, 3, 5 and infinity. Within every range we divide the polar cosine $\mu$ relative to the sign-specific stream direction into four equal intervals. All comparisons use the same cells. The two signs contribute equal numbers.

The transformed coordinates describe each particle in its reference stream frame, but the sampled daughters are still counted on the model wall's common time slice. Relativistic transformation of the momentum-volume element gives the thermal count density
$$
dP=\frac{e^{-x}\,\gamma' p'}{2e^{1/\theta}K_2(1/\theta)}
(1+\beta v'\mu)\,dx\,d\mu,\qquad
\gamma'=1+\theta x,\quad p'=\sqrt{\gamma'^2-1},\quad v'=p'/\gamma'.
$$
Here $\beta$ is the absolute speed of the reference stream in the wall frame, in units of light speed. The factor containing $\mu$ accounts for the lab counting measure. Integrating this formula gives the thermal count in each of the 20 cells. The reported discrepancy is half the sum of the absolute differences of observed and target fractions. Equivalently, it is the largest count difference across any union of these cells.

| Example on fixed mean histories | Dedicated discrepancy | Shared discrepancy |
| --- | ---: | ---: |
| Uniform mode pair fitted to the two local lateral pressures | 46.77 percentage points | 55.71 percentage points |
| Same early pulse fitted to the pooled lateral pressure | 34.66 points | 34.49 points |
| Wavenumbers 2 and 40 selected for their relatively close transverse-momentum CDF | 49.38 points | 47.38 points |
| Distributed wavenumber 515 with the early wavenumber-2 field | 25.89 points | 26.75 points |

For example, in the wavenumber-2 and 40 test, the single cell $3\le x<5$ and $0.5\le\mu\le1$ contains 18.31 percentage points more dedicated daughters than the reference predicts, and 19.57 points more shared daughters. These are ordinary broad count differences visible after averaging over all positions and transverse azimuths. The 8,192-pair run and an independently scrambled 4,096-pair run with a larger orbit time step agree within 0.22 percentage points on the 20-cell discrepancies in the table. The numerical check applies to these selected fields and bin sizes. The source-bank theorem above has a wider field scope and a finer resolution.

```sh
OPENBLAS_NUM_THREADS=1 python orbits/coarse.py --history history \
    --fields output/primary/traces --scan output/scan/traces \
    --pairs output/PAIRS.json --output output/COARSE.json --jobs 2
```

[Coarse count comparison](coarse.py) · [All measured and thermal cell counts](COARSE.json)
