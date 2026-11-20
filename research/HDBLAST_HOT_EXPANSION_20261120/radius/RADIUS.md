# Can lower initial spatial curvature keep the wall hot?

A closed wall starts at rest and receives one finite pulse. Increasing its starting radius lowers its spatial curvature, but it also changes the mass available to supply the pulse. In the published preparation, increasing the initial radius from 2 to 2.01 reduces the peak fraction of wall energy in internal particle motion from **91.49% to 15.6%**. Near 2.02455, the pulse exhausts the positive enclosed mass of the simulated bulk during launch. Material at the wall has scarcely begun to move at this point.

## Initial conditions and source normalization

All radii and times use the preceding calculation's bulk curvature length. The initially cold neutral parents have rest energy density $D_i=10^{-4}$. Write $a_i$ for the initial radius, with $H_i=0$. The closed-wall constraint, initial mass and density available to charged daughters are

$$
H^2+a^{-2}=e(2+e)+\mu/a^4,\qquad
\mu_i=a_i^2-a_i^4 D_i(2+D_i),\qquad
D_*(a)=\frac{D_i}{0.95}\left(\frac{a_i}{a}\right)^3 .
$$

Initial positive mass alone gives $a_i<[D_i(2+D_i)]^{-1/2}\simeq70.71$. The finite pulse puts a much tighter limit on the studied branch. The wall and every radial element of its Gaussian-normal collar start from these data. The conservation checks use $a_i^3D_i$ as initial comoving rest energy. All runs use the same compact source waveform and duration in source proper time. The local dose is $\mathcal D=2\int I^2 d\tau_s$. Per unit solid angle on the unit three-sphere, the emitted proper energy is $\mathcal E_s=2\int r_s^3 I^2 d\tau_s$, where $r_s$ is the source radius. Equal local dose and equal total emitted energy therefore describe different physical preparations.

The code measures mass on every collocation point and at every step during launch, and every 32 steps after time 0.025. A detected negative mass terminates the run. A rebound requires prior contraction, a drop below $(1-10^{-5})a_i$, and a negative-to-positive zero of $H$ after $H<-10^{-4}/a_i$. All reported rebounds comfortably exceed these cuts. Subsequent maxima require a recorded rebound. Runs with no relevant event stop at the requested time, contraction below $0.8a_i$, or expansion past $100a_i$. The widening of the expansion cutoff allows prolonged scalar trapping to be distinguished from exit. An additional time step limit resolves field damping when a weak pulse lets the scalar roll back past its entrance.

## Fixed local dose

Here the pulse has $\mathcal D=0.123$ and the entire material model is inherited from the preceding study. The hot fraction is the intrinsic random energy density of charged particles divided by total added wall energy density. It is measured on the expanding branch after scalar exit. The largest continuous radius ratio with fraction above 90% is reported. A dash indicates that the corresponding event was absent during the run.

| Initial radius | Minimum bulk mass | Rebound radius | Exit radius | Peak hot share | Growth above 90% |
|---:|---:|---:|---:|---:|---:|
| 2 | 0.09654 | 1.881 | 4.921 | 91.49% | 1.00777× |
| 2.0025 | 0.08698 | 1.884 | 5.079 | 76.70% | — |
| 2.005 | 0.07736 | 1.887 | 5.359 | 49.00% | — |
| 2.01 | 0.05793 | 1.893 | 6.370 | 15.58% | — |
| 2.015 | 0.03826 | 1.898 | 9.288 | 2.26% | — |
| 2.02 | 0.01835 | 1.904 | — | — | — |
| 2.024 | 0.00223 | 1.908 | — | — | — |
| 2.025 | crossed zero | — | — | — | — |
| 3 | crossed zero | — | — | — | — |
| 4 | crossed zero | — | — | — | — |


The source-mass balance also explains the scaling. If $u_s$ is the proper radial speed of the source, $r_{s,x}$ its radial derivative in the inward normal coordinate, and $I_s,O_s,B_s$ the two traveling amplitudes and magnetic field, its exact equation in this model is

$$
-\frac{d\mu_s}{d\tau_s}
=4r_s^3\!\left[(u_s-r_{s,x})I_s^2+(u_s+r_{s,x})O_s^2-\tfrac12 u_s B_s^2\right].
$$

During a short launch toward the wall with little reflected flux, taking $u_s\simeq 0$ and $r_{s,x}\simeq-r_s$ yields the scale estimate $\Delta\mu_s\simeq 2r_s\mathcal E_s\simeq2a_i^4\mathcal D$. Holding the local fluence fixed therefore makes the mass cost grow roughly as the fourth power of radius while the initial mass grows roughly as the second power over this range. The full integrations retain the motion, reflection and magnetic terms.

The dimensionless ratio $(\mu_i-\min\mu)/a_i^4$ lies between 0.24376 and 0.24378 at the seven positive-mass initial radii from 2 to 2.024 in the table. Meanwhile the available $\mu_i/a_i^4=a_i^{-2}-D_i(2+D_i)$ decreases with radius. Treating the mass cost measured at radius 2 as constant predicts a crossing at 2.02458. Direct integration over $0\le\tau\le0.014$ gives the launch threshold near **2.02455**. Three mesh bisections with 16, 24 and 32 radial cells are tabulated in the JSON. Independent runs at 96 cells and a smaller time step preserve positive mass at 2.02454 and cross zero at 2.02456. This is a numerical bracket for the local crossing, with no assertion of global monotonicity in a different family of sources. The negative-mass runs could retain a positive mass at the wall while the inner source already failed. For example, with initial radius 3 the measured wall mass was 8.98 when the source mass crossed zero. A wall-only check would miss this exclusion.

## Fixed total emitted energy

The source dose at each radius was adjusted until the total incoming proper energy per unit solid angle matched the original 0.978037. The tabulated doses were obtained by starting at $0.123(2/a_i)^3$ and twice multiplying by the ratio of target to measured emitted energy over $0\leq\tau\leq0.015$, with 32 cells and CFL 0.08. The frozen table supplies the resulting amplitudes for direct replay. The unmodified scalar action, late gauge coupling, neutral precursor, collision and decay parameters are used for these controls.

| Initial radius | Local dose | Emitted energy | Minimum bulk mass | Largest scalar value | Rebound radius |
|---:|---:|---:|---:|---:|---:|
| 2.02 | 0.119382316 | 0.978037 | 0.1377 | 0.24136 | 1.917 |
| 2.1 | 0.106250903 | 0.978037 | 0.3104 | 0.00987 | — |
| 3 | 0.036442453 | 0.978037 | 3.1297 | 0.00685 | — |
| 4 | 0.015373954 | 0.978037 | 8.1416 | 0.00584 | — |

With the same finite energy spread over a larger surface, the wall at radii 2.1, 3 and 4 cannot even climb the scalar entrance to its plateau at $\psi=0.012$ before the integration reaches the contraction limit. Starting at 2.02 yields a rebound and an extended plateau, without completing the exit before the calculation stops.

## Explicitly retuned preparation

The retained comparison at radius 4 changes the local dose to 0.03075 and the entrance calibration to 0.03. An initial inverse-square reduction of both the local dose and the plateau height led to a rebound without a completed exit. For the revised comparison, the plateau height is 0.0266955, the late gauge increment is 111985.95 and the scalar decay threshold is 0.85054. The scalar threshold, plateau and gauge increment were selected after exploring trajectories. These parameters would require an independent physical origin. The pulse width, parent density, charged mass ratios, decay coefficient 100 and collision rate 6 retain their original values. This retuned pulse emits 1.956188 per unit solid angle, about twice the original total. The outcomes combine changes in geometry, energy supply and material parameters.

| Charge-to-mass parameter | Scalar gate | Cells before / after remesh | Peak hot share | Growth above 90% |
|---:|---:|---:|---:|---:|
| 0 | 0.85054 | 32 / 6 | 0.00% | — |
| 2.5e+07 | 0.85054 | 32 / 6 | 79.51% | — |
| 1e+08 | 0.85054 | 32 / 6 | 91.54% | 1.02846× |
| 6e+08 | 0.85054 | 24 / 6 | 91.73% | 1.03804× |
| 6e+08 | 0.85054 | 32 / 6 | 91.73% | 1.03797× |
| 6e+08 | 0.85054 | 48 / 6 | 91.73% | 1.03793× |
| 6e+08 | 0.85054 | 64 / 6 | 91.73% | 1.03791× |
| 6e+08 | 0.85054 | 48 / 10 | 91.73% | 1.03793× |

The 64 / 6 run starts at rest, contracts to 3.9239, and rebounds. Its minimum sampled bulk mass is 0.33345. When the scalar has completed its exit at radius 12.52, the expansion speed in curvature units is $aH=2.59$. When the matter first reaches 90% internal energy, the radius is 28.43 and $aH=0.290$. The first turnaround follows at radius 29.51. The delay between scalar exit and heating has consumed most of the remaining expansion.

With zero electric current, the peak internal energy fraction is $1.68\times10^{-5}$, supplied by daughter injection and its subsequent evolution. A family of conducting cases transfers energy into random particle motion under the assumed collision closure. The intrinsic fraction is a model energy ratio. A temperature measurement, equilibration rate derived from microscopic scattering, and a realistic long radiation era are distinct tests. The hot-expansion factor in this family falls short of two. The tests constrain these specified source and material choices.

## Reproduction and numerical checks

From the capsule root, run

```sh
python -I -B radius/radius.py --output /tmp/RADIUS.json
```

The runner checks the SHA-256 of the bundled preceding C++ solver, applies the radius patch in a temporary directory and compiles with a C++17 compiler. Python uses the standard library. The system tools are GNU-compatible `patch` and `g++`. Each retained parameter set, event, energy ledger and mass bracket is in the JSON. Rows ending at negative mass report the first excluded state sampled by the diagnostic, so its magnitude depends on the time step.

With the original radius and 64 early cells, absolute differences from the archived primary in minimum mass, bounce, scalar-exit and maximum radii, peak fractional internal energy and duration in e-folds are all below 1e-13. At 64 / 6, the large-radius comparison has peak absolute Friedmann residual 1.7e-11, surface-energy residual 1.1e-09, sheet-energy residual 6.7e-09, gas-energy residual 5.8e-10 and independently evolved tube-mass residual 6.4e-06. The tube-mass residual can be compared with minimum bulk mass 0.333.

The launch scan varies radial resolution and time step. For late-time transport the primary large-radius comparison uses six radial cells after time 0.08. At fixed 48 early cells, increasing the late resolution from 6 to 10 changes the hot radius ratio by 9.5e-12 and the peak intrinsic fraction by 5.2e-13. Extra digits in the archive support comparisons across runs. The full collar is finite and driven by a prescribed source boundary. Its local positive-mass checks do not construct a regular bulk interior beyond that boundary.
