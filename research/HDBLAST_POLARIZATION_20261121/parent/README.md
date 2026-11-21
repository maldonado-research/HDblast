# Can an early higher-dimensional pulse leave a field for particles born later?

Yes, in this HDblast model a patterned pulse can leave an electric field on the expanding shell before charged particles appear. The source must emit some of its wave in a second direction of polarization. The simulation starts that component at zero everywhere, launches it at the outer edge of a five-dimensional region, follows its passage to the shell, and measures its strength when particles form. The unforced component stays zero. This tests a specific step in the proposed origin of a hot universe. [Question recorded before evaluating the new runs](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7066509412) · [Previous benchmark](https://zenodo.org/records/24415138)

## Where does the energy come from?

Of the 39 angular scales tested consecutively on the round shell, mode 7 left the strongest electric field when particle formation began. Launching this pattern at 1% of the original pulse amplitude gives a field of about 0.00997 in the model's units. A rotation of the local source polarization by at most 1.40 degrees could supply this component at fixed total launch power, to first order in the waveform used in the simulation. How a parent event would produce this angular structure is open.

The field gets its energy mainly from a changing property of the wall, the coefficient that determines how much energy an electromagnetic field stores. For a unit launch in mode 7, this change supplied about 14,499 units of comoving energy. The net flow through the radial boundary removed 269, leaving 14,230 stored when particles appeared. Holding the coefficient at its value just before its final decrease made the electric field at that moment about **1,300 times weaker** on the same geometry. The scalar and geometry that set the coefficient follow recorded histories. A complete calculation would also subtract the energy used by the additional field from its source. [Five-dimensional propagation, angular spectrum and controls](geometry/README.md) · [When the pulse acted](transfer/README.md) · [Shaping a wave at the wall](control/README.md)

## Does that field make a thermal gas?

The earlier benchmark describes positive and negative particles moving in opposite directions. A hot gas in motion can already have strong pressure along its direction of travel. We compare each measured quantity with the thermal prediction in the same frame.

| At 1% of whole-run absolute electrical work, two sampled single-polarization waves in a local patch | Separate interactions | Shared interaction |
|---|---:|---:|
| Mean particle energy above the moving-gas prediction | 2.21% | 0.67% |
| Difference in the largest principal pressure fraction | 0.075 percentage points lower | 0.0053 percentage points lower |
| Weakest principal pressure fraction divided by its thermal prediction | 1/56,400 | 1/35,300 |
| Largest difference in cumulative counts of sideways momenta | 13.1 percentage points | 12.9 percentage points |

The two cases use different input waves chosen in the sampled family. The field can give individual particles large sideways kicks. It puts nearly all the sideways pressure in one direction. A thermal particle stream should exert equal pressure in its two sideways directions. The test therefore distinguishes two models that would look nearly identical if one measured only their total energy. The local calculation allows the currents of the new particles to drive fields in a finite set of spatial modes. A separate test tracks particles around the entire spherical wall through a smooth field whose direction changes with position. In the two global mode-5 examples, the pressure across the local pattern is below 0.04% of its moving-gas value. Both tests use collisionless trajectories, fixed mean birth and expansion histories, and the stated electromagnetic forcings. These pressure diagnostics refer to the modeled particle stream. [Current-coupled local results and frame comparison](orbits/README.md) · [Global charged-particle trajectories](sphere/README.md)

![Electric field remaining versus angular scale and comparison of energy and principal pressures with a moving thermal gas](prebirth-response.png)

## How strong must an additional pattern be?

We also tested all superpositions within a finite band of angular patterns, using the same timing for every component at the outer source. This covers **44,200 independent real angular coefficients** with either orientation of curl. Particles begin at uniformly distributed positions and can move at any allowed speed after birth. The bound uses the short distance they can travel by the comparison time.

For a source change with 1% of the original RMS amplitude, any such pattern in the prescribed collisionless model leaves a discrepancy of at least **8.09 or 2.28 percentage points** in the counts below some sideways momentum, on the finer time grid. To bring every count within one percentage point of the moving thermal target, the added source would need **at least about 2.7% or 1.5%** of the original amplitude. These are necessary conditions. The coarser histories give 2.66% and 1.47%. The result concerns the chosen wave timing and angular scales, and uses fields propagated through the five-dimensional region with a prescribed mean wall evolution. The arithmetic was cross-checked by direct cohort integration and separate numerical integration of the thermal angular distribution. [Derivation, grid checks and count witnesses](sphere/README.md)

The concrete next test is to send two independently polarized components, evolve the particles and their currents, and measure both lateral pressures in the same stream while charging the field energy to its source. An astronomical prediction requires the model to specify the original angular emission and connect the resulting matter to a later observable. The calculations here supply simulation benchmarks for those steps. [Primary literature and comparisons](SOURCES.md)

## Reproduce the results

Use Python 3.12, NumPy 2.3.5, SciPy 1.17.0 and a C++17 compiler available as g++. From this directory, verify the study manifest and rebuild the results in a new or empty output directory. The replay reconstructs four reference wall histories, checks their compact records, propagates the additional fields and recomputes the particle and bound comparisons. Allow several gigabytes of temporary disk space. The source for the figure uses Matplotlib. [Figure code](figure.py)

```sh
python -I -B -S study.py
python -I -B study.py --verify --output-dir ../prebirth-replay --jobs 2
```
