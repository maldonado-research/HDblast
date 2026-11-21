# Could a different time profile store more field at birth?

This calculation places the new input directly at the wall. The comparison fixes the time-dependent upper limit on its amplitude to the measured inward axial characteristic. On each recorded time cell the additional component can reverse sign or have any smaller amplitude. We ask for the largest combined electric and magnetic amplitude present when particle production begins. Such a chosen waveform must still be produced by the bulk and source to constitute a physical solution.

The response is linear on the recorded wall geometry. Propagate the two terminal field components backwards through the transpose response, or equivalently multiply each cell's forward source column by the transfer map from that cell to birth. For a cell, midpoint sampling of the source and the exact exponential of the homogeneous midpoint matrix give its terminal vector. All sums of these vectors with individual factors in the interval from −1 to 1 form a centrally symmetric polygon. Its farthest vertex gives the maximum terminal norm under this cellwise envelope. Sorting the segments by direction visits every vertex. A brute-force sign enumeration on a separate small problem and projection of the maximizing vertex onto its own support direction check the optimization. The lowest physical mode driven with the recorded signs reproduces the reference terminal electric and magnetic vector to relative accuracy better than one part in a million.

| Comoving wavenumber | Unchanged waveform, stored spatial RMS | Best sign-shaped waveform, stored spatial RMS | Fraction of maximum |
|---:|---:|---:|---:|
| 2 | 0.14643 | 0.14849 | 98.6% |
| 5 | 0.61407 | 0.62449 | 98.3% |
| 10 | 0.57881 | 0.59988 | 96.5% |
| 20 | 0.000712 | 0.01196 | 6.0% |
| 50 | 0.002316 | 0.002316 | 99.99% |

With the favorable mode 5, shaping the signal's sign leaves little improvement under this particular envelope. At mode 20 much of the response cancels unless the driving sign changes. The table uses the time-refined recorded geometry and a standing wave whose spatial peak is normalized to the reference characteristic. A cosine has spatial RMS equal to its peak divided by the square root of two. The source here acts at the inner boundary. The separate true-collar calculation applies its source at the outer tube, propagates it, and obtains a different waveform at the wall. [Reproduction code](optimal.py) · [All results](OPTIMAL.json) · [Outer-tube propagation](../geometry/README.md)
