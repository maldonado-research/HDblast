# Tidal curvature on the constrained slice

A zero of the expansion rate at a time-symmetric bounce makes normalization by that rate singular. We instead use the fixed positive curvature scale $h^2=u_*$ of the homogeneous solution. The electric Weyl tensor measures the directional part of gravitational tidal acceleration. In spherical symmetry its magnetic counterpart vanishes. The difference $D=(E^r{}_r-E^\theta{}_\theta)/h^2$ determines the scalar $C_{abcd}C^{abcd}/h^4=16D^2/3$.

With $y=\cos\chi$ and the conformal factor $\psi(y)$ from the initial constraint, the trace-free spatial Gauss equation at $K_{ij}=0$ gives
$$
D=\psi^{-4}(1-y^2)\left[-2\frac{\psi_{yy}}{\psi}+6\left(\frac{\psi_y}{\psi}\right)^2-32 A^2 y^2\right].
$$
The last term is half the trace-free scalar gradient contribution with its sign fixed by $R_{ij}=\partial_iq\,\partial_jq+3u\,g_{ij}$. A second evaluation uses only first derivatives and the angular spacetime sectional curvature $(1-R_s^2)/R^2$, where $R=h^{-1}\psi^2\sin\chi$ is the areal radius and $s$ proper radial distance. Subtracting the four-dimensional Ricci traces gives
$$
-\frac23 D =
\psi^{-4}\left[1+4yp-4(1-y^2)p^2+\frac{64}{6} A^2 y^2(1-y^2)\right]-\frac{u}{u_*},
\qquad p=\psi_y/\psi .
$$
These expressions extend smoothly to the poles. Their equality also tests the constraint solution independently of the spectral equation residual.

Writing $\beta=\mu^2/h^2$, the conformal solution through second order yields
$$
D=A^2(1-y^2)\left[\frac{4(\beta-24)}{35}+\frac{16(\beta-10)}{7}y^2\right]+O(A^3).
$$
The code compares the two finite-amplitude formulas and checks the sign-averaged profile against the analytic coefficient. It computes the RMS using the proper volume $\psi^6 d\mu_{S^3}$, where $d\mu_{S^3}$ is round volume. The uniform control has zero Weyl curvature. On both branches the ripple gives nonzero tidal curvature already on the initial slice. This invariant quantifies the anisotropy of curvature in the spherical models. Determining gravitational-wave power would require independent radiative degrees of freedom.

[Solver](tidal.py) · [Both branches and signed amplitudes](TIDAL.json) · [Conformal initial data](../initial/INITIAL.json)
