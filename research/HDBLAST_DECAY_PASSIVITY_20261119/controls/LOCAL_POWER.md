# Local daughter power on the recorded B1 walls

Consider a real scalar initially in the ground state of its instantaneous physical Hamiltonian. A bath that damps it toward a squeezed Gaussian supplies energy to this scalar. If the daughters alone pay this cost, their local power is negative. Assigning the daughters positive energy from the scalar's occupation in the bath's particle basis gives the opposite sign, and an additional reservoir would have to supply the difference. This test quantifies the difference for one proposed bath on the HDblast shell.

## Readouts

The table gives the positive occupation assignment per physical volume and proper time in units of $H_0^5$. The last column gives the gap between that assignment and the negative daughter power required to balance the physical parent. The walls come from the fine reciprocal B1 calculation with its original decaying scalar. Each row evaluates a hypothetical physical-Hamiltonian vacuum at that point. The state at each readout is specified independently. No trajectory of that test state or daughter density is inferred here.

| Fixed wall and readout | $Q_+$ at $K=8$ | $Q_+$ at $K=12$ | $Q_+$ at $K=14$ | $2Q_+$ at $K=8$ |
|---|---:|---:|---:|---:|
| Common preparation, $s=0.769028$ | 0.888678 | 2.711226 | 4.055728 | 1.777356 |
| Quiet at $s=2.10$ | 0.00853304 | 0.02806640 | 0.04384125 | 0.01706609 |
| Near positive at $s=2.10$ | 0.00855392 | 0.02813473 | 0.04394764 | 0.01710785 |
| Near negative at $s=2.10$ | 0.00843623 | 0.02774955 | 0.04334789 | 0.01687247 |
| Far positive at $s=2.10$ | 0.00852961 | 0.02805515 | 0.04382373 | 0.01705922 |

The coarse calculation gives common-preparation $Q_+(8)=0.888277$ and quiet $Q_+(8)=0.00852687$. On both grids and at all three cutoffs the local gap is largest for the near positive pulse, followed by quiet, far positive and near negative at $s=2.10$. On the fine grid with $K=8$, the pulse changes in $Q_+$ relative to quiet are $+0.2447\%$, $-1.1345\%$ and $-0.04026\%$ for near positive, near negative and far positive, respectively. Those percentages describe a property of the local test state and geometry. The earlier published radiation density used an evolved state and an instantaneous-physical-vacuum bath. The sign test here concerns a proposed change of bath. [Published source study](https://doi.org/10.5281/zenodo.24355836)

## Oscillator and energy convention

Write $s=H_0\tau$, $h=d\ln a/ds$, $a=\exp(A_b-A_{*,\mathrm{quiet}})$ and $\mu=100(\phi-0.5)$. The signed mass yields $M^2=\mu^2$, $\dot M^2\equiv d(M^2)/ds=2\mu\dot\mu$, with $\dot\mu=100\,d\phi/ds$. For dimensionless comoving $k$, the physical momentum is $p=k/a$, $E^2=p^2+M^2$, and $\Omega=aE$. All frequencies and physical momenta in these equations are in units of $H_0$. The archive supplies $h$ and $d\phi/ds$ directly. Both logarithmic scale factors at the reference crossing come from the corresponding original quiet classical grid. Linear interpolation uses the two rows bracketing $s=2.10$ in proper time.

The physical Hamiltonian per mode is $H_b=\Pi^2/(2a^3)+a^3E^2\chi^2/2$. Choose $u=a\chi$ and the local zeroth-order WKB data $u=(2\Omega)^{-1/2}$, $u'=(-i\Omega-\Omega'/(2\Omega))u$, where a prime is a conformal-time derivative. These are instantaneous Gaussian data for a minimally coupled scalar, using $W=\Omega$ at lowest adiabatic order. The model assigns damping locally to this specific target. Set

$$
 C=\frac{2hM^2+d(M^2)/ds}{4},
 \qquad \ell=h+\frac{1}{2}\frac{d\ln\Omega}{ds}
            =h+\frac{C}{E^2},
 \qquad n_B=\frac{\ell^2}{4E^2}.
$$

The target has the same field variance as the physical vacuum and an extra real momentum slope $-a^3\ell\chi$. The overlap occupation $n_B$ is symmetric between these two pure Gaussians. With an attenuator $\dot V|_{\rm bath}=-\Gamma(V-V_B)$ and $\Gamma=y^2M^2/(8\pi E)$, the dissipator changes the physical parent energy by $+\Gamma E n_B$ in the initial physical vacuum. Assigning $-d\langle H_b\rangle_{\rm bath}/ds$ to the daughters therefore gives $-\Gamma E n_B$. The positive occupation assignment gives $+\Gamma E n_B$ per mode. These expressions refer only to bath power. Work from the explicit time dependence of the Hamiltonian belongs separately in the full energy balance.

Integrating the positive assignment over physical phase space of one real scalar, with $y=1$ for the test,

$$
 Q_+(K)=\frac{y^2M^2}{64\pi^3}\int_0^{K/a}
      \frac{p^2}{p^2+M^2}\left(h+\frac{C}{p^2+M^2}\right)^2dp .
$$

For $P=K/a$, $M=|\mu|>0$, $L=P^2+M^2$, and $t=\arctan(P/M)$, write $Q_+=y^2 M^2(h^2 I_1+2hCI_2+C^2 I_3)/(64\pi^3)$. The three primitives of $p^2/(p^2+M^2)^j$ from zero to $P$ are

$$
 I_1=P-Mt,\qquad I_2=\tfrac12(t/M-P/L),\qquad
 I_3=\tfrac18(t/M^3+P/(M^2 L))-P/(4L^2).
$$

All proposed test baths use $y=1$. For another $y$ the local power scales as $y^2$ on the same wall. The data include historical $y=0,1$ controls with matter feedback switched off and on. These labels identify which previously recorded geometry supplies $a,h,\phi,\dot\phi$. The K=12 and 14 calculations extend the integral on the original K=8 histories. They keep the geometry fixed. At preparation $M/H_0$ in dimensional language is approximately 46.09, and at the common readout it is 50.97–50.98. Even at $K=14$, the physical $P$ on the primary fine wall is only 26.33 at preparation and about 11.33 at readout. An ultraviolet asymptotic estimate would miss the regime of the table.

## Numerical checks and reproduction

Composite Simpson quadrature at 128 and 256 intervals agrees with the elementary primitive to maximum relative differences $1.30\times10^{-10}$ and $8.08\times10^{-12}$ across the 120 slice/band entries. A second calculation constructs $(\chi,\Pi)=(u/a, au'-a^2 h u)$, checks the Wronskian, computes the physical energy difference from its two canonical variances and obtains the occupation from the negative frequency coefficient. Integrated in comoving $k$ with the $a^{-3}$ volume conversion, these two signed rates agree with the analytic values to $8.19\times10^{-12}$ relative. Replacing linear time interpolation by four-point local polynomial interpolation changes the cut power by at most $9.60\times10^{-6}$ relative. Corresponding coarse and fine entries differ by at most $7.49\times10^{-4}$ relative, measured against the fine grid.

For these massive slices the integrand is $c\,p^2+O(p^4)$ near zero, where $c=y^2(h+C/M^2)^2/(64\pi^3)$. An independent integral from zero to $10^{-4}M$ agrees with $c(10^{-4}M)^3/3$ to $1.32\times10^{-8}$ relative. The calculation includes the full infrared interval. A slice at the zero-mass crossing requires its own limiting analysis.

In a checkout containing the two historical study directories, run

```sh
python -I -B local_power.py /path/to/HDblast > reproduced.json
```

The script uses the Python standard library and authenticates the wall archive and two scale-factor references before reading them. Every modeling stratum and both numerical grids appear in [POWER.json](POWER.json). The finite-volume-independent units represent a density rate in the three spatial dimensions of the shell. The result tests a specified open-system energy ledger. It supplies no measurement of a microscopic daughter bath or the origin of the observed Universe.
