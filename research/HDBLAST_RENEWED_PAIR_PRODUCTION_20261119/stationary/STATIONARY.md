# What can switching alone produce?

We held the scale factor and each parent oscillator's proper frequency fixed at the beginning of the late gate. We then applied an external pulse with the same duration between its mapped endpoints and the same unit-peak sine-to-the-fourth profile in proper time. On the quiet wall this reproduces precisely the prescribed gate as a function of proper time. On other walls this control uses their actual elapsed duration and drops the small internal time distortion of the scalar clock. In a stationary vacuum with the interaction off before and after the pulse, every emitted pair and its accompanying parent excitation draw energy from the switch. This isolates one possible source of particles in a finite-duration calculation. [Registered control](https://github.com/maldonado-research/HDblast/issues/157#issuecomment-7017450151)

| Frozen quiet-shell gate | Pairs created from the vacuum | Proper energy in the two daughters | Proper energy in the parent | Work supplied by the switch |
|---|---:|---:|---:|---:|
| Full second window | 0.03647 | 2.19568 | 1.83473 | 4.03041 |
| First half at unit peak | 0.20190 | 28.40695 | 10.15708 | 38.56404 |
| Later half at unit peak | 0.19509 | 27.36979 | 9.99111 | 37.36090 |

Values are leading coefficients per normalized comoving volume through daughter frequency 1920. They carry the same factor $y^2$ as the other pair calculations. Shortening the primary pulse by a factor of two at the same starting wall time increases the stationary-vacuum switching work by a factor of 9.57. A short switch can produce rare, energetic states. The table concerns a deliberately frozen control. The native rolling shell has an excited, phase-sensitive parent and exchanges energy with its scalar field.

For a proper gate duration $T$, write $G_T(\nu)=\int_0^T dt\,\sin^4(\pi t/T)e^{i\nu t}$. If $x=\nu T$, its centered real Fourier amplitude is
$$
 e^{-ix/2}G_T(\nu)/T=
 \frac{48\pi^4\sin(x/2)}{x(x^2-4\pi^2)(x^2-16\pi^2)}.
$$
The apparent poles have removable limits, evaluated with shifted sinc functions. This transform also agrees with a direct quadrature of the proper-time gate. For frozen parent energy $E_q$ and canonical factor $I=a^3 E_q$, the two amplitudes are $D=G_T(\omega/a+E_q)/\sqrt{2I}$ and $A=G_T(\omega/a-E_q)/\sqrt{2I}$. The vacuum contribution is $|D|^2$. When the retained covariance is used instead, the ordered probability includes $n|A|^2+(n+1)|D|^2+2\operatorname{Re}(D^* A\langle bb\rangle)$. We propagate the calibrated free covariance up to switch-on to get $n$ and $\langle bb\rangle$. For the full late gate, the frozen retained-parent calculation gives 45.21 proper daughter-energy units through 1920, including 38.73 in the parent-annihilation term. Vacuum-only creation is therefore insufficient to describe the full source.

For vacuum creation the integrand is positive. Each transition costs $\omega/a+E_q$, the sum of final parent and daughter energies. Multiplying by that energy and summing the spin-summed massless pair measure gives the work in the table. For $|x|\ge8\pi$ the transform obeys $|G_T|\le 48\pi^4(16/15)(4/3)T/|x|^5$. Integrating this envelope over the entire tail bounds the missing quiet work above cutoff 1920 by $1.30\times10^{-7}$ for the full gate, and by $3.43\times10^{-5}$ for either half gate. The analytic tail bound here applies to the stationary vacuum. The all-frequency claim for the *rolling* source uses its own Fourier-norm bound. [Frozen calculation](stationary.py) · [Twelve gate and wall combinations](STATIONARY.json)
