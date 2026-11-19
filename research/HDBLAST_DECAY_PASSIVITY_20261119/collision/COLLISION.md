# Work in a finite collision with a squeezed oscillator

This two-mode example realizes attenuation toward a squeezed vacuum and measures the energy of a physical ancillary oscillator. Both oscillators have the same bare quantum energy $E$. The ancilla begins in its physical vacuum. The interaction vanishes at the endpoints, so the increase of the sum of bare energies measures the net work supplied by the external control.

Let $A$ annihilate the ground state of the physical system Hamiltonian and $C$ annihilate the ancilla vacuum. Define $B=A\cosh r+A^\dagger\sinh r$, with $\sinh^2r=1/4$. Apply

$$
U(\theta)=\exp[\theta (B C^\dagger-B^\dagger C)],\qquad \eta=\cos^2\theta=0.64.
$$

In the infinite oscillator space this transformation gives $U^\dagger C U=C\cos\theta+B\sin\theta$. The classical control can exchange energy with the pair through the pair creation and annihilation terms in the generator. Work is positive when the control supplies energy. Its division among a pump and individual switching operations would require the drive waveform.

| Initial system state | $n_A$ | $n_B$ | System energy change / $E$ | Ancilla energy change / $E$ | Control work / $E$ |
|---|---:|---:|---:|---:|---:|
| Physical vacuum | 0 | 0.25 | +0.09 | +0.09 | +0.18 |
| Opposite squeezing phase in the target basis | 1.25 | 0.25 | −0.36 | +0.09 | −0.27 |
| Target vacuum | 0.25 | 0 | 0 | 0 | 0 |

The first two states deliver the same energy to the ancilla. One requires positive net work, while the other allows energy to flow back into the external control. The third state is stationary. These results are for this explicitly driven collision with a physical vacuum ancilla.

## Numerical test

The calculation starts in the physical $A,C$ number basis, projects and normalizes the analytic squeezed-vacuum coefficients there, directly constructs the sparse generator and applies its matrix exponential to each state. The three system preparations are $|0_A\rangle$, $S(2r)|0_A\rangle$ and $S(r)|0_A\rangle$, where $S(r)=\exp[r(A^2-A^{\dagger 2})/2]$. No covariance evolution enters the exponentiation.

With 16, 32, 64, 96 and 128 levels per oscillator, the largest difference between any computed energy change or net work and the infinite-oscillator value is $1.71\times10^{-2}$, $2.28\times10^{-4}$, $2.73\times10^{-8}$, $2.79\times10^{-12}$ and $1.42\times10^{-15}$ in units of $E$. A second readout integrates $\langle[H_A+H_C, B C^\dagger-B^\dagger C]\rangle=2E\sinh r\langle A^\dagger C^\dagger+AC\rangle$ over the unitary path using eight Gauss–Legendre nodes. The two work readouts agree within $2.1\times10^{-15}E$ at the largest dimension. [Numerical results](COLLISION.json)

Reproduce with Python, NumPy and SciPy. The optional `--eta` argument changes the transmissivity and `--dimensions` selects the Fock cutoffs.

~~~sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python -I -B collision.py
~~~
