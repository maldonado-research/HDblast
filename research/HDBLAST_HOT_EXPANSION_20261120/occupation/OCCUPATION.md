# When would the omitted modes carry appreciable energy?

This calculation asks what the parity-odd surface pulse would do to **one isolated quantized Maxwell field** on the prescribed homogeneous wall history. Assume the modes are in their conformal vacuum at the entrance to the scalar step, the wall description is valid up to the highest counted wavenumber, and coupling to the bulk and currents can be neglected for these perturbations during the step. These assumptions go beyond the classical evolution. In a fully coupled calculation scattering, plasma effects, fluctuations already present and the response of the scalar and geometry can all modify the answer.

For each physical transverse harmonic on a round unit three-sphere, the curl eigenvalue has magnitude $k=2,3,\ldots$ and degeneracy $k^2-1$ per helicity. Its proper energy per excitation is $k/(aL)$ and the proper spatial volume is $2\pi^2(aL)^3$. The two helicity representations of $SU(2)\times SU(2)$ have spins $(k/2,(k-2)/2)$ and their exchanges, giving $(k+1)(k-1)$ states each. [Closed Einstein universe electromagnetic spectrum, Sec. II](https://arxiv.org/html/1110.4504v1#S2)

Let $M_{kh}$ be the real, unit-determinant transfer matrix in free, unit-energy oscillator coordinates at the two ends of the pulse. Starting from covariance $I/2$ gives final occupation above the outgoing vacuum

$$
n_{kh}=\frac{\operatorname{tr}(M_{kh} M_{kh}^T)-2}{4}
      =\frac{\langle G_{kh}\rangle_{\rm phase}-1}{2}.
$$

A free rotation gives zero. Subtracting the final zero-point energy defines excitation energy at the free endpoint. For a selected set of distinct harmonics $\mathcal K$, the dimensionless energy density in these excitations is

$$
e_{\mathcal K}=\frac{1}{2\pi^2 a_f^4 S_0}
  \sum_{k\in\mathcal K,h=\pm1}(k^2-1) k\, n_{kh},
\qquad S_0=\sigma L^4
$$

in natural units, with wall tension $\sigma$ and length unit $L$. All terms are nonnegative. We use the 18 tested wavenumbers above the homogeneous $k=2$ multiplet. There is no interpolation between sampled integer shells, multiplicity for additional gauge fields, or extrapolation past the tested band. This gives a **lower bound on the free test-field excitation energy** under the stated initial-state and cutoff assumptions. It supplies no lower bound on the energy of a backreacting five-dimensional solution.

At the end of the parity-odd scalar exit, $a_f=8.12451$ and the coherent electromagnetic density is 0.0984662 in wall-tension units. The 18 shells give $e_{\mathcal K}/e_{\rm coherent}=2.9022\times10^{29}/S_0$. Their contribution alone reaches 1% at $S_0=2.9022\times10^{31}$. A smaller action increases their relative importance in this estimate. A larger action suppresses this finite sum while uncounted modes still require a calculation. The selected shell with the largest contribution is $k=8192$. If the field theory ceases to apply at lower momentum, a sum with only the verified modes below that cutoff is required. The released table gives each contribution separately. The same prescribed geometry and kinetic history with the parity-odd coupling turned off give $3.6827\times10^9/S_0$ against the same coherent reference. This control isolates the linear mode response to the added coupling. [Mode transfer](../helicity/HELICITY.json) · [Wall-action convention](https://doi.org/10.5281/zenodo.24406309) · [Finite sum](OCCUPATION.json)

For the two-stage parity-even background the two helicities have equal occupation. Summing the nine sampled shells with $2<k\le1024$ gives $e_{\mathcal K}/e_{\rm coherent}=9.9613\times10^7/S_0$, using the coherent comoving electromagnetic energy 491.6665 at $\tau=6.67$. This separate lower bound uses the combined transfer across both changes in $K$ and the freely evolving phase between them. The restriction to $k\le1024$ avoids estimating very small high-frequency excitation by subtracting near-equal numbers. Numerically we use $n=(G_{\max}-1)^2/(4G_{\max})$, which is equivalent to the covariance expression for a determinant-one matrix. The two lower bounds use different backgrounds and wavenumber selections. Neither computes the total excitation or the nonlinear fate. [Two-stage transfer and endpoint](../helicity/EVEN.json)

The coherent denominator sums the background gauge triad. The numerator counts one field. Three decoupled, identically prepared fields would triple that numerator. The $k=8192$ shell contributes 91.36% of the parity-odd finite sum. Its physical momentum at the entrance is about $1031/L$. Including $k=65536$ requires a valid surface description to about $8245/L$ already at entry. The unknown cutoff and wall action must be specified before using either number to judge a physical model.
