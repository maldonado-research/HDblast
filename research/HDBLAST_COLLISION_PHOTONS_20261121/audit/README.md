# A closed-spectrum check on collision radiation

The spinor calculation integrates amplitudes at selected incoming pair masses and interpolates across the event population. We independently evaluated analytic photon spectra at the actual mass of every recorded encounter in all eight elastic simulations. The total energies from the two calculations agree within $4.3\times10^{-9}$ in relative terms. The largest fractional difference for any of their five energy bands is $1.3\times10^{-7}$. At the fine 7.05 endpoint the direct result is **$1.101101848\times10^{-14}$ of the final gas energy**, and 77.88% comes from radiation off the incoming charged legs. The top 1% of encounters by incoming invariant mass contribute 12.46% of that emitted energy. They carry 3.69% of the energy summed over incoming pairs. The largest single contribution amounts to 0.034% of the predicted emission. These two computations share the assumed interaction, elastic histories and daughter charge.

## Analytic spectra

Write $s$ for the incoming invariant mass squared in daughter-mass units, $x=2E_\gamma^*/\sqrt{s}$ for the pair-frame photon fraction, $q=s(1-x)$, and $b_z=\sqrt{1-4/z}$. Define $F(s,x)=\alpha^{-1}dN_{\rm FSR}/dx$. The angle-integrated vector result for radiation from the outgoing unit-mass fermion pair is

$$
xF(s,x)=\frac{[(1+b_s^2)(3-b_s^2)-2x(3-b_s^2)+2x^2]\log[(1+b_q)/(1-b_q)]
 -2[(3-b_s^2)(1-x)+x^2]b_q}{\pi b_s(3-b_s^2)}.
$$

This is Eq. (5) of Bystritskiy and collaborators divided by their Born cross section and by $\alpha$. Their Born normalization uses a pointlike vector source with light incoming electrons. After the complete angular average, the tensor of the opposite incoming pair is proportional to the projector transverse to the neutral vector's momentum. Its mass-dependent overall normalization cancels in this final-state radiation ratio. [Bystritskiy et al., equations (5) and (8)](https://arxiv.org/abs/hep-ph/0505236)

To obtain the incoming-leg result, exchange the vector virtualities before and after emission and reverse the outgoing photon momentum in the charged-pair radiation tensor. The arguments of $F$ become $(q,-sx/q)$. The spin-averaged elastic matrix element, averaged over angles, is proportional to

$$
B(s)=\frac{(s+2)^2}{(s-8300^2)^2+0.30^2(s+2)^2 b_s^2}.
$$

The spin and transverse-vector projectors for the crossed radiation tensor then match those of the corresponding outgoing tensor. The phase-space factor for an emitted photon at fixed $s,x$ is proportional to $s x b_q/b_s$ before division by $B(s)$. For the crossed arguments the corresponding factor is $-s x b_s/b_q$, before division by $B(q)$. The differential initial-state spectrum follows,

$$
I(s,x)=-\frac{b_q^2}{b_s^2}\frac{B(q)}{B(s)}
          F\!\left(q,-\frac{s x}{q}\right).
$$

The negative photon fraction is an analytic continuation of the kernel. All emitted photons in the integrals have $0<x<1-4/s$. For a highly relativistic incoming pair this expression also approaches the familiar initial-state splitting function $[1+(1-x)^2](\log s-1)/(\pi x)$ times the recoil-to-original vector Born cross-section ratio at the appropriate virtualities. The calculation retains the finite masses when applying the exact radiator.

The elastic angular factor implicit in $B(s)$ is the same $3-2b_s^2+b_s^4\cos^2\theta$ used to generate the event histories. With the full final-state angles integrated, the initial–final interference is odd under exchanging the outgoing charges and drops out of the inclusive photon-energy moment. Exchanging the incoming beams makes the remaining pair-frame photon direction symmetric, so the inclusive energy transforms to the wall by multiplying by the Lorentz factor of the incoming pair. Conditioning a photon on a selected outgoing charge or a restricted angular acceptance changes that argument. The benchmark includes the registered neutral $s$-channel amplitude. A complete same-flavor vector theory can also exchange a mediator in the crossed elastic $t$ channel and interfere with the annihilation amplitude. Including it would change the underlying event rates and angular distribution. [BABAR charge-asymmetry measurement](https://arxiv.org/abs/1508.04008) · [Registered elastic interaction](https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_SCATTERING_20261121/live)

## Checks and reproduction

The code integrates the energy moment $(x/2)[I+F]$ from zero to the two-body massive recoil endpoint. For positive lower cuts it integrates the photon number spectrum as well. During hard recoil it resolves the return to the neutral pole in $\sqrt{\log(q/4)}$, splitting the integration at $q=M^2/2,M^2,2M^2$ when accessible. It uses 48 and 24 Gauss nodes per segment, and sums contributions directly at every measured $s$. The relative change in integrated energy between these orders stays below $6.2\times10^{-7}$ across the eight ensembles. [Generator](radiators.py) · [Results](RADIATORS.json)

At the largest incoming mass in the fine late run, $s=1.46782\times10^{11}$, the mean pair-energy fraction emitted per unit $\alpha$ is 74.653070595 with 96 quadrature nodes. Doubling from 48 to 96 changes it by less than $4\times10^{-11}$ in absolute units. The direct spinor integration gives 74.653031806, a fractional difference of $5.2\times10^{-7}$. Checks of the separated incoming and outgoing differential spectra cover threshold, resonance, and high-energy kinematics, with photon cuts as close as 99.99% of the massive endpoint.

~~~sh
python -I -B audit/radiators.py --hard hard --artifacts "$results/events" \
    --reference "$results/HARD.json" --output "$results/RADIATORS.json"
~~~
