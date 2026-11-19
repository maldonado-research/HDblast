# How the initial phase changes the count

The free scalar Hamiltonian is diagonal at the instant of preparation in the chosen oscillator basis. Its mean energy depends on the occupation $n(q)$. A second independent datum, $c(q)=\langle a_{\mathbf q}a_{-\mathbf q}\rangle$, describes the phase correlation of opposite momenta. We use zero first moments and fix the mode phase by taking its initial field amplitude real and positive. All arms of an experiment use exactly the same preparation.

With $F_\pm(\omega)=\int d\eta\,f(\eta)v_q(\eta)e^{\pm i\omega\eta}$ and $K=(\omega^2-q^2)/(8\pi^2)$ above $\omega=q$, define $D_q=\int d\omega\,K(|F_+|^2+|F_-|^2)$ and $Z_q=\int d\omega\,K F_+ F_-$. The pair count in excess of the vacuum on each wall, divided by the squared coupling, is

$$
N_{\rm exc}/y^2=\int dq\,w(q)\left[n(q)D_q+2\,\mathrm{Re}\big(c(q) Z_q\big)\right].
$$

Here $w(q)=3q^2/8^3$ normalizes the uniform band $0\le q\le8$ and $n=1$ throughout it. A common shift of the origin of conformal time cancels in the product $F_+ F_-$. Rotating the initial oscillator convention requires transforming the phase of $c$ as well. The calculation includes the reference vacuum's switched creation channel in each *total*, then subtracts that reference when forming the excess. The vacuum term can vary with the wall.

For a homogeneous pure two-mode squeezed state $|c|=\sqrt{n(n+1)}$. At a single momentum and frequency the total intensity in units of $K y^2$ is $|\sqrt{n+1}\,\overline{F_-}+e^{i\theta}\sqrt n\,F_+|^2$. The matched vacuum intensity is $|F_-|^2$. Positivity of the total is manifest even when a subtracted excess is negative. The independent calculation integrates this square directly. [Independent numerical output](../independent/SQUEEZE.json)

For pulse minus quiet, let $d=\int w n\Delta D_q\,dq$. A single common phase sweeps out the interval $d\pm2\sqrt{2}\,|\int w\Delta Z_q\,dq|$. If the phase can depend on $q$, its sharp radius becomes $2\sqrt{2}\int w|\Delta Z_q|\,dq$. Choosing a phase opposite the argument of $\Delta Z_q$ attains the upper endpoint. The opposite sign attains the lower. Each arm has its own extrema. The code evaluates the explicit states as an independent check on each bound. [Numerical intervals](REFINED.json)

We selected $\theta=-0.5$ radians after inspecting the fine, wide-window excess contrast, then applied that same phase to both windows, both wall meshes and both energy ranges. For the fine wide window the quiet excess falls from $0.29790965\,y^2$ without coherence to $0.06862725\,y^2$ with this phase. The latter preparation has a total of $0.1310363\,y^2$ pairs per mean initial parent from the modeled parent band. Its matched vacuum contributes $0.0624091\,y^2$. The three pulse-minus-quiet excess contrasts reverse sign. For this phase the raw band counts keep their original signs, and converting the excess count to a local number density changes the sign of the distant-positive contrast through the final volume factor. [Transfer calculation](../independent/SQUEEZE.json) · [Primary evolution](coherence.py)

The input energy match holds at preparation. Pressure, the scalar force and later free energy can all depend on correlations. The calculation uses externally supplied wall and mass histories. Determining the phase produced during a higher-dimensional disturbance would require a dynamical preparation model.
