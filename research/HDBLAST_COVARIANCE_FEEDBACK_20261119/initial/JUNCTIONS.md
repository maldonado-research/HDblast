# What can the retained state match at the intervention?

At the historical proper time $s_i=1.30371703$ on the fine wall ($T=1.11$), the retained scalar has positive field–momentum covariance at all 24 momentum nodes. Keeping the two diagonal covariances fixed leaves essentially no freedom to change the mixed covariance while preserving the first time derivative of the scalar force. The reflected and erased preparations share the same energy, pressure, force and instantaneous junctions. Their subsequent evolution includes a change in the time derivatives of the boundary data at this corner. [Numerical diagnostics](CORNER.json)

## First derivatives of the junctions

Use a canonical real-mode covariance with entries $X=\langle\chi^2\rangle$, $P=\langle\Pi^2\rangle$ and $C=\langle\{\chi,\Pi\}/2\rangle$, with $\Pi=a^3d\chi/ds$ in the native normalization. Primes here denote derivatives in $s=H_0\tau$. On the closed free parent sector, $X'=2C/a^3$, $P'=-2a^3 E_q^2 C$, and $C'=P/a^3-a^3E_q^2X$. The physical frequency satisfies $E_q^2=q^2/a^2+G^2(\phi_b-\phi_*)^2$. Write $\Delta C=C_{\mathrm{new}}-C_{\mathrm{ret}}$ and $M_j=\int q^{2j}\Delta C\,q^2dq/(2\pi^2)$ over the parent band. At a fork with $X$ and $P$ equal modewise, the instantaneous subtraction is common. At that instant

$$
\Delta \rho'=0,\qquad
\Delta p'=-\frac{2}{a^3}\left(m^2 M_0+\frac{2M_1}{3a^2}\right),
\qquad
\Delta j'=\frac{2G^2(\phi_b-\phi_*)}{a^3}\,M_0.
$$

Away from the mass crossing, matching $p'$ and $j'$ is equivalent to $M_0=M_1=0$. In the bulk conformal coordinates, differentiating the signed junctions at fixed initial endpoint geometry gives $\Delta B_{TZ}=-(b/2)(ds/dT)^2\Delta p'$ and $\Delta\phi_{TZ}=-(b/2)(ds/dT)^2\Delta j'$, where $b=6.49\times10^{-6}$. The $A_{TZ}$ junction stays equal since $\Delta\rho'=0$. Reflecting $C$ gives, on the fine mesh, $\Delta p'=126.82957$ and $\Delta j'=574.77734$ in the native hatted variables, and changes the required coordinate velocity slopes by $-0.00087318$ and $-0.00395713$. Erasing $C$ halves these values. An exactly smooth continuation of the retained bulk data across the intervention would require the differentiated boundary data to agree as well.

## An obstruction from purity

Under closed Hamiltonian evolution $XP-C^2=1/4$ for a state initially pure in these variables. At fixed $X,P$, the uncertainty inequality for a new covariance gives $|C_{\mathrm{new}}|\le\sqrt{XP-1/4}=|C_{\mathrm{ret}}|$. With $C_{\mathrm{ret}}>0$, every admissible $\Delta C$ is nonpositive. The positive quadrature weights give $M_0=0$ only when every $\Delta C$ is zero. In a continuous band the analogous implication holds almost everywhere under the stated sign hypothesis. More generally, at most one sign change of the retained $C$ across nonzero momentum nodes suffices when both derivatives are matched. Choose a real affine function of $q^2$ with the sign of the retained $C$ on those nodes. Its product with $\Delta C$ is nonpositive; its weighted integral vanishes by the two moment constraints. For nodes with $C=0$, purity already forces $C_{\mathrm{new}}=0$.

| Retained fork | Smallest $C$ on native nodes | Maximum absolute purity error | Upper bound on the relative change at matched $j'$ |
|---|---:|---:|---:|
| Coarse | 0.00309858 | $4.31\times10^{-11}$ | $8.63\times10^{-7}$ |
| Fine | 0.00307982 | $1.42\times10^{-12}$ | $2.83\times10^{-8}$ |

The last column accounts for the small numerical determinant drift. At exactly matched $j'$, upward and downward changes in $C$ have equal total weight. Thus $\int |\Delta C|\,d\nu\le 2\int[\sqrt{XP-1/4}-C]_+\,d\nu$. The table divides this bound by the size $2\int C\,d\nu$ of a full reflection. An independent mode evolution using the archived quiet wall reproduced the 24 fine-node covariances to $4.56\times10^{-9}$ maximum absolute difference. Across 257 equally spaced momenta including both endpoints the smallest sampled fine $C$ was 0.00306766. Sampling tests the sign across the band. The verified obstruction applies to the positive quadrature nodes of the native experiment. [Historical fork covariances](FORK.json) · [Independent mode solver](corner.py)

In the instantaneous oscillator basis the anomalous correlation is $z=R+iC$, with $R=(a^3 E_q X-P/(a^3 E_q))/2$. On a known state, an equal number rotation of the opposite momenta by $\arg z$ sends $z$ to $e^{-2i\arg z}z=z^*$. The rotation commutes with the instantaneous free Hamiltonian. A classical random choice between this rotation and the identity can produce any intermediate mixed covariance while keeping the diagonals. Its mixture need not be Gaussian. A Gaussian replacement with the same second moments gives the same free stress and leading-order pair readout used here. Such rotations describe how an imposed preparation can preserve instantaneous energy. Establishing a physical mechanism on the shell for these momentum-dependent operations requires dynamics for the intervention.

## Two mixed preparations with equal differentiated sources

Two newly prepared states can share the first source derivatives with one another even on this cut. Put $t=(q/8)^2$, choose $h(t)=t^2-b_1t-b_0$ orthogonal to $1,t$ under the positive weight $C_{\rm ret}\,d\nu$, and set $C_\pm=\pm0.9 C_{\rm ret}h/\max_{[0,1]}|h|$. Their shared $X,P$ and their vanishing $M_0,M_1$ relative to $C=0$ give the same instantaneous stress, force and first derivatives as the erased state. Both differ from the retained incoming derivative. The factor 0.9 puts the branches inside the allowed covariance interval.

With $M_j$ formed from $C_+-C_-$, generally $M_2\ne0$. To find the first separation on a shared smooth background, integrate the three covariance equations against $q^{2j}d\nu$. For example $P_j'=-2a^3m^2 C_j-2a C_{j+1}$. Equality of all initial $X_j,P_j$ and the two matched $C_j$ gives

$$
\Delta p=\Delta p'=\Delta p''=0,\qquad
\Delta p'''=\frac{16M_2}{3a^7},\qquad
\Delta j^{(k)}=0\ (0\le k\le 4),\qquad
\Delta j^{(5)}=\frac{32G^2(\phi_b-\phi_*)M_2}{a^7}.
$$

On the fine wall $M_2=10.26050$, $\Delta p'''=140.63286$ and $\Delta j^{(5)}=-1.75373\times10^6$. A direct time integration of the covariance differences on the common quiet history at $\Delta s=0.00075$ gives 0.9963 and 0.9939 times the predicted leading pressure and force changes, respectively. Halving $\Delta s$ from 0.0015 moved both ratios closer to one. This control compares two newly prepared ensembles on a prescribed history, independently of the native bulk evolution of the retained, reflected and erased branches.

To reproduce both meshes from this directory, run the following command against HDblast containing the published, hash-pinned closed-parent wall archive. The native checkpoints used to extract the initial covariances for this check are identified by the file digests in FORK.json.

~~~sh
OPENBLAS_NUM_THREADS=1 python -I -B corner.py /path/to/HDblast > rerun.json
~~~
