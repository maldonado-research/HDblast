# How negative must the wall pressure be for expansion to accelerate?

A radial light pulse presses on a wall even if it reflects completely and leaves no energy behind. We can test the sign of the wall acceleration directly from the bulk normal pressure. The calculation uses the five-dimensional Einstein equation and the Israel junction conditions. The role of bulk stress and flux in cosmological wall evolution has long been established. Here we apply these laws to the two ordinary retained interiors used in this study. [Binétruy, Deffayet and Langlois](https://arxiv.org/abs/hep-th/9905012) · [Jennings et al.](https://arxiv.org/abs/hep-th/0412281)

## Normal pressure in a spherical spacetime

In units of AdS length one, write the metric as $h_{ab}dx^a dx^b+r^2d\Omega_3^2$ and let $D$ denote the connection of the two-dimensional Lorentzian metric $h$. Define the local spherical mass parameter $m=r^2[1+r^2-(Dr)^2]$. In an orthonormal frame consisting of the wall velocity $U$ and spacelike normal $n$, the radial part of the Einstein tensor is

$$G_{ab}=-\frac{3}{r}D_aD_b r+
\left[\frac{3}{r}D^2r-\frac{3}{r^2}\{1-(Dr)^2\}\right]h_{ab}.$$

With $D^2r=-D^2_{UU}r+D^2_{nn}r$, the normal projection cancels the normal Hessian. Substitution of $G_{nn}-6=\kappa_5^2T_{nn}$ gives

$$ \frac{D^2_{UU}r}{r}=-1-\frac{m}{r^4}-\frac{\kappa_5^2}{3}T_{nn}. $$

Here $D^2_{UU}r$ denotes the Hessian on two copies of $U$. On the wall $r=a(\tau)$, its proper-time acceleration also includes $(D_UU)\cdot Dr$. Choose each side's normal from the retained interior toward the wall, and let $b_i=n_i\cdot Dr/a>0$, $Q_i=m_i/a^4+\kappa_5^2T_{nn,i}/3$. Since $D_UU=K^\tau{}_{\tau,i} n_i$, the result on either side is

$$ K^\tau{}_{\tau,i}=(\ddot a/a+1+Q_i)/b_i. $$

## The pressure threshold

Write the extra surface energy and pressure in units of the tuned tension as $e=\rho/\lambda$ and $p=P/\lambda$. The angular and temporal junction conditions for an isotropic wall give $(b_1+b_2)/2=1+e$ and $(K^\tau{}_{\tau,1}+K^\tau{}_{\tau,2})/2=1-2e-3p$. With $B_h=2/(1/b_1+1/b_2)$ and $\bar Q=(Q_1/b_1+Q_2/b_2)/(1/b_1+1/b_2)$, one obtains

$$\frac{\ddot a}{a}=(1-2e-3p)B_h-1-\bar Q.$$

At specified local data the exact acceleration threshold is $p<[1-2e-(1+\bar Q)/B_h]/3$. The harmonic mean satisfies $B_h\le 1+e$. If $\bar Q\ge0$, positive acceleration requires $1-2e-3p>0$ and hence

$$ p<-\frac{e+2e^2+\bar Q}{3(1+e)}. $$

For reflection-symmetric geometry this last inequality is also sufficient for positive local acceleration. For unequal $b_i$ it is only necessary. For instance, at rest with $a=2$, $e=0.05$, $b_1=0.99$, $b_2=1.11$ and $\kappa_5^2 T_{nn,i}/3=(0.05,0.08)$, the necessary ceiling is $-0.08579$ while the exact threshold is $-0.08706$. Accelerating expansion and a local bounce have different initial data. A smooth local bounce additionally has $\dot a=0$ and passes from contraction to expansion.

For ordinary wall radiation $p=e/3$ with $e\ge0$, nonnegative masses and nonnegative bulk normal pressures imply $\ddot a\le0$. The sign survives one-sided irradiation, transmission or reflection whenever these assumptions continue to hold. Take incoming and returning radial null momenta proportional to $U+n$ and $U-n$ in a wall rest frame. For nonnegative wall-frame intensities $I,R$, their contribution to $T_{nn}$ is $I+R$ and the net outward energy current $-T_{Un}$ is $I-R$. Equal beams double the pressure relative to the incident beam and cancel the net current. This local argument leaves the reflection law, the wall’s microscopic response and the time evolution of the two streams unspecified.

## Matched example at rest

Join identical massless vacuum AdS interiors. With $a(0)=2$, $\dot a(0)=0$ the angular condition fixes $e=\sqrt{1.25}-1\simeq0.118034$. The critical extra pressure is $-0.0434984$. An additional positive vacuum surface energy with $p=-e$ obeys wall conservation at constant $e$. Its exact solution $a(\tau)=2\cosh(\tau/2)$ contracts and bounces, with $\ddot a/a=1/4$. Ordinary wall radiation with the same initial density and radius has $\ddot a/a\simeq-0.277864$ at rest and begins contracting. The vacuum energy is present in the initial state. A [separate calculation](SCALAR_CAPTURE.md) supplies a candidate mechanism by transferring the pulse energy to a canonical scalar.

The replay computes the full five-dimensional Christoffel symbols and Ricci tensor for an independently specified double-null metric at two points and two observer velocities. The resulting effective stress can have negative normal pressure. These cases check the geometric identity for general stress. A separate check in retarded Vaidya recovers the analytic positive null pressure for both inward and outward wall velocities. The replay also checks the two-stream flux and pressure algebra, the unequal-side example and both matched equations of state. [Computation](normal_pressure.py) · [Numerical results](JUNCTION.json)
