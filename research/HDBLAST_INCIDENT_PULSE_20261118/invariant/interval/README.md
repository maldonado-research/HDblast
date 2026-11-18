# A certified nonzero gluing amplitude

On the fixed initial slice prescribe the quiet areal profile $A=A_0$ and $A_T=1$, and let the other metric coefficients respond. For the specified compact scalar perturbation, exactly one amplitude in $0<\eta\le0.003$ matches the quiet geometry at both ends. The certified amplitude lies in **(0.00201, 0.00203)**. These statements use the declared background, slice, profile and wall normalization. They concern the initial constraints. Subsequent energy transfer requires evolving the data.

## The endpoint whose sign is certified

Write $t=-Z$, $E=\exp B_0$, $q=A_{0,Z}$, $r=\phi_{0,Z}$. The bump is $f=\exp[1-1/(1-x^2)]$ on $|x|<1$, zero elsewhere, with $x=6-50t$. Put $g=f_Z$, $v=\eta g$ and $h=\eta f$. The scalar data are $\phi=\phi_0+h, p=\phi_T=-v$. Define the planar mass defect $D=\mu-\mu_0$, with $\mu_0=E^2(1-q^2)$. Over the pulse it satisfies

```text
D_t = a D + b
a = −[q(v²+(r+v)²)+2v(r+v)]/[3(1−q²)]
b = −(2/3)E²(q+1)(rv+v²) − (2/3)E⁴q[U(φ₀+h)−U(φ₀)].
```

The bump is supported on $0.10\le t\le0.14$. At its first edge $D=0$. The endpoint function is $F(\eta)=D(0.14)$. As the linear mass equation depends smoothly on amplitude, $F$ is continuous and twice differentiable. When $F=0$, uniqueness of the homogeneous mass-defect equation beyond the support gives exact quiet data there. Flatness of the bump at each edge gives matching derivatives. The code verifies $q>2.7418$ and $\mu<-191.74$ up to the far pulse edge throughout the enclosing amplitude interval $[.001,.003]$, so the metric reconstruction stays on its real, regular branch.

| Amplitude | Outward enclosure for $F$ |
|---:|---:|
| 0.001 | [0.01252, 0.01263] |
| 0.00201 | [0.0001660, 0.0003764] |
| 0.00203 | [−0.0003324, −0.0001199] |
| 0.003 | [−0.03634, −0.03601] |

For the whole interval $[0,.003]$ the program also encloses the second parameter derivative by $-24848<F''<-24416$. It integrates $V=D_\eta$, $S=D_{\eta\eta}$ alongside $D$, with initially zero derivatives:

```text
V_t = aV + a_η D + b_η
S_t = aS + 2a_η V + a_ηη D + b_ηη.
```

The endpoint signs give existence by the intermediate value theorem. Strict concavity and $F(0)=0$ make $F(\eta)/\eta$ strictly decreasing for positive amplitude. The nonzero root is unique on the declared interval and has negative derivative since the earlier endpoint value is positive.

## A scalar test of physical distinctness

For every $\eta\in[.001,.003]$, the value of $\phi$ at the pulsed event $Z=-.13$ lies strictly inside the quiet range on the interval $-.13\le Z\le-.129$. The proof checks quiet values at $t=.13$ and at the first mesh point at or above $t=.129$. It verifies $r>2.788$ throughout that comparison interval. On this monotone quiet branch there is one point with the same scalar value. Its invariant $X=(\nabla\phi)^2$ satisfies $X<0.282892$. At the pulsed event the invariant obeys $X>0.290630$, using

```text
X_pulse = (1+D/μ₀)(r²+2r ηg)/E²,     X_quiet = r²/E².
```

The invariant differs at equal field value for all possible roots in the bracket. The conclusion compares to this quiet branch with its declared normalization.

## Validated integration and replay

The background starts at the matter-free wall using the decimal inputs $\phi_b=0.0006170897672752272$, $\delta=0.1$, $c=0.5975949350280132$, $d=-3.106933495673783$. Define $W=1-\phi+\phi^3/3$, $U=W_\phi^2/2-2W^2/3$, and $\sigma=2W+\delta(1+c\phi+d\phi^2/2)$. Set $E_b=(\sigma^2/36-\sigma_\phi^2/48+U/6)^{-1/2}$, $q_b=E_b\sigma/6$, $r_b=-E_b\sigma_\phi/2$. This enforces the quiet static first integral at the wall prior to interval evaluation. The program uses second-order Taylor endpoint enclosures and verifies a Picard box for each entire step of

```text
E_t=−Eq,   q_t=−1+q²+r²/6+E²U/3,
φ₀,t=−r,   r_t=−E²U_φ+3qr.
```

It then uses whole-step interval coefficients and Picard boxes for the mass and both parameter derivatives. The mean-value integral over each step bounds its endpoint update. The potential difference is evaluated by its factored polynomial to retain small differences. Each elementary arithmetic interval is rounded outward by one adjacent binary64 number. The square root uses IEEE 754 arithmetic. The only exponential uses a bounded alternating series after range reduction, with the bump derivative maximum handled analytically near the edges. No nonvalidated ODE solver or platform exponential library enters this calculation.

Compile in IEEE binary64, round-to-nearest mode, with floating-point contraction and fast-math disabled. The default run uses 16,384 steps per 0.01 of radial distance and takes a few seconds. A failed check returns a nonzero exit status. The lines labeled defect[0,1,2,3,4] in the output correspond to 0.001, 0.003, the full interval [0.001,0.003], 0.00201 and 0.00203, respectively.

```sh
c++ -std=c++17 -O2 -fno-fast-math -ffp-contract=off bracket.cpp -o bracket
./bracket
```

[Replay source](bracket.cpp) · [Recorded interval endpoints](BRACKET.txt)
