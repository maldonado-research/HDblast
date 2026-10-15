# Wider-state feasibility from the retained analytic box

15 October 2026 UTC. Planning-only exact-rational deduction. The proposed
ell=10^-6 state is feasible for the same linear first-variation method. No
new waveform, response-coefficient interval, finite-amplitude prediction, or
change to the frozen checkpoint is supplied.

## The fixed-box collar test has unused room

Use the authenticated shell intervals and the same outer rectangle

    rho in [129,131], eta in [-9/100000,-8/100000],
    p in [-14/100000,-12/100000], a in [111/1000,112/1000].

The inward ODE is

    rho_z=-rho^2 a, eta_z=-rho p,
    p_z=-rho(U_eta-4ap), a_z=1/rho+rho p^2/3.

The new checker reconstructs the four derivative intervals independently
and exactly reproduces the retained rational intervals. It uses
U_eta=W_eta*(W_etaeta-4W/3), with W=1/3+eta^2+eta^3/3. No reference-only
producer is imported or executed.

For an initial component interval [l_i,u_i], outer interval [a_i,b_i], and
derivative enclosure [f_i,g_i], strict first-exit inclusion is equivalent to

    l_i+L*min(f_i,0)>a_i,
    u_i+L*max(g_i,0)<b_i.

For each negative f_i the lower-face depth threshold is (l_i-a_i)/(-f_i);
for each positive g_i the upper-face threshold is (b_i-u_i)/g_i. Their exact
minimum is

    L_star = 5669686020285051/123672167032900272500
           = 0.000045844478643094814...,

and is attained by the eta upper-face ratio. The other directional ratios
are approximately 0.000480991870426 (rho lower), 0.000146098392826
(p upper), and 0.08031022023164 (a upper). Full rational values are in the
receipt.

Every 0<=L<L_star passes this strict closed-box test. Equality gives zero
eta upper margin and fails that strict test. Thus L_star is its supremal
admissible depth, not an attained maximum or the actual maximal depth of
the regular background. Failing this sufficient test does not establish
an actual background exit or singularity. A convenient larger certified
collar is L=1/25000=4e-5.

For the concrete proposal

    ell=1/1000000, T=3ell, Z=L_energy=8ell=1/125000,

all eight first-exit margins are strictly positive. The smallest margin is
the eta upper margin, about 6.94067738314e-7. Smooth ODE continuation in
the compact rho>0 rectangle extends the retained actual background over
this collar. This is an analytic interval argument, not a validated
numerical background trajectory.

## Coefficients and the pointwise feasibility constants

The outer rectangle has not changed. Therefore the same coarse constants
m=1/1000, M=10^6, B=1000 and Gamma=250 remain valid throughout the longer
collar. The checker independently re-evaluates the derivative-free potential

    mathcal W=-35*tanh(t)^2/4
      +rho^2*(2h^2-10ah+35a^2/4-U_etaeta-p^2/6), h=U_eta/p,

and the true boundary formulas. Its intervals reproduce, to display precision,

    mathcal W in [11698.043620265384,101464.88345514954],
    mu in [.00205977121593046,.00227166178146915],
    beta in [202.071852682969,247.373483657537],
    gamma in [202.086050569103,247.367804503084].

These decimals summarize exact rational endpoints and do not decide gates.
They also independently verify the useful improved rounded constants

    m=1/500, M=102000, B=250, Gamma=250.

For either set of constants, let U=1+T, r=MT^2/2, s=BT^2/(2m),
d=(1-r)(1-s)-Tr/m. The same positive matrix majorant gives

    A_star=U*(1-s+T/m)/d, V_star=U/(m*d),
    eta=ell*(MT*A_star+(1+Gamma*T)*V_star),
    V_tail=ell/m^2+ell/m+(r/m)*A_star+s*V_star,
    eta_tail=ell*(MT*A_star+V_tail+Gamma*T*V_star).

Both determinants are strictly positive. Conditional on the same smooth
linear coupled problem with the newly specified wider initial tangent,
these are candidate rigorous bounds on ||H-Hhat||_infinity and on the
post-pulse |H| respectively, where H=ell*D(ell*x) and
Hhat=Q'(x)+ell*exp(-ell*x)*Q(x).

| Constants | Exact global eta | Exact tail eta_tail |
| --- | --- | --- |
| Original coarse | 2007523995553919/1999981973040500000 | 9538004522513419/1999981973040500000000 |
| Improved rounded | 4005463686978011869/7999991822494065500000 | 14919547389560282619/15999983644988131000000000 |

The coarse values are about 0.001003771045247 and 4.76904524695e-6;
safe bounds are eta<0.001004 and eta_tail<4.77e-6. The improved values are
about 0.000500683472665 and 9.32472665010e-7; safe bounds are
eta<0.000500684 and eta_tail<9.32473e-7. The support-aware L2 error is
respectively below 0.0010038 or 0.000500685. These are new-width error
budgets, not a transfer of the narrow-width leading response coefficient.

## A concrete shape-resolution gate passes

The separate fresh rational calculation in `shape/shape_gate.py` verifies
the unchanged normalized bump geometry and the two endpoint gaps

    M_Q-Q'(1.29)>0.029681983025445975,
    M_Q-Q'(1.32)>0.044578913596622392,
    M_Q=max Q'>4.235640422134888218.

Both gaps exceed 1/50. Thus the conservative sufficient gate

    2*eta+ell<1/50

excludes global positive peaks outside (1.29,1.32), and by reflection
excludes global negative minima outside (1.68,1.71). The original coarse
bound already passes: 2*(0.001004)+10^-6=0.002009<0.02. The improved
bound also passes. Samples at x=5/4 and 7/4 retain opposite signs, with
magnitudes greater than 3.7479 under the coarse error budget. A convenient
coarse amplitude summary is

    4.2346<max H<4.2367, -4.2367<min H<-4.2346.

These certify feasibility of resolving the dominant signed shape. They do
not identify which lobe attains the absolute maximum, or prove unique
actual extrema/zeros. The physical scalar would have reversed signs because
the retained gain remains negative. No new continuous waveform export,
integrated response coefficient, or physical-unit calibration is made here.

## Assumptions and reproducibility

The retained analytic branch theorem, exact shell box, inward ODE and smooth
continuation remain inputs. For an exact nonlinear-family interpretation,
the localized inverse, reflection, gauge, whole-correction support control,
and finite-regularity evolution hypotheses must apply to the newly specified
wider tangent; the narrow-state result alone does not establish that fact.
The same proposed causal geometry still has strict margins:

    Z-(L_corr+c_*T_evol)=5ell/8>0,
    T_evol-T=ell/2>0,

with L_corr=3ell, c_*=5/4, T_evol=7ell/2. No effective nonlinear epsilon
radius is computed. First-exit inclusion is only a background estimate and
does not replace these separate family hypotheses.

Both mathematical inputs are SHA-256 authenticated by the checker. It uses
only Python's standard library, performs every decision with exact rational
arithmetic, and retains its checks under python -O. The shape verifier is
self-contained. To replay after relocating the planning package, provide
an extracted frozen waveform capsule containing `inputs/previous/archive`:

```sh
python -I -B check_feasibility.py --source-root /absolute/path/to/extracted/capsule
python -I -O -B check_feasibility.py --source-root /absolute/path/to/extracted/capsule
python -I -B shape/shape_gate.py
python -I -O -B shape/shape_gate.py
```

The default source-root supports the original task layout only; the explicit
option removes that dependency. Output is JSON on stdout. No archived or
reference-only producer is executed, and no frozen scientific file is edited.
