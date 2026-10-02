# Matched absolute source for a smooth flat-space crossing

This applies standard quantum field theory to a controlled prerequisite for the HDBLAST matter calculation. It is not an external mathematical novelty claim or a new cosmological solution.

## One operator and one covariant prescription

Let x=M²=4tanh²(t), y_j=x+jLambda², c_j=(1,-3,3,-1). Each minimally coupled real field solves fddot+(k²+y_j)f=0 with Wronskian i. The PV identities Sum c_j=Sum c_j j=Sum c_j j²=0 cancel the requisite ultraviolet moments. Regulators are auxiliary mathematical subtractions, not physical negative-energy particles, and Lambda is not the five-dimensional gravity cutoff.

The common one-loop local action has -V(x)+F(x)R. In the adopted convention:
V_PV=Sum c_j y_j²ln(y_j/mu²)/(64pi²),
F_PV=Sum c_j y_j[ln(y_j/mu²)-1]/(192pi²).
The same scale mu cancels out of the PV sums. Subtract the quadratic Taylor polynomial C_V of V_PV at r=4, and the linear polynomial C_F of F_PV there. This fixes V,Vx,Vxx,F,Fx=0 at that reference. These are chosen finite EFT matching conditions, not facts inferred from the hypothesis.

## Why flat pressure needs curvature matching

Varying the curvature term gives
T_munu[FR]=-2[F G_munu+(g_munu box-nabla_mu nabla_nu)F].
Even when R=0, homogeneous F(t) gives rho_F=0 and p_F=2 Fddot. Therefore subtracting C_F R adds -2F'_PV(r)xddot to pressure.

The high-momentum dynamical mode pressure contains -xddot/(24k³), whose integral diverges logarithmically. Vacuum-potential matching alone cannot remove it. An energy ledger rho_dot=xdot*S at H=0 contains no pressure and cannot diagnose this failure.

The independently derived second-order ultraviolet terms are
q2=xddot/(16omega^5)-5xdot²/(64omega^7),
e2=xdot²/(64omega^5).
These provide large-positive-momentum checks, not a heavy-mass replacement at x=0.

## Stable mode observables and analytic vacuum restoration

For each mode e=(|fdot|²+omega²|f|²)/2,
p=(|fdot|²-(k²/3+y)|f|²)/2, q=|f|².
Define E_j,P_j,Q_j by integrating e-omega/2, p-k²/(6omega), q-1/(2omega), respectively, against k²dk/(2pi²). Then
rho=Sum c_j E_j+V_m,
p=Sum c_j P_j-V_m-2F'_PV(r)xddot,
S=Sum c_j Q_j/2+V_m',
where V_m=V_PV-C_V. The external mass source obeys rho_dot=xdot*S. For phi=tanh(t), J_phi=8phi*S.

A hard three-momentum cutoff on bare vacuum pieces gives
rho0_PV(K)+p0_PV(K)=K³ Sum c_j sqrt(K²+y_j)/(12pi²),
with leading tail -Lambda^6/(32pi²K²). This can grow along K proportional to Lambda. Numerically integrate the dynamic differences and restore the exact continuum static potential instead. The implementation uses an equivalent constant-reference form to keep the infrared cancellation smooth through x=0.

## Regulator limit and a second pressure representation

For fixed r>0:
V_R=[x²ln(x/r)-1.5x²+2rx-.5r²]/(64pi²),
V_R'=[xln(x/r)-x+r]/(32pi²).
Use the continuous xlnx=0 limit at zero.

Heavy-regulator dynamical energy/current decouple in the formal limit, so rho_R=E_physical+V_R and Q_R=Q_physical+2V_R'. Pressure has the convergent representation
p_R=integral k²dk/(2pi²)[p_physical-k²/(6omega)+xddot(2k²/3+r)/(16(k²+r)^(5/2))]-V_R.
The subtraction uses the fixed positive reference mass. Its analytic integral through K is
xddot/[48pi²] [asinh(K/sqrt(r))-v+v³/6], v=K/sqrt(K²+r).
At large K this becomes xddot[ln(2K/sqrt(r))-5/6]/(48pi²). The constant is fixed by the stated curvature matching.

The formal local curvature remainder F_R=[xln(x/r)-x+r]/(192pi²) has a singular derivative expansion at x=0. Do not add 2F_Rddot there. Exact modes and polynomial counterterms evaluated at positive r avoid that invalid approximation.

## A pressure test independent of the energy ledger

Let Q_R=2S and
A_Lambda=-4V_m+2xV_m'-Lambda² Sum_(j>=1)c_j j Q_j.
Then -rho+3p=Q_Rddot/2-xQ_R+A_Lambda after momentum cutoff removal. Evaluate Q_Rddot from resolved finite differences rather than substituting the mode equations. At finite K add
D_K=1/2 d²/dt²[Q0_PV(K)-2V'_PV(x)],
where Q0_j(K)=[Ksqrt(K²+y_j)-y_j asinh(K/sqrt(y_j))]/(8pi²), with continuous value at y_j=0. The code evaluates this correction stably. Its leading tail is 35Lambda^6*xddot/(256pi²K^6).

The regulator limit is
A_infinity=(x-r)²/(32pi²)+xddot/(96pi²).
This is a matched, scheme-dependent trace relation. It is not a universal cosmological prediction. Fixed K/Lambda alone does not establish either limit; fixed-Lambda momentum comparisons are included.

## Exact incoming state

For z=(1+tanh(t))/2, a+b=1, ab=4 and c=1-iomega_infinity,j,
f_in=e^(-iomega_infinity,j*t) 2F1(a,b;c;z)/sqrt(2omega_infinity,j).
Use the convergent small-z series and its derivative at the numerical starting time. Time reflection and complex conjugation give the exact out state. Its occupation is cosh²(pi sqrt(4-1/4))/sinh²(piomega_infinity,j). Omitting an overall constant mode phase changes no observable.

The intended state is the asymptotic scattering vacuum; choosing a fresh plane wave at a finite time defines a different state and can introduce an unwanted ultraviolet tail. Series truncation and integration are numerical approximations, checked separately. At k=0,x=0 the instantaneous particle basis is singular although the exact field mode is regular; quadrature avoids that measure-zero endpoint.

## Limits of the result

This is the absolute expectation value for a specified Gaussian field and finite matching convention on an external smooth background. It includes its nonadiabatic evolution, but excludes scalar loops, matter interactions, gravity, shell junctions and backreaction. The external scalar supplies the energy. Production does not establish thermal radiation.

Curvature-squared terms have zero first metric variation on this exactly flat background; that does not remove their matching requirements in an FRW/shell extension. The numerical results will provide finite regulator/cutoff evidence, not a mathematical regulator-limit proof.
