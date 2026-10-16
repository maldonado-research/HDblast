# Continuous approximation errors for the actual retained master and gauge

Unpublished numerical checkpoint, 16 October 2026. The retained smooth branch,
master, gauge and their previously proved integral/reflection identities are
conditional authenticated premises. The programs read pinned JSON copies and
never execute earlier research modules. They add actual approximation-error
certificates; they do not evaluate the second response.

## Master comparison and derivative errors

Use ell=10^-6, x=t/ell, r=z/ell and D={0<=x<=3,r>=0,x+r<=8}.
Q(r)=exp(4-1/((r-1)(2-r))) on (1,2), zero otherwise. Write
Qo(s)=sign(s)Q(abs(s)), Qe(s)=Q(abs(s)), and
Wbar=(Qo(r+x)+Qo(r-x))/2. All these extensions are smooth because of the
fixed zero collar and flat packet endpoints. The actual master has initial
velocity ell V(r), V=exp(-ell r)Q(r); Wbar has zero velocity.

The retained exact decomposition is W=Wbar+Wvelocity+Ly-Kf, where
f=ell^2 mathcalW W and y=W(x,0) is the actual dynamic shell trace.
For every mixed derivative of total order k=0,1,2,3, respectively,

* Wvelocity is bounded by 3ell, ell, ell(Q1+ell),
  ell(Q2+2ell Q1+ell^2).
* Ly is bounded by y,y',y'',y''' from the retained dynamic trace equation.
* Kf is bounded by (9/2)||f||, retained first_error, maximum second_error,
  maximum third_error.

The last two rows come from differentiated actual integral formulas. In
particular second_error includes the reflected shell_f term, and third_error
includes shell_fx, with the characteristic joins already justified in the
retained proof. The retained free_second/free_third also include the bump
and lift; we do not reuse those totals. Thus no trace term is counted twice.
The bound Ek is the sum of the three displayed rows. It is an error bound
for W-Wbar and its actual derivatives, not a differentiation of a C0 estimate.

## Forced gauge comparison and complete true residual

Set c=6ell/130, Pe(s)=s Qe(s), which is odd, and
P=c/4[Pe(r+x)+Pe(r-x)-r(Qo(r+x)+Qo(r-x))].
Then P has zero initial displacement and velocity, P(x,0)=0, and
(Pxx-Prr)=c Wbar_r. The Pe terms solve the homogeneous wave equation;
differentiating the factor r in the remaining terms supplies the source.
Let V=(P,0,0) in the registered raw gauge channels (f,b,d). It satisfies
all three homogeneous physical gauge boundary rows exactly. The bound
for every mixed derivative of total order k is

  pk=c[5Qk+k Q(k-1)], k>=1; p0=5c.

Indeed |Pe^(k)|<=2Qk+kQ(k-1), |r|<=8 and differentiating r contributes
at most jQ(k-1) for j radial derivatives. The sum of the two characteristic
terms gives c/2[(2+8)Qk+(k+j)Q(k-1)] and j<=k gives the stated cap.

For the actual U equation Uxx-Urr=A Ux+B Ur+D U+(S,0,0),
S=Ps W_r+P0 W, E=U-V satisfies the SAME homogeneous mixed boundary
problem with zero initial value/velocity and true forcing residual

  R=A Vx+B Vr+D V+[(Ps-c)Wbar_r+Ps(W-Wbar)_r+P0 W]e_f.

Every mixed derivative of R through total order two is enclosed by the
ordinary product rule, with the actual coefficient intervals and Ek. Ps-c
is subtracted before taking absolute values. This is a complete true PDE
residual. Its source part alone is also reported; it is an approximation
error for S-c Wbar_r. The temporal/spatial conversions and the true master
are inherited in Ps,P0. This does not yet evaluate the common-gauge
Einstein-scalar quadratic sources.

## Independent derivative-error budget

For theta=3/2 use Ni(E)=max(||Ei||/theta,||Ei_x||,||Ei_r||).
M=3(A*+B*+theta D*) has nonnegative exact inverse Rm=(I-M)^-1.
All following vectors are nonnegative componentwise majorants. Let r0,rx,
rr,rxx be the corresponding forcing-residual bounds, and
Lx=Ax*+Bx*+theta Dx*, Lxx=Axx*+Bxx*+theta Dxx*.

  n0=Rm 3r0,
  n1=Rm [2r0+3(Lx n0+rx)],
  h1=A* r0+rx,
  gw=max(2r0/3+2h1, rr+h1),
  n2=Rm [gw+3(2Lx n1+Lxx n0+rxx)].

These are bounds for N(E),N(Ex),N(Exx). The stabilized actual initial data
are (Ex,Exx)=(0,R0) and (Exx,Exxx)=(R0,A(0)R0+Rx0).
In particular gw retains the radial derivative of R0 and the nonzero
A(0)R0 term. These data cannot be obtained by setting higher iterated
initial jets to zero. The retained differentiated mixed-wave construction
proves the required convergence; its linear estimates apply to E because
V is a smooth explicit function with exact matching boundary/initial rows.

Exx is bounded by min(n1,theta n2); Exxx and Exxr by n2. The remaining
radial derivatives use the true PDE, including the complete Rr in Errr.
The JSON lists each of the ten actual mixed derivative errors separately.
This bound includes the shell, reflected transition and packet joins. The
comparison fields themselves are explicitly defined continuous functions;
the error certificate does not differentiate independent interval selections.

## Limits and next step

These are fully defined comparison functions and rigorous errors relative to
the actual retained first solutions. No numerical fixed-point iterations have
been performed, and the comparisons are deliberately inexpensive. The master
third error is dominated by the actual dynamic trace jet, while the gauge
third error is dominated by the source-coefficient freezing discrepancy.
Refinement should evaluate the actual differentiated integral iterations,
using the correlated background approximation and its errors as inputs.
No second waveform, common-gauge quadratic-source values, finite-amplitude
remainder, invariant force, or cosmological inference follows from these files.
