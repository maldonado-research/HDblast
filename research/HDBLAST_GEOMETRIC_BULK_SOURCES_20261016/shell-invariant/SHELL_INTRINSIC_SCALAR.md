# An intrinsic material-shell scalar at quadratic order

This is a diagnostic of the registered conditional shell model. The derivation
uses familiar perturbative gauge covariance and does not establish an exact
nonlinear family, a physical force, or a cosmological origin mechanism.

## Definition and embedding convention

Let X_e embed the material shell into the perturbed spacetime, gamma_e=X_e^*g_e
be its Lorentzian induced metric, and psi_e=X_e^*phi_e be the scalar pulled back
to that same shell. Define the intrinsic scalar

    I_e = gamma_e^{mu nu} (d_mu psi_e)(d_nu psi_e).

Assume the background pullback psi_0 is constant on the shell, as in the
retained static-radial background. Use second parameter derivatives:
psi_e=psi_0+e eta_mat+e^2 omega_mat/2+o(e^2). Then d psi_0=0 gives

    I_0=0,  I'_0=0,
    I''_0=2 gamma_0^{mu nu}(d_mu eta_mat)(d_nu eta_mat).       (1)

The first and second variations of the induced metric, and omega_mat, first
enter higher powers of e. Equation (1) is exact at quadratic order; it does
not require solving for the second-order metric or scalar.

The field and embedding must transform jointly. In a normal displacement
description eta_mat=chi_b+p_b zeta, a plus-pullback gauge change has
chi_b -> chi_b+p_b f and zeta -> zeta-f, so eta_mat is unchanged. The registered
first solution is represented with fixed material shell, zeta=0, and its gauge
vector obeys f|_b=0. Thus eta_mat=v(t,0)Y in that representative. This
identification is not asserted for an additional unspecified physical shell
displacement. Tangential gauge relabellings cannot change eta_mat because
psi_0 is constant. More generally the second-order gauge law for I has only
Lie derivatives of I_0 and I'_0, both zero. Therefore I''_0 is invariant under
the joint perturbative field/embedding gauge transformations and shell
relabellings covered by these premises.

## Explicit quadrupole and exact angular consequences

For gamma_0=rho_b^2[-dt^2+c(t)^2 kappa], c=cosh(t), Y=n1^2-n2^2 on the unit
three-sphere and S=n1^2+n2^2, one has |D Y|^2=4(S-Y^2). Hence

    I''_0 = (2/rho_b^2)[-v_t^2 Y^2 + (4v^2/c^2)(S-Y^2)].     (2)

The coefficient of e^2 in I_e is one half of (2). These coefficients use
the model's units; no SI force or finite-amplitude magnitude is inferred.

Probability-normalized sphere moments are <Y^2>=1/6 and <|D Y|^2>=4/3.
Consequently, on the chosen background time slices,

    <I''_0> = [-v_t^2+8v^2/c^2]/(3 rho_b^2).                 (3)

This average is a diagnostic in the fixed background foliation. No new
physical clock or observational sampling prescription is introduced.

At any time with v_t nonzero, (2) is strictly negative at n=(1,0,0,0).
At any time with v nonzero, it is strictly positive at
n=(1/sqrt(2),1/sqrt(2),0,0). Thus when both are nonzero the invariant
quadratic coefficient has both signs on the same angular section. This is
a statement about the Lorentzian squared gradient at perturbative order,
not attractive versus repulsive force. At S=0 it vanishes identically.
For S>0 and v*v_t nonzero its remaining zero set satisfies

    Y^2/S = 4v^2/(c^2 v_t^2+4v^2).                          (4)

The right side is strictly between zero and one. Formula (4) describes the
coefficient's angular nodal set without evaluating a square root. A strictly
negative mean in (3) does not imply negativity at every angle.

## Two simple mean zeros in each of two certified time brackets

Write M(t)=<I''_0> and assume v does not vanish on a closed time bracket.
Define H=c v_t/v and N=v v_tt-v_t^2+theta v v_t, where theta=tanh(t).
Direct differentiation and factorization give

    H_t = c N/v^2,
    M = [v^2/(3 rho_b^2 c^2)] (8-H^2).                     (5)

If N<0 throughout the bracket, H is strictly decreasing. If its left
endpoint is greater than sqrt(8) and its right endpoint is less than
-sqrt(8), the intermediate value theorem and injectivity give exactly one
root H=+sqrt(8) and exactly one root H=-sqrt(8). They are the only M zeros
in that bracket. Both are simple since H_t<0 and H is nonzero there.
M is positive strictly between these roots and negative on the two
remaining portions of the bracket.

The rational actual-jet certificate verifies these conditions on both
scaled-time brackets [5/4,11/8] and [13/8,7/4] at r=0, with t=10^-6 x.
The former has v<0 throughout, the latter v>0. Its endpoint comparisons
use the signs of H and exact squared rational inequalities against 8,
not a floating square root. Its N enclosure covers each entire closed
bracket, not just a sample. Therefore the mean quadratic coefficient has
four distinct simple zeros in the union of these two brackets, exactly
two in each. No count is asserted outside those brackets, nor is this
a statement about zeros of the pointwise scalar at every angle.

These mean zeros are separate from the previously studied coordinate shell
source coefficients and from initial second-order scalar-profile zeros.

## Numerical scope

Any certified numerical assertion about (2)-(4) must use enclosures of the
actual shell v and ordinary t derivative v_t, with frame/time conversion
retained. The predecessor C3 error certificate and the new actual jet oracle
can supply these; the free comparison alone is insufficient. Shell f and
tangential time derivatives vanish by the physical boundary rows, not by
silently discarding an uncertainty. The accompanying numerical report must
authenticate those inputs and distinguish complete time intervals from points.

Even a certified nonzero I''_0 is a result for a formal perturbative coefficient.
Interpreting its sign for a finite perturbation requires an actual nonlinear
family and control of its remainder. Neither is supplied by this diagnostic.
