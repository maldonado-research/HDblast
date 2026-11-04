# Standard local conservation identities for the retained pure-tension shell

4 November 2026. Supporting derivation. These are standard algebraic and
Codazzi/Ward consequences, not a new physical law or a construction of an
incoming state. The primary comparison is Maeda–Wands equations (2.11), (2.15),
(2.19), (2.26) and Langlois–Rodriguez-Martinez equations (73)–(75). Their normal
conventions must be converted before applying a sign.

## Actual shell stress and a local exchange ledger

In the retained normalization of the bulk action R/2 (kappa5 squared equals 1)
and one-copy outward rows,

```
S_ab=-sigma(phi)*gamma_ab,
Kout_ab=(sigma/6)*gamma_ab,
nout(phi)=-sigma_phi/2,
Ric_AB=partial_A phi*partial_B phi+(2/3)V*g_AB.
```

For every unit future shell observer u, with h_ab=gamma_ab+u_a*u_b,

```
rho_shell=S_ab*u^a*u^b=sigma,
p_shell=(1/3)*h^ab*S_ab=-sigma,
rho_shell+p_shell=0.
```

This is the actual stress from variation of the shell action. It is distinct
from effective scalar/Weyl terms obtained by projecting the bulk Einstein equation.
With the usual retained Codazzi sign,

```
D_a(Kout^a_b)-D_b Kout=Ric(nout,e_b),
-(1/2)*sigma_phi*D_b phi=nout(phi)*D_b phi.
```

Equivalently, for the two identical interiors,

```
D_a S^a_b=2*nout(phi)*D_b phi=-sigma_phi*D_b phi.
```

Writing dot(phi)=u(phi), Theta=D_a u^a, the observer energy equation is

```
dot(rho_shell)+Theta*(rho_shell+p_shell)=sigma_phi*dot(phi).
```

For a compact shell section transported by the stated normal proper-time flow,
`E_sigma=int sigma dV3` obeys

```
dE_sigma/dtau-int sigma*Theta dV3=int sigma_phi*dot(phi) dV3.
```

The right side is energy transferred from the interiors into changing tension;
the second term is the negative-pressure work term. It is not a radiation-
production rate. With no independent shell matter term, rho and p keep the pure-
tension relation even when sigma changes.

## What an actual matter extension changes

Suppose a new local shell action term L_m is supplied, depending on chi and its
tangential derivatives, phi algebraically, and the induced metric, with no
derivatives of phi or extrinsic-curvature terms. Define tau by its metric variation
and S=-sigma*gamma+tau. In this same outward convention the formal rows become

```
Kout_ab=(1/2)*(S_ab-(S/3)*gamma_ab),
2*nout(phi)=-sigma_phi+partial L_m/partial phi.
```

On the matter equations, diffeomorphism invariance then gives

```
D_a tau^a_b=(partial L_m/partial phi)*D_b phi.
```

For an independently specified conformal matter action
`S_m[chi,tilde_gamma]`, define

```
tilde_gamma_ab=A(phi)^2*gamma_ab, a(phi)=partial_phi log A,
J_m=(1/sqrt(-gamma))*delta S_m/delta phi.
```

The same stress convention `tau_ab=-2/sqrt(-gamma)*delta S_m/delta gamma^ab`
and covariant metric variation give

```
delta_phi tilde_gamma_ab=2a*tilde_gamma_ab*delta phi,
sqrt(-tilde_gamma)*tilde_tau_trace=sqrt(-gamma)*tau_trace,
J_m=a*tau_trace, D_a tau^a_b=a*tau_trace*D_b phi.
```

In the algebraic-phi setting above, J_m is simply `partial L_m/partial phi`.
Writing the metric instead as `tilde_gamma=Omega^(-2)*gamma` gives
`a=-Omega_phi/Omega`; the opposite `Omega^2` definition gives the opposite
coefficient. Matter measure and stress definition must be matched before a
printed conformal coefficient is transferred. A strictly traceless radiation
fluid has zero classical exchange from either conformal coupling. This statement
does not include trace anomalies or a specified quantum state. Calling a fluid
radiation does not supply a particle-production mechanism.

For example, L_m=-|Dchi|^2/2-U(phi,chi) would give the shell scalar source
`-U_phi/2` in the normal boundary law and matter exchange `-U_phi*D_b phi`.
A coupling U=g^2*phi^2*chi^2/2 is an additional model choice. Its zero classical
chi solution with zero classical data stays zero. Quantum production requires a
specified state, pump history and mode/renormalization calculation, with subsequent
backreaction and thermalization analysis. None is included in the retained model.

These identities do not determine the incoming bulk state or the free second
metric preparation. In particular a local initial-shell energy cannot distinguish
compact interior preparations sharing the initial collar. An interior/past state
or compatible generating process remains necessary.


## Distinction from the accompanying preparation diagnostics

The companion incoming/preparation ledger uses normalized linear packet
components and an auxiliary estimating energy. Incoming/outgoing combinations
for the principal wave operator, and their positive diagnostic weights, are not
the physical Noether energy or the stress-tensor flux used in the Ward identity
above. A relative canonical constraint generator can vanish on shell while an
observer-defined matter or tension energy has a nonzero value. These quantities
must not be substituted for one another.

Likewise, a forced estimate on a future cone and partial null-constraint
integration do not by themselves construct a retarded past solution. A past
domain, nontrivial first-state or source, constraint-compatible second data and
physical boundary/corner conditions still have to be supplied. No shell matter
sector, particle population or radiation-generation mechanism is implemented by
this supporting calculation.
