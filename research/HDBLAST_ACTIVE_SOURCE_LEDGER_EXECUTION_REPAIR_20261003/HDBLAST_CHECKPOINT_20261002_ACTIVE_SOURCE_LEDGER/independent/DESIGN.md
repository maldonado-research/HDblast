# Independent active-source route (prospective)

No physical source, phase, mode, stress, or ledger evaluation has occurred during this preparation. The fixed proposed interval is [-4.5,-3.5], with retained observation anchors 1 and 3. The two sources, coarse/fine inherited momentum rules, and K=64,128,256 give twelve cases. Older FAIL receipts remain unchanged.

This route starts from the actual raw u,w at anchor 1, without projecting their normalization defect. It evaluates the exact analytic bump and its directional derivatives at fixed quadrature nodes. The ODE is u'=w, w'=2ik w-epsilon*g. Exact oscillatory propagation is combined with independent fixed-order Gaussian source quadrature. No polynomial source surrogate is used. Source nodes and phases use the pinned native long-double arithmetic; final weighted ledger/flow endpoint terms are reduced in separate 80- and 100-decimal-digit contexts. Those contexts check reduction arithmetic only: they are not full high-precision source/phase calculations and do not certify the inherited initial state.

The independently expanded bare stress terms, with U=u/epsilon, W=w/epsilon, are

```
r_k = [(2k^2+3L^2)Re(U)-L Re(W)-k Im(W)+L h'+2L^2 h]/(2k)
p_k = [(2k^2/3-L^2)Re(U)-L Re(W)-k Im(W)+L h'-2L^2 h]/(2k).
```

The complete metric subtraction is differentiated directly using pair-valued directional algebra applied to the published W0/W2/W4 minimal-scalar inventory. Exact symbolic regrouping converts that independently expanded inventory to weighted inverse-frequency moments. These discrete moments and the rationalized baseline are computed once for each grid/quadrature geometry, then reused across both sources. Both stresses remain separately defined, and no conservation formula defines a contact. The fixed physical mass, operator variation, scale/prefactor contacts, and the rationalized finite-band baseline are included. Pressure is not defined by conservation. For every dense time node, this route constructs

```
R_k=r_k-delta(subtraction_R)_k-4h R0_k,
P_k=p_k-delta(subtraction_P)_k-4h P0_k,
F_k=L*(R_k-3P_k)-3h'*(R0_k+P0_k).
```

The measured final stored R is never used to define the canonical integral. Instead, I is computed by a separately weighted direct quadrature of these F values at dense local nodes. Dense modes use explicit local Duhamel kernels, not a finite difference or a conservation-based stress. The fixed three numerical controls are GL16/GL16, GL16/GL24, and GL24/GL24 for local source quadrature/dense F quadrature. GL16/GL16 and GL16/GL24 share their endpoint trajectory, isolating dense ledger quadrature; GL16/GL24 and GL24/GL24 hold ledger quadrature fixed and compare forcing quadrature. These fixed control differences are empirical and supply no rigorous error bound. Quadrature constants come from binary64 NumPy leggauss nodes/weights cast to native long double, which is an explicit algorithm distinction from the primary route. The raw inherited momentum nodes/weights and scalar pi/epsilon preserve their native binary80 values.

At anchor 3, raw stored modes are compared with the freshly propagated modes. The resulting signed density projection is

```
e_flow,k=[(2k^2+3L_b^2)Re(delta_u)/k-Im(delta_w)-L_b Re(delta_w)/k]/(2epsilon),
```

with a separate coefficient triangle bound. Identical endpoint source contacts cancel in this difference. This bounds a projection of the measured flow defect only; it does not bound the inherited anchor-1 state or the exact-source quadrature error. Freshly evaluated source contacts and stored contacts are compared as separate profiles.

The inherited global Simpson prefixes and endpoint stresses remain unchanged, with the original long-double subtraction followed by exact-ratio conversion. The signed decomposition is D_S=DeltaR-S, D_cont=DeltaR-I, E_Q=I-S, D_S=D_cont+E_Q. D_cont includes trajectory, contact-reading, and direct-ledger discretization differences. A small empirical residual does not prove a total error certificate. The old endpoint/refinement gates 2e-6 and 1e-6 are preserved; the old experiment remains FAIL irrespective of this new diagnosis.

Momentum nodes and weights are exactly the inherited values; this route establishes no continuum-momentum error bound. Momentum blocks cap memory. The complete route, provenance through output serialization, is limited to 900 seconds and 262144 KiB. Only fabricated-source, fabricated-grid, and exact-symbolic controls may run before the public registration GO.
