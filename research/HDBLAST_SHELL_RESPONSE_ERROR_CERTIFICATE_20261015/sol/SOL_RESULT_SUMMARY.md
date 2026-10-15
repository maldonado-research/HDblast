# Sol result: conditional physical derivative-trace certificate

The complete derivation is in `RESIDUAL_TO_PHYSICAL_RESPONSE.md`.

For the inherited coupled master, put
\(\gamma=C_b/2\), \(R=w_z-\gamma w_b\). The physical quadratic coefficient is
exactly

\[
\mathcal C_T=\frac{3\pi^2}{\rho_b^2p_b^2}\int_0^T R(t)^2dt.
\]

The hyperbolic-cosine weight cancels. With error orientation \(e=w-u\),
the signed residual equations are
\(e_{tt}-e_{zz}+We=f\),
\(-e_z+\mu e_{tt}+\beta e=g\).
The boundary characteristic system gives
\(\mu v_e'+v_e+\beta q_e=j_+^e+g\).

Let \(a=1/\mu_-\), \(|W|\le M\), \(|\beta|\le B\),
\(|f|\le F\), \(|g|\le G\). With the initial norms specified in the full note,
the nonnegative system

\[
E'=J_0+I+V,\quad I'=ME+F,\quad
V'=aBE+aI+a(b_++G),\quad (E,I,V)(0)=(A_0,0,V_0)
\]

bounds the master error, accumulated incoming characteristic error and boundary
velocity error. In particular

\[
|R-\widehat R|\le b_++I+V+|\gamma|E+
\delta_\gamma|u_b|.
\]

No potential/residual time derivatives enter this estimate. It uses a uniform
initial displacement derivative error, initial velocity error, and continuous
residual bounds. The approximate boundary acceleration enters its residual; a
coherent coefficient normalization avoids multiplying acceleration by an
independent coefficient uncertainty.

If the resulting observation norm error is \(\varepsilon_R\), a rigorous
approximate norm interval is \([N_-,N_+]\), and the physical prefactor is
\([K_-,K_+]\), then

\[
\mathcal C_T\in
[K_-\max(0,N_--\varepsilon_R)^2,\ K_+(N_++\varepsilon_R)^2].
\]

The certificate proves positivity when \(N_->\varepsilon_R\). A vanishing lower
endpoint is inconclusive.

The script `rational_majorant.py` supplies exact-rational matrix-exponential
majorant bounds for rational constant input enclosures. Its saved controls
exercise zero error, pure boundary forcing, simultaneous bulk/boundary forcing,
a nonnilpotent positive matrix, and rejection of omitted boundary residuals and
floating inputs. These are operator controls, not a physical waveform.

The independent `algebra/` audit passes 35 exact identities and rejects 9 wrong
formulas, including incorrect trace/weight choices, boundary acceleration signs,
and omitted coefficient residuals. It also verifies a high-frequency family
whose field supremum tends to zero while its derivative observation does not.

No numerical physical background, waveform, physical coefficient interval, or
effective nonlinear amplitude radius has been produced in this branch. The
retained exact-data, full coupled reconstruction, physical gauge/clock mapping,
and finite-time exact-family hypotheses remain conditional inputs. No previous
source or publication artifact was edited.
