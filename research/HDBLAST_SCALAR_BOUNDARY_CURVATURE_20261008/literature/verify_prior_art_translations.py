"""Exact convention translations, not a novelty test or branch-existence proof.

Equations transcribed from the fetched primary PDFs documented in the adjacent
receipts. No original HDBLAST numerical producer is imported or executed.
"""
import hashlib
import json
from pathlib import Path
import sympy as s

U, V, sigma, sigma_phi, lam, lam_phi, W, W_phi, W_D_phi = s.symbols(
    "U V sigma sigma_phi lam lam_phi W W_phi W_D_phi", real=True
)
A_y, phi_y, phi_D_y, H2 = s.symbols("A_y phi_y phi_D_y H2", real=True)
p, M, m, k, c, g, q = s.symbols("p M m k c g q", real=True)
mass_D2, mu2, beta, robin_D = s.symbols("mass_D2 mu2 beta robin_D", real=True)
z, eta = s.symbols("z eta", real=True)
checks = {}

def equal(name, left, right):
    residual = s.simplify(s.cancel(s.together(left - right)))
    checks[name] = {"pass": residual == 0, "residual": str(residual)}
    if residual != 0:
        raise RuntimeError(name + ": " + str(residual))

# DFGK Eqs. (6), (8), (18): phi_H=sqrt(2)*phi_D,
# U_H=2*V_D, sigma_H=2*lambda_D, W_H(phi_H)=W_D(phi_H/sqrt(2)).
equal("DFGK_bulk_constraint", (-V / 3 + phi_D_y**2 / 6).subs(
    {V: U / 2, phi_D_y: phi_y / s.sqrt(2)}), -U / 6 + phi_y**2 / 12)
equal("DFGK_superpotential", 2 * (W_D_phi**2 / 8 - W**2 / 3).subs(
    W_D_phi, s.sqrt(2) * W_phi), W_phi**2 / 2 - 2 * W**2 / 3)
equal("DFGK_shell_curvature", (lam**2 / 9 - lam_phi**2 / 24 + V / 3).subs(
    {lam: sigma / 2, lam_phi: sigma_phi / s.sqrt(2), V: U / 2}),
    sigma**2 / 36 - sigma_phi**2 / 48 + U / 6)

# Ghosh et al. (1807.09794) Eq. (2.19), d=4:
# varphi_G=sqrt(2)*phi_H, V_G=2*U_H, R_B=12*H_H^2.
equal("GKNW_bulk_constraint_d4", 12 * A_y**2 - (s.sqrt(2) * phi_y)**2 / 2
    + 2 * U - 12 * H2, 12 * (A_y**2 - phi_y**2 / 12 + U / 6 - H2))

# George Eq. (43), x=-y, lambda=c, beta=-g,
# m_G^2=M^2+k^2, mu_G^2=m^2+k^2, omega^2=M^2+k^2+p^2.
R_D = (p * (p**2 - mu2 + mass_D2) - s.I * (
    robin_D * (p**2 - mu2 + mass_D2) + beta**2)) / (
    p * (p**2 - mu2 + mass_D2) + s.I * (
    robin_D * (p**2 - mu2 + mass_D2) + beta**2))
mapping = {robin_D: c, beta: -g, mass_D2: M**2 + k**2, mu2: m**2 + k**2}
D = m**2 - M**2 - p**2
R_H = ((s.I*p+c)*D-g**2) / ((s.I*p-c)*D+g**2)
equal("George_reflection_matrix", R_D.subs(mapping), R_H)
equal("George_bound_cubic", (q**3 + robin_D*q**2 + (mu2-mass_D2)*q
    + robin_D*(mu2-mass_D2)-beta**2).subs(mapping),
    (q+c)*(q**2+m**2-M**2)-g**2)
equal("George_stability_margin_k0", (mu2*(robin_D+s.sqrt(mass_D2))-beta**2).subs(
    mapping).subs(k, 0).subs(s.sqrt(M**2), M), m**2*(c+M)-g**2)

# Mintchev-Pilo Eq. (34), for z>M^2 and eta=c>0 (no bound atom).
equal("Mintchev_Pilo_continuum_density", (s.sqrt(z-M**2) / (
    s.pi*(z+eta**2-M**2))).subs(eta,c), s.sqrt(z-M**2)/(s.pi*(c**2+z-M**2)))

source = Path(__file__)
print(json.dumps({"status": "PASS_EXACT_PRIMARY_SOURCE_TRANSLATIONS",
    "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
    "sympy": s.__version__, "count": len(checks), "checks": checks,
    "limits": "Algebraic convention checks; no branch proof, priority finding, or empirical cosmology inference."}, indent=2))
