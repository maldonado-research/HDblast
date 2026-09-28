#!/usr/bin/env python3
"""HDBLAST Chat 9 - literature cross-check script (numpy + stdlib only).

Purpose: apply two PRIOR-ART methods to the registered HDBLAST shell and compare them.

 (A) Moduli-space approximation (MSA) of Brax, van de Bruck, Davis, Rhodes, hep-th/0209158, eqs. (13), (25), (42),
     specialised here to ONE Z2 brane in a BPS bulk with general superpotential W (kappa_5 = 1):
         L_Jordan = f(z) R + a^2 W (dz)^2 - a^4 V(z),   f = int_z^inf a^2,  a = e^A, A' = -W/3, phi' = W_phi,
         V = t (1 + c phi)     (detuning of the tension sigma_t = 2W + t(1 + c phi))
     Einstein frame (M_pl = 1):  L = R/2 - (1/2) G_zz (dz)^2 - V_E,
         G_zz = (3/2) a^4/f^2 - a^2 W/f,     V_E = a^4 V/(4 f^2).
     At a stationary point of V_E:  mu^2 = m^2/H^2 = 3 (ln V_E)_zz / G_zz.
     For the registered wall (phi = -tanh z, brane at phi_b = 0, a = 1, f = I_plus):
         stationarity  <=>  c = 2/I_plus - 4/3 = c_star      (this IS the registered c)
         (ln V_E)_zz = -2/I^2 + 4/I - 28/9,   G_zz = 3/(2 I^2) - 1/I.

 (B) Frolov-Kofman hep-th/0309002 scalar perturbation equations (8a),(8b) with their boundary condition (14),
     transcribed to the registered chart ds^2 = dy^2 + rho^2 dS_4(unit), harmonic Box_gamma Q = mu^2 Q:
         Y = rho^2 Phi,  Z = dphi/phi',     Y' = (2/3) rho^2 phi'^2 Z,     Z' = [1 - (3/2)(mu^2+4)/(rho^2 phi'^2)] Y/rho^2
     regular cone:  Y ~ y^p,  p = (5 + sqrt(9 - 4 mu^2))/2   (equivalently dphi ~ y^gamma, mu^2 = -gamma(gamma+3),
                    the Himemoto-Sasaki gr-qc/0010035 regularity condition)
     shell (bulk on y < y_b, tension sigma_t(phi)):   dphi' - phi' Phi = -(1/2) sigma_t'' dphi
         <=>  Z phi'^2 [U_phi/phi' - 4 rho'/rho + sigma_t''/2] = (3/2)(mu^2+4) Y/rho^4.

 Output: mu^2 from (A) (t -> 0 limit) and from (B) at several t (c = c_star fixed), to see whether (B) -> (A) as t -> 0.
 Floating-point diagnostic, NOT a certificate.
"""
import math, sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "background"))
import hdblast_background as bg

def W(p):  return 1 - p + p**3/3
def W1(p): return p*p - 1
def W2(p): return 2*p
def U(p):  return 0.5*W1(p)**2 - (2/3)*W(p)**2
def U1(p): return W1(p)*W2(p) - (4/3)*W(p)*W1(p)
def U2(p): return W2(p)**2 + W1(p)*2 - (4/3)*(W1(p)**2 + W(p)*W2(p))

def msa(I):
    lnV_zz = -2/I**2 + 4/I - 28/9
    G_zz = 1.5/I**2 - 1/I
    return lnV_zz, G_zz, 3*lnV_zz/G_zz

def shoot(mu2, phi_h, v_b, t, c, dv=2e-4, v0=1e-2):
    """integrate background + FK perturbation from the cone to the shell; return boundary mismatch (normalised)."""
    rho, phi, s = [float(x) for x in bg.series(phi_h, v0)]
    p = (5 + math.sqrt(9 - 4*mu2))/2
    Y = 1.0                                  # overall scale is arbitrary (linear problem)
    Z = 1.5*p*Y/v0/(rho*rho*s*s)             # from Y' = (2/3) rho^2 s^2 Z with Y = (y/v0)^p
    n = max(1, int(math.ceil((v_b - v0)/dv))); h = (v_b - v0)/n
    lam = mu2 + 4
    def f(rho, phi, s, Y, Z):
        rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
        return (rp, s, U1(phi) - 4*(rp/rho)*s,
                (2/3)*rho*rho*s*s*Z,
                Y/(rho*rho) - 1.5*lam*Y/(rho**4*s*s))
    for _ in range(n):
        y0 = (rho, phi, s, Y, Z)
        k1 = f(*y0)
        k2 = f(*[a + h/2*b for a, b in zip(y0, k1)])
        k3 = f(*[a + h/2*b for a, b in zip(y0, k2)])
        k4 = f(*[a + h*b for a, b in zip(y0, k3)])
        rho, phi, s, Y, Z = [a + h/6*(b1 + 2*b2 + 2*b3 + b4) for a, b1, b2, b3, b4 in zip(y0, k1, k2, k3, k4)]
        sc = abs(Y) + abs(Z)
        if sc > 1e100: Y /= sc; Z /= sc
    rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
    B = U1(phi)/s - 4*rp/rho + 0.5*(2*W2(phi))          # sigma_t'' = 2 W''(phi)   (the detuning t(1+c phi) is linear)
    lhs = Z*s*s*B; rhs = 1.5*lam*Y/rho**4
    return (lhs - rhs)/(abs(lhs) + abs(rhs)), B, rho

def eigen(phi_h, v_b, t, c, lo, hi, dv):
    flo = shoot(lo, phi_h, v_b, t, c, dv)[0]; fhi = shoot(hi, phi_h, v_b, t, c, dv)[0]
    if flo*fhi > 0: return None
    for _ in range(60):
        mid = 0.5*(lo + hi); fm = shoot(mid, phi_h, v_b, t, c, dv)[0]
        if flo*fm <= 0: hi, fhi = mid, fm
        else: lo, flo = mid, fm
        if hi - lo < 1e-10: break
    return 0.5*(lo + hi)

if __name__ == "__main__":
    I = bg.I_plus(); c = 2/I - 4/3
    lnV_zz, G_zz, mu2_msa = msa(I)
    print("I_plus = %.13f   c_star = %.13f" % (I, c))
    print("(A) MSA:  (ln V_E)_zz = %.8f   G_zz = %.8f   mu^2_MSA = 3 (ln V_E)_zz / G_zz = %.6f" % (lnV_zz, G_zz, mu2_msa))
    print("    growth exponent  sqrt(9/4 - mu^2) - 3/2 = %.5f  (perturbation ~ exp(that * H t))" % (math.sqrt(2.25 - mu2_msa) - 1.5))
    out = dict(I_plus=I, c_star=c, lnV_zz=lnV_zz, G_zz=G_zz, mu2_MSA=mu2_msa, FK=[])
    dv = 4e-4
    guess = None
    for t in [1e-3, 3e-4, 1e-4]:
        if guess is None: g = (math.log10(8.7855e-7), 8.6)
        else: g = guess
        sol = bg.solve_shell(t, c, guess=g, dv=dv)
        # continuation guess for next (smaller) t: alpha ~ t^1.8, y_b shifts by -(l/2) dln t with l = 9/5
        guess = (math.log10(sol['alpha']) + 1.8*math.log10(0.3), sol['v_b'] + 0.9*math.log(1/0.3))
        print("\n t = %g: phi_h+1 = %.6e  y_b = %.8f  rho_b = %.6f  phi_b = %.4e  h/(t/6I) = %.8f  resid = %.1e %.1e" % (
            t, sol['alpha'], sol['v_b'], sol['rho_b'], sol['phi_b'], sol['h']/(t/(6*I)), *sol['residual']))
        # scan for sign changes of the boundary mismatch
        grid = [-12 + 0.5*i for i in range(0, 29)]      # mu^2 in [-12, 2]
        vals = [shoot(m, sol['phi_h'], sol['v_b'], t, c, dv)[0] for m in grid]
        roots = []
        for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
            if fa*fb < 0:
                r = eigen(sol['phi_h'], sol['v_b'], t, c, a, b, dv)
                if r is not None: roots.append(r)
        Bb = shoot(-5.0, sol['phi_h'], sol['v_b'], t, c, dv)[1]
        print("   shell bracket B = U_phi/phi' - 4 rho'/rho + sigma''/2 = %.6e   (B/t = %.5f)" % (Bb, Bb/t))
        print("   (B) FK eigenvalues mu^2 in [-12,2]:", ["%.6f" % r for r in roots])
        out["FK"].append(dict(t=t, roots=roots, B=Bb, y_b=sol['v_b'], rho_b=sol['rho_b'], phi_b=sol['phi_b']))
    json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "prior_art_crosscheck_output.json"), "w"), indent=1)
