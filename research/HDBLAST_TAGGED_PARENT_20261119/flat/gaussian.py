"""Leading vacuum-subtracted Yukawa rates of a zero-momentum scalar in 3+1 dimensions."""
import json
import math

ROOT_PI = math.sqrt(math.pi)


def simpson(f, end, n):
    h = end/n
    return h/3*math.fsum((f(0), f(end), 4*math.fsum(f(j*h) for j in range(1,n,2)),
                           2*math.fsum(f(j*h) for j in range(2,n,2))))


def branch(u, power, creating, n):
    # z = omega*T. Pull out exp(-u²) in the very small creation tail.
    def integrand(z):
        exponent = -z*(z + 2*u) if creating else -(z-u)**2
        return (z/u)**power*math.exp(exponent)

    if creating:
        # q=2uz resolves the endpoint layer at large u. The discarded tail
        # after q=64 is bounded by the corresponding e^-q moment.
        return math.exp(-u*u)*simpson(lambda q: integrand(q/(2*u)),64,n)/(2*u*ROOT_PI)
    return simpson(integrand,u+14,n)/ROOT_PI


def moment_check(u, n):
    d, r = (branch(u,2,up,n) for up in (False,True))
    ed, er = (branch(u,3,up,n) for up in (False,True))
    return dict(annihilation=d, stimulated_creation=r, inclusive_pairs=d+r,
                signed_depletion=d-r, daughter_energy=ed+er, external_work=(ed-d)+(er+r))


def calculation(u):
    g = math.exp(-u*u)/(u*ROOT_PI)
    erf = math.erf(u)
    count = 1 + 1/(2*u*u)
    lost = count*erf + g
    energy = (1+3/(2*u*u))*erf + (1+1/(u*u))*g
    work = erf/(u*u) + g/(u*u)
    # erfc avoids inferring the creation probability from a rounded difference
    # of nearly equal O(1) numbers. The direct positive integral cross-checks
    # the remaining polynomial cancellation, strongest here at u=8.
    r = (count*math.erfc(u)-g)/2
    closed = dict(annihilation=count-r, stimulated_creation=r, inclusive_pairs=count,
                  signed_depletion=lost, daughter_energy=energy, external_work=work)
    grids = {str(n): moment_check(u,n) for n in (4096,8192)}
    for key, exact in closed.items():
        if exact <= 0 or not math.isclose(grids["8192"][key],exact,rel_tol=2e-9):
            raise ValueError(f"Failed direct quadrature u={u} {key}")
    if not math.isclose(energy-lost,work,rel_tol=3e-14):
        raise ValueError("Signed energy identity failed")
    bound = math.exp(-u*u)/(4*ROOT_PI*u**5)
    if r > bound:
        raise ValueError("Creation tail upper bound failed")
    return dict(mT=u, analytic=closed, quadrature=grids, stimulated_creation_upper_bound=bound)


if __name__ == "__main__":
    print(json.dumps(dict(order="y^2",number_unit="Gamma0 sqrt(pi) T",
                          energy_unit="m Gamma0 sqrt(pi) T",Gamma0="y^2 m/(8 pi)",
                          data=[calculation(u) for u in (1,2,4,8)]),indent=2,allow_nan=False))
