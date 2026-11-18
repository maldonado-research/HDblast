"""Integrate quiet and perturbed absolute B1 fields with independent fixed-grid RK4."""
import argparse
import hashlib
import json
import math
from pathlib import Path

PIN = "bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f"


def potential(x):
    superpotential = 1 + x*(x*x/3-1)
    derivative = (x-1)*(x+1)
    return (derivative**2/2-2*superpotential**2/3,
            2*x*derivative-4*superpotential*derivative/3)


def field_pulse(z, parameters, coefficients):
    pulse = derivative = 0.0
    width = parameters["half_width"]
    for center, amount in zip(parameters["centers"], coefficients):
        x = (z-center)/width
        if abs(x) < 1:
            d = 1-x*x
            bump = math.exp(1-1/d)
            pulse += amount*bump
            derivative -= 2*amount*x*bump/(width*d*d)
    return pulse, derivative


def integrate(summary, parameters, coefficients, nodes):
    rho, phi = summary["rho_b"], summary["phi_b_static"]
    p = summary["params"]
    c, d, det = p["c"], p["d"], p["delta"]
    b = math.log(rho)
    q = rho*(2*(1-phi+phi**3/3)+det*(1+c*phi+d*phi**2/2))/6
    r = -rho*(2*(phi**2-1)+det*(c+d*phi))/2
    state = (b,q,phi,r,b,q,1.0)
    spacing = -0.20/nodes

    def derivative(z, values):
        bg, slope, x, grad, geom, q, velocity = values
        U, Up = potential(x)
        bg_kin = math.exp(2*bg)
        h, j = field_pulse(z, parameters, coefficients)
        total_grad, speed = grad+j, -j
        acc = (velocity**2 - q*q - (speed*speed + total_grad**2)/6
               - math.exp(2*geom)*potential(x+h)[0]/3)
        return (slope, 1-slope*slope-grad*grad/6-bg_kin*U/3,
                grad, bg_kin*Up-3*slope*grad, q, acc, -speed*total_grad/3)

    def moved(values, step, func):
        return tuple(x + step*d for x, d in zip(values,func))

    grid = []
    for index in range(nodes+1):
        z = index*spacing
        h,j = field_pulse(z,parameters,coefficients)
        grid.append((state[4], state[2]+h, state[6], -j))
        if index == nodes:
            break
        k1 = derivative(z,state)
        k2 = derivative(z+spacing/2, moved(state,spacing/2,k1))
        k3 = derivative(z+spacing/2, moved(state,spacing/2,k2))
        k4 = derivative(z+spacing, moved(state,spacing,k3))
        state = tuple(v + spacing*(a+2*b+2*c+d)/6
                      for v,a,b,c,d in zip(state,k1,k2,k3,k4))
    h_max = m_max = 0.0
    for index in range(2, nodes-1):
        B, x, u, ptime = grid[index]
        def first(col):
            return (grid[index-2][col]-8*grid[index-1][col]
                    +8*grid[index+1][col]-grid[index+2][col])/(12*spacing)
        Bz, r, uz = (first(col) for col in (0,1,2))
        Bzz = (-grid[index+2][0]+16*grid[index+1][0]-30*B
               +16*grid[index-1][0]-grid[index-2][0])/(12*spacing**2)
        terms = [-2*math.exp(2*B)*potential(x)[0], 6*u*u, -12*Bz*Bz,
                 6*Bz*Bz, -6*Bzz, -ptime*ptime, -r*r]
        momenta = [-3*uz, -3*u*Bz, 3*u*Bz, -ptime*r]
        h_max = max(h_max,abs(math.fsum(terms))/sum(map(abs,terms)))
        m_max = max(m_max,abs(math.fsum(momenta))/sum(map(abs,momenta)))
    return {"intervals":nodes,"outer_absolute_field_difference":[state[4]-state[0],state[5]-state[1],state[6]-1],
            "max_scaled_H_original_finite_difference":h_max,
            "max_scaled_M_original_finite_difference":m_max}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary",type=Path)
    args = parser.parse_args()
    raw = args.summary.read_bytes()
    if hashlib.sha256(raw).hexdigest() != PIN:
        raise ValueError("Unexpected B1 source summary")
    summary = json.loads(raw)
    profiles = json.loads(Path(__file__).with_name("PROFILES.json").read_text())
    result = {}
    for name in ("quiet_connected","second_completion"):
        result[name] = [integrate(summary,profiles,profiles[name],n) for n in [256,512,1024,2048]]
        if max(abs(v) for v in result[name][-1]["outer_absolute_field_difference"]) > 5e-12:
            raise ValueError("Independent absolute-field verification failed")
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
