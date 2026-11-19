"""Closed-sphere photon production, independently soluble pulse and ultraviolet tail."""
import argparse
from fractions import Fraction
import json
from math import comb, pi
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import solve_ivp


def original_modes(frequencies, kinetic, start, end, max_step):
    """Evolve original A, P, accumulated work with incoming vacuum data."""
    k = np.asarray(frequencies, dtype=float)
    initial_I, _ = kinetic(start)
    field = 1/(initial_I*np.sqrt(2*k))
    momentum = -1j*initial_I*np.sqrt(k/2)
    def rhs(t, state):
        field, momentum, _ = state.reshape(3, -1)
        I, pump = kinetic(t)
        electric = abs(momentum/I)**2
        magnetic = abs(k*I*field)**2
        return np.concatenate((momentum/I**2, -(k*I)**2*field, pump*(magnetic-electric)))
    solved = solve_ivp(rhs,(start,end), np.concatenate((field,momentum,0*k)),
                       rtol=4e-13,atol=5e-14,method="DOP853",max_step=max_step)
    if not solved.success:
        raise ArithmeticError(solved.message)
    field,momentum,work = solved.y[:,-1].reshape(3,-1)
    I,_ = kinetic(end)
    beta = (np.sqrt(k)*I*field-1j*momentum/(I*np.sqrt(k)))/np.sqrt(2)
    number = abs(beta)**2
    return number, float(max(abs(work.real-k*number))), float(max(abs(
        1j*(field.conj()*momentum-momentum.conj()*field)-1)))


def fast_scattering(k, pump, end, tolerance):
    k = np.asarray(k, dtype=float)
    def rhs(t, values):
        alpha,beta = values.reshape(2,-1)
        phase = np.exp(2j*k*t)
        return pump(t)*np.concatenate((phase*beta,phase.conj()*alpha))
    sol = solve_ivp(rhs,(0,end),np.concatenate((1+0j*k,0j*k)),
                    atol=tolerance*0.05,rtol=tolerance,method="DOP853",
                    max_step=min(end/24,3/max(k)))
    if not sol.success:
        raise ArithmeticError(sol.message)
    return sol.y[:,-1].reshape(2,-1)


def exact_sech():
    g, tau = .2, .18
    def couplings(sign):
        def kinetic(eta):
            t=eta/tau
            return np.exp(sign*2*g*np.arctan(np.exp(t))), sign*g/(tau*np.cosh(t))
        return kinetic
    k=np.array([2,3,4,7,12,16],dtype=float)
    def occupation(k):
        z=np.exp(-2*pi*tau*k)
        return np.sinh(pi*g)**2*4*z/(1+z)**2
    exact=occupation(k)
    positive, work, wronskian=original_modes(k,couplings(1),-35*tau,35*tau,tau/12)
    negative, neg_work, neg_wr=original_modes(k,couplings(-1),-35*tau,35*tau,tau/12)
    upper=96
    all_k=np.arange(2,upper+1,dtype=float)
    values=occupation(all_k)
    energy = np.dot(all_k*(all_k**2-1),values)/pi**2
    # Bounding each later (k+1)^3/k^3 by the first such ratio gives a geometric tail.
    next_k = upper+1
    ratio = np.exp(-2*pi*tau)*((next_k+1)/next_k)**3
    tail=(4*np.sinh(pi*g)**2/pi**2*next_k**3*np.exp(-2*pi*tau*next_k)/(1-ratio))
    checks=dict(max_absolute_occupation_error=float(max(abs(positive-exact))),
                max_relative_occupation_error=float(max(abs(positive/exact-1))),
                max_reciprocal_difference=float(max(abs(positive-negative))),
                max_work_residual=max(work,neg_work), max_wronskian_residual=max(wronskian,neg_wr))
    if (checks["max_relative_occupation_error"]>1e-7 or checks["max_work_residual"]>1e-9
            or checks["max_reciprocal_difference"]>1e-10 or checks["max_wronskian_residual"]>1e-9):
        raise ArithmeticError(checks)
    return dict(g=g,tau=tau,log_endpoint_ratio=pi*g,full_sphere_a4_density=float(energy),
                highest_k=upper,analytic_tail_upper=float(tail),
                total_count_on_sphere=float(2*np.dot(all_k**2-1,values)),
                squeeze_upper_per_oscillator=float(np.sinh(pi*g)**2),
                samples=[dict(k=int(j),exact=float(n),canonical_ode=float(o))
                         for j,n,o in zip(k,exact,positive)],checks=checks)


def l1_upper(poly, bisections=4):
    """Exact upper bound to the integral of |polynomial| in rational Bernstein form."""
    deg=poly.degree()
    bern=[sum(poly.coef[j]*Fraction(comb(i,j),comb(deg,j)) for j in range(i+1))
          for i in range(deg+1)]
    pieces=[bern]
    for _ in range(bisections):
        children=[]
        for values in pieces:
            front,back=[],[]
            while len(values)>1:
                front.append(values[0])
                back.append(values[-1])
                values=[(a+b)/2 for a,b in zip(values,values[1:])]
            children.extend([front+values, values+back[::-1]])
        pieces=children
    return sum(sum(map(abs,piece)) for piece in pieces)/((deg+1)*len(pieces))


def exact_derivative(poly):
    return Polynomial([Fraction(i)*q for i,q in enumerate(poly.coef) if i])


def exact_primitive(poly):
    return Polynomial([Fraction(0)]+[q/Fraction(i+1) for i,q in enumerate(poly.coef)])


def compact_pulse():
    """Log I has a beta(5,5) ramp over an interval, extended constantly on both sides."""
    step=Fraction(1,4)
    width=2
    # Polynomials in y=eta/width. Exact field-coefficient arithmetic until output.
    u=Polynomial([Fraction(0)]*4+[630*step*(-1)**i*comb(4,i) for i in range(5)])
    potential=exact_derivative(u)+u*u
    J=exact_primitive(potential)
    residual=(-exact_derivative(exact_derivative(potential))+2*potential**2
              +J*exact_derivative(potential)-potential*J**2*Fraction(1,2))
    A=l1_upper(potential)/width
    B=l1_upper(residual)/width**3
    def pump(t):
        y=t/width
        return 630*float(step)/width*(y*(1-y))**4
    k=np.arange(2,129,dtype=float)
    a,b=fast_scattering(k,pump,width,3e-12)
    _,b2=fast_scattering(k,pump,width,5e-13)
    numbers=abs(b)**2
    shells=k*(k*k-1)*numbers/pi**2
    aprxA,aprxB=float(A),float(B)
    rows=[]
    for upper in (8,16,32,64,128):
        # Fractions A,B are outer bounds, printed values are rounded for display.
        tail=np.exp(2*aprxA/(upper+1))*aprxB**2/(128*pi**2*upper**2)
        rows.append(dict(highest_k=upper,sum_a4_density=float(shells[k<=upper].sum()),
                         rigorous_tail_bound_rounded=float(tail)))
    bound=np.exp(aprxA/k)*aprxB/(2*k)**3
    checks=dict(max_abs_occupation_refinement=float(max(abs(numbers-abs(b2)**2))),
                max_canonical_residual=float(max(abs(abs(a)**2-numbers-1))),
                max_beta_over_mode_bound=float(max(abs(b)/bound)),
                max_number_over_squeeze_bound=float(max(numbers)/np.sinh(float(step))**2))
    if (checks["max_abs_occupation_refinement"]>1e-10 or checks["max_canonical_residual"]>1e-9
            or checks["max_beta_over_mode_bound"]>1 or checks["max_number_over_squeeze_bound"]>1):
        raise ArithmeticError(checks)
    return dict(log_endpoint_ratio=float(step),support_duration=width,
                L1_potential_upper=aprxA,L1_second_residual_upper=aprxB,
                A_rational=str(A),B_rational=str(B),partial_sums=rows,checks=checks)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=dict(energy="late physical excess above constant-I out vacuum times a^4",
                sphere="unit radius harmonics k>=2 with 2(k^2-1) real transverse oscillators",
                sech=exact_sech(),compact=compact_pulse())
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
