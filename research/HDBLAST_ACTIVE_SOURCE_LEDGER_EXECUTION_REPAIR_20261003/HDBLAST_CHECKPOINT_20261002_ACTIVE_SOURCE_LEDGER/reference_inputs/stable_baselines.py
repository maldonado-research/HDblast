"""Exact rationalized minimal fixed-r baseline differences, before momentum integration.

These are baseline bare minus complete subtraction identities, not substitutes
for either operator or a Ward definition. verify_stable_baselines.py proves
the identities symbolically, including the relation omega^2=k^2+2a0^2.
"""
import numpy as np
LD = np.longdouble


def horner(argument, coefficients):
    out = np.zeros_like(argument)+LD(coefficients[-1])
    for coefficient in reversed(coefficients[:-1]):
        out = out*argument+LD(coefficient)
    return out


def baseline_differences(k, omega, a0):
    ratio = k/omega
    difference = 2*a0*a0/(omega+k)
    variance = difference**3*horner(ratio,(16,23,21,15,5))/(32*k*omega**3)
    density = difference**4*horner(ratio,(384,493,-204,-1188,-1940,-1990,-868,308,420,105))/(1024*k*omega**2)
    pressure = -difference**4*horner(ratio,(384,1536,2944,3039,252,-5740,-15260,-19250,-8652,3612,4620,1155))/(3072*k*omega**2)
    return variance,density,pressure
