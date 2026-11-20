"""Compare the exact spherical Maxwell scattering matrix with the frequency-domain mesh."""
import math

import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve
from scipy.special import jv, yv, jvp, yvp


def basis(radius, omega):
    x = omega*radius
    return np.array([[jv(2,x), yv(2,x)], [omega*jvp(2,x), omega*yvp(2,x)]])


def exact_channel(radius, gap, omega, even, kinetic):
    kappa = 4/radius**2
    # Even fields are continuous, odd fields vanish at the common wall.
    value = 1. if even else 0.
    slope = kinetic*(omega**2-kappa)/2 if even else 1.
    coefficients = np.linalg.solve(basis(radius,omega),[value,slope])
    port_q,port_slope = basis(radius-gap,omega)@coefficients
    speed = -1j*omega*port_q
    incoming = (speed-port_slope)/2
    outgoing = (speed+port_slope)/2
    return outgoing/incoming, -1j*omega*value/incoming


def mesh_channel(mesh, radius, gap, omega, even, kinetic, elements):
    mass, angular, gradient, spatial = mesh(radius,gap,elements)
    center = len(mass)//2
    port_weight = (radius-gap)/radius
    correction = (-omega**2*mass).astype(complex)
    correction[center] += kinetic*(4/radius**2-omega**2)
    correction[[0,-1]] -= 1j*omega*port_weight
    drives = np.zeros(len(mass))
    drives[0] = 1.
    drives[-1] = 1. if even else -1.
    speed = -1j*omega*spsolve((spatial+diags(correction)).tocsc(),2*port_weight*drives)
    return speed[0]-1.,speed[center]


def pulse_transmission(radius, gap, pulse, fft_length, step=1/128):
    """Flux spectrum of the compact incident speed times the exact power transmission."""
    start,end,signal = pulse(1.,20.)
    sampled = np.array([signal(t) for t in np.arange(start,end,step)])
    spectrum = abs(np.fft.rfft(sampled,fft_length))**2
    spectrum[1:-1] *= 2
    omega = 2*math.pi*np.fft.rfftfreq(fft_length,step)[1:]
    # The fundamental matrix inverse uses W(J2,Y2)=2/(pi*r). Its real,
    # common normalization cancels from both reflected characteristic ratios.
    jb, yb = jv(2,omega*radius), yv(2,omega*radius)
    djb,dyb = omega*jvp(2,omega*radius),omega*yvp(2,omega*radius)
    xs = omega*(radius-gap)
    js,ys,djs,dys = jv(2,xs), yv(2,xs), omega*jvp(2,xs), omega*yvp(2,xs)
    def reflection(value, slope):
        cj,cy = dyb*value-yb*slope, jb*slope-djb*value
        speed, deriv = -1j*omega*(cj*js+cy*ys), cj*djs+cy*dys
        return (speed+deriv)/(speed-deriv)
    even = reflection(1.,(omega**2-4/radius**2)/2)
    odd = reflection(0.,1.)
    power = abs((even-odd)/2)**2
    # The finite-curvature zero-frequency limit has zero transmission.
    return float(np.dot(spectrum[1:],power)/spectrum.sum())


def pulse_check(run, pulse, radius):
    short = pulse_transmission(radius,1.,pulse,32768)
    spectral = pulse_transmission(radius,1.,pulse,131072,step=1/256)
    row = run(radius=radius, width=1., channel='one', g=0., cap=20.,
              gap=1., elements=10, horizon=84., tol=2e-10)
    t = row['readout']
    transmitted,incident = t['outgoing_right'],t['incident_left']
    remaining = sum(t[k] for k in ('surface_electric','surface_magnetic','bulk_electric',
                                   'bulk_angular','bulk_radial_gradient','scalar'))
    # The unmeasured outgoing energy on the right is bounded by all energy
    # still inside the domain after the incoming pulse has ended.
    uncertainty = 3e-7
    if abs(spectral-short)>1e-7 or not (
            transmitted/incident-uncertainty <= spectral <= (transmitted+remaining)/incident+uncertainty):
        raise ArithmeticError('Broadband exact Maxwell transmission disagrees with the pulse')
    return dict(radius=radius,readout_time=t['t'],incident=incident,outgoing_right=transmitted,
                leftover_energy=remaining,spectral_fraction=spectral,
                shorter_transform_fraction=short)


def calculate(mesh, run, pulse):
    gap,kinetic = 1.,1.
    points = []
    # omega=2/R is the zero of the surface reactive term. Lower and higher
    # frequencies test both signs and the magnetic curvature term.
    for radius in (2.,4.,16.):
        for omega in (.45,2/radius,1.7,4.):
            phases, gains = [],[]
            errors = []
            for even in (True, False):
                reflected, wall = exact_channel(radius,gap,omega,even,kinetic)
                phases.append(reflected)
                gains.append(wall)
                coarse = mesh_channel(mesh,radius,gap,omega,even,kinetic,4)
                fine = mesh_channel(mesh,radius,gap,omega,even,kinetic,12)
                coarse_error = max(abs(coarse[0]-reflected),abs(coarse[1]-wall))
                fine_error = max(abs(fine[0]-reflected),abs(fine[1]-wall))
                errors.append((coarse_error,fine_error))
            # A one-sided incident wave is half even and half odd.
            left, right = (phases[0]+phases[1])/2, (phases[0]-phases[1])/2
            residual = max(abs(abs(z)**2-1) for z in phases)
            if (max(x[1] for x in errors)>1e-8 or residual>2e-12 or
                    abs(abs(left)**2+abs(right)**2-1)>2e-12):
                raise ArithmeticError('Radial Maxwell phase or flux verification failed')
            points.append(dict(radius=radius,omega=omega,
                               even_phase=float(np.angle(phases[0])),odd_phase=float(np.angle(phases[1])),
                               even_wall_speed_gain=float(abs(gains[0])),
                               one_side_reflectance=float(abs(left)**2),
                               one_side_transmittance=float(abs(right)**2),
                               coarse_mesh_error=float(max(x[0] for x in errors)),
                               refined_mesh_error=float(max(x[1] for x in errors))))
    result = dict(gap=gap,surface_kinetic=kinetic,coarse_elements=4,fine_elements=12,points=points)
    result['pulse_transmission']=[pulse_check(run,pulse,r) for r in (2.,4.)]
    return result
