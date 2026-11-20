"""Compare asymptotic reflection phases for two photon weights on the central fixed cut."""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import hankel1


THROAT_SLOPE = 5 / 9
TANGENTIAL_MOMENTUM = 1.


def warp(y):
    return -(y+(2/3)*math.log(math.cosh(y))+math.tanh(y)**2/6)/3


def impedance_derivative(y, selected):
    t = math.tanh(y)
    warp_y = -(1+t-t**3/3)/3
    return 1.5*warp_y if selected else -2*t+warp_y/2


def jost_reflection(masses, depth, selected, *, forward=True, rtol=3e-11):
    """Return reflection in an outgoing Hankel basis and its large-z phase."""
    masses = np.asarray(masses)
    s = 1.5 if selected else (1+4/THROAT_SLOPE)/2
    nu = s + .5
    asym = ((2/3)*math.log(2.)-1/6)/3
    r = math.exp(THROAT_SLOPE*depth-asym)/THROAT_SLOPE
    dr_dz = THROAT_SLOPE*r*math.exp(warp(depth))
    x = masses*r
    prefactor = np.sqrt(math.pi*x/(2*dr_dz))
    wave = prefactor*hankel1(nu,x)
    momentum = (dr_dz*masses*prefactor*hankel1(nu-1,x)
                -(s*THROAT_SLOPE+impedance_derivative(depth,selected))*math.exp(warp(depth))*wave)
    def rhs(y,state):
        u,v = state.reshape(2,-1)
        log_j_y = impedance_derivative(y,selected)
        dz_dy = math.exp(-warp(y))
        return np.r_[log_j_y*u + dz_dy*v,-log_j_y*v - dz_dy*masses**2*u]
    if forward:
        # Outgoing Jost determinant with the real solution (j_b,0) at the wall.
        initial = np.r_[np.ones_like(masses),np.zeros_like(masses)]
        sol = solve_ivp(rhs,(0.,depth),initial,method='DOP853',rtol=rtol,atol=rtol*1e-3)
        if not sol.success:
            raise ArithmeticError(sol.message)
        regular, regular_momentum = sol.y[:,-1].reshape(2,-1)
        determinant = regular*momentum-regular_momentum*wave
    else:
        sol = solve_ivp(rhs,(depth,0.),np.r_[wave,momentum],method='DOP853',rtol=rtol,atol=rtol*1e-3)
        if not sol.success:
            raise ArithmeticError(sol.message)
        determinant = sol.y[:,-1].reshape(2,-1)[1]
    # Multiplying the Hankel wave by exp(i(pi*nu/2+pi/4)) makes its
    # asymptotic phase m*r. The offset r-z is common to the two weights.
    plane_determinant = np.exp(1j*(math.pi*nu/2+math.pi/4))*determinant
    return -plane_determinant.conj()/plane_determinant


def relative(masses,depth,**kwargs):
    selected=jost_reflection(masses,depth,True,**kwargs)
    baseline=jost_reflection(masses,depth,False,**kwargs)
    return selected/baseline


def delays(masses, depth, step, **kwargs):
    offsets=np.r_[-2.,-1.,0.,1.,2.]*step
    grid=(masses[None,:]+offsets[:,None]).ravel()
    phase=relative(grid,depth,**kwargs).reshape(5,-1)
    centered=np.angle(phase*phase[2].conj())
    if np.max(np.abs(centered))>0.5:
        raise ArithmeticError('Frequency steps do not resolve the phase branch')
    derivative=(centered[0]-8*centered[1]+8*centered[3]-centered[4])/(12*step)
    return phase[2],derivative*np.hypot(masses,TANGENTIAL_MOMENTUM)/masses


def calculate():
    masses=np.array([.25,.5,1.,2.,4.])
    phase,delay=delays(masses,9.,.0008)
    finer_phase,finer_delay=delays(masses,10.,.0004,forward=False,rtol=1e-11)
    phases=np.angle(phase)
    if (max(abs(phase-finer_phase))>3e-7 or max(abs(delay-finer_delay))>2e-6 or
            max(abs(abs(phase)-1))>1e-12):
        raise ArithmeticError('Independent reflection propagation differs')
    # Compact bump in m describes the energy-flux spectrum of each incoming packet.
    def packet(central,half_width,order):
        nodes,weights=np.polynomial.legendre.leggauss(order)
        spectrum=np.exp(-1/(1-nodes**2))*weights
        _,result=delays(central+half_width*nodes,9.,.0008)
        return float(result@spectrum/spectrum.sum())
    packets=[]
    for central in (.5,1.,2.):
        half_width=.15*central
        value=packet(central,half_width,24)
        refined=packet(central,half_width,40)
        if abs(refined-value)>1e-6:
            raise ArithmeticError('Reflection delay quadrature differs')
        packets.append(dict(center_mass=central,half_width=half_width,mean_relative_delay=refined))
    return dict(reference='selected minus earlier weight on the same central cut',
                k=TANGENTIAL_MOMENTUM,asymptotic_reflected_energy_fraction=1.,
                rows=[dict(mass=float(m),relative_phase_radians=float(angle),
                           relative_delay=float(dt),independent_phase_radians=float(np.angle(check)),
                           independent_relative_delay=float(dcheck))
                      for m,angle,dt,check,dcheck in zip(masses,phases,delay,finer_phase,finer_delay)],
                packets=packets)
