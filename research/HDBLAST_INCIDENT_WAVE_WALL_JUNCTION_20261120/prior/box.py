"""Independent finite-volume spectral and driven-oscillator tests on a conformal box.

Use table(phi) for the resolution study, grid_result(phi, length, spacing, k_values,
widths) for a single box. Energy coefficients are E_cont / (d^2 E_in) for one
retained half-line and unit warp at the wall. Lengths and times use L0.
"""
import math

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.linalg import eigh_tridiagonal


def primitive(y):
    logcosh = np.logaddexp(y, -y) - np.log(2.)
    return y + (2/3)*logcosh + np.tanh(y)**2/6


def radial_system(phi, length, spacing, selected=False):
    """Mass-lumped linear elements for -(w u')' = m^2 w u, w=exp(A) G."""
    if abs(phi) >= 1 or min(length, spacing) <= 0:
        raise ValueError("A finite wall and a positive box size and mesh are required")
    intervals = math.ceil(length / spacing)
    dz = length / intervals
    zhalf = np.linspace(0., length, 2*intervals + 1)
    wall = -np.arctanh(phi)
    reference = primitive(wall)

    def proper_position(z, y):
        return np.exp((reference - primitive(y))/3)

    trajectory = solve_ivp(proper_position, (0., length), [wall], method="DOP853",
                           rtol=4e-13, atol=2e-15, t_eval=zhalf)
    if not trajectory.success:
        raise ArithmeticError(trajectory.message)
    y = trajectory.y[0]
    warp = np.exp((reference - primitive(y))/3)
    weight = warp**3 if selected else warp / np.cosh(y)**4
    mass = dz * weight[::2]
    mass[[0, -1]] *= .5
    springs = weight[1::2] / dz
    diagonal = (np.r_[0., springs] + np.r_[springs, 0.]) / mass
    off = -springs / np.sqrt(mass[1:] * mass[:-1])
    return mass, diagonal, off, y[-1]


def row_key(k, width):
    return f"{k:g}/{width:g}"


def grid_result(phi, length, spacing, k_values=(.25, 1., 2.), widths=(.5, 1., 2.),
                time_check=False, selected=False):
    if min(*k_values, *widths) <= 0:
        raise ValueError("Positive incidence and pulse width required")
    mass, diagonal, off, far_y = radial_system(phi, length, spacing, selected=selected)
    # The last retained temporal frequency is at least nine inverse pulse widths
    # from the incident one, so the omitted finite-box weight is exponentially small.
    upper = max(k_values) + 9/min(widths)
    eigen, vectors = eigh_tridiagonal(diagonal, off, select="v",
                                      select_range=(-1., upper**2), lapack_driver="stebz")
    zero_vector = np.sqrt(mass / mass.sum())
    overlap = abs(vectors[:, 0] @ zero_vector)
    pole = vectors[0, 0]**2 / mass[0]
    if abs(eigen[0]) > 2e-8 or 1 - overlap > 2e-9 or abs(pole*mass.sum()-1) > 2e-8:
        raise ArithmeticError("Finite-box Neumann zero mode is inconsistent")
    squared_masses = eigen[1:]
    wall_weight = vectors[0, 1:]**2 / mass[0]
    norm = (1+phi)**2*(2-phi)/3
    far_phi = -math.tanh(far_y)
    exterior = (1+far_phi)**2*(2-far_phi)/3
    gb = (1-phi**2)**2
    if selected:
        # A zero-mode normalisation is an integral of G(y) dy in proper distance.
        # Fix the local gauge field convention by setting the dielectric to one at the wall.
        yb = -np.arctanh(phi)
        ab = primitive(yb)
        def integrand(y):
            return np.exp(2*(ab-primitive(y))/3)
        norm = quad(integrand,yb,np.inf,epsabs=3e-14,epsrel=3e-13)[0]
        exterior = quad(integrand,far_y,np.inf,epsabs=1e-16,epsrel=3e-13)[0]
        gb = 1.
    box_integral = norm - exterior
    factor = gb**2/(4*norm)
    result = {}
    for k in k_values:
        omega = np.sqrt(k*k + squared_masses)
        for width in widths:
            def kernel(freq):
                return 2*np.pi*(width*freq)**2*np.exp(-(width*freq)**2)
            weight = kernel(omega-k) + kernel(omega+k)
            # The full discrete wall evaluation squared is 1/mass[0]. Each
            # term of the Gaussian spectrum decreases beyond its own maximum.
            tail = 2*factor*kernel(np.hypot(k,upper)-k)/mass[0]
            result[row_key(k, width)] = dict(energy=factor*float(wall_weight @ weight),
                                             discarded_box_tail_bound=float(tail))
            if selected:
                # A reflecting finite box has a normalizable constant whose wall
                # residue 1/N_L exceeds the half-line residue 1/N. Correct JUST
                # this known geometrical difference; report the raw box separately.
                # X(0)=0 for a closed pulse, while X(2k) is generically nonzero.
                pole_bias = factor*kernel(2*k)*exterior/(box_integral*norm)
                record = result[row_key(k,width)]
                record['excess_box_pole_contribution'] = float(pole_bias)
                record['corrected_energy'] = record['energy'] + pole_bias
    time_result = None
    if time_check:
        # q_j'' + omega_j^2 q_j = f(t), vanishing data at t=-10T. Two
        # quadratures of the incident massless wave implement phase averaging.
        # The common boundary factor is reinserted when the modal energies are read.
        k = 1.
        width = 1.
        omega = np.sqrt(k*k+squared_masses)
        n = len(omega)
        def forced(t,state):
            q,v = state.reshape(2,2,n)
            derivative = -t/width**2*np.exp(-.5*(t/width)**2)
            force = derivative*np.array([math.cos(k*t), math.sin(k*t)])[:,None]
            return np.r_[v.ravel(),(force-omega**2*q).ravel()]
        sol = solve_ivp(forced,(-10*width,10*width),np.zeros(4*n),method="DOP853",
                        rtol=2e-10,atol=5e-13)
        if not sol.success:
            raise ArithmeticError(sol.message)
        q,v = sol.y[:,-1].reshape(2,2,n)
        measured = factor*float(wall_weight @ (2*np.sum(v**2+omega**2*q**2,axis=0)))
        reference_energy = result[row_key(k,width)]['energy']
        if abs(measured/reference_energy-1) > 5e-8:
            raise ArithmeticError("Integrated finite-pulse energy disagrees with its Gaussian spectrum")
        time_result = dict(k=k,width=width,energy=measured,
                           relative_difference=abs(measured/reference_energy-1))
    return dict(length=length,spacing=length/(len(mass)-1),nonzero_modes_kept=len(squared_masses),
                continuum_cutoff=upper,zero_mass_squared=float(eigen[0]),
                zero_overlap=float(overlap),pole_error=abs(float(pole*mass.sum())-1),
                mass_quadrature_relative_error=float(mass.sum()/box_integral-1),
                time_check=time_result,entries=result)


def table(phi, k_values=(.25,1.,2.), widths=(.5,1.,2.), selected=False):
    """Return box-size and mesh studies; extrapolation assumes second-order spatial error."""
    if 1. not in k_values or 1. not in widths:
        raise ValueError("Include k=T=1 for the independent forced-work check")
    if selected:
        # The selected law has a slower-decaying throat than sech^4(y).
        # Compare at fixed mesh through three sizes for the low-mass spectrum.
        sizes=((48.,.04),(96.,.04),(96.,.02),(96.,.01),(192.,.04))
    else:
        sizes=((24.,.04),(48.,.04),(48.,.02),(48.,.01))
    rows=[grid_result(phi,length,spacing,k_values,widths,time_check=(i==0),selected=selected)
          for i,(length,spacing) in enumerate(sizes)]
    estimate={}
    for name in rows[-1]['entries']:
        label='corrected_energy' if selected else 'energy'
        vals=[r['entries'][name][label] for r in rows]
        fine=(4*vals[3]-vals[2])/3
        coarser=(4*vals[2]-vals[1])/3
        estimate[name]=dict(energy=fine,relative_finest_change=abs(vals[3]/fine-1),
                            relative_extrapolation_change=abs(coarser/fine-1),
                            relative_box_change=abs(vals[1]/vals[0]-1))
        if selected:
            raw=[r['entries'][name]['energy'] for r in rows]
            estimate[name].update(relative_larger_box_change=abs(vals[4]/vals[1]-1),
                                  raw_relative_box_change=abs(raw[1]/raw[0]-1),
                                  raw_relative_larger_box_change=abs(raw[4]/raw[1]-1))
    return dict(phi=phi,grids=rows,estimates=estimate)
