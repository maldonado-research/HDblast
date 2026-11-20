"""Check the spherical-wall normal-pressure law against the full five-dimensional curvature."""
import math

import numpy as np


def metric_jets(point, geometry):
    """Five-metric, first derivatives, areal radius and its first two coordinate derivatives."""
    u, v, chi, theta, _ = point
    dtype = point.dtype
    g = np.zeros((5, 5), dtype=dtype)
    dg = np.zeros((5, 5, 5), dtype=dtype)
    dr = np.zeros(5, dtype=dtype)
    ddr = np.zeros((5, 5), dtype=dtype)
    if geometry == 'double_null':
        sigma = .13*u - .09*v + .07*u*v
        g[0, 1] = g[1, 0] = -np.exp(2*sigma)
        dg[:2, 0, 1] = dg[:2, 1, 0] = 2*g[0, 1]*np.array([.13+.07*v, -.09+.07*u])
        radius = 2 + v - u/3 + u*v/10 + u*u/20 + v*v/25
        dr[:2] = [-1/3+v/10+u/10, 1+u/10+2*v/25]
        ddr[:2, :2] = [[.1, .1], [.1, .08]]
    elif geometry == 'vaidya':
        radius = v
        dr[1] = 1
        mass = .7 - .18*u + .02*u*u
        g[0, 0] = -(1 + v*v - mass/(v*v))
        g[0, 1] = g[1, 0] = -1
        dg[0, 0, 0] = (-.18+.04*u)/(v*v)
        dg[1, 0, 0] = -2*v - 2*mass/v**3
    else:
        raise ValueError(geometry)
    s, c = np.sin(chi), np.cos(chi)
    t, d = np.sin(theta), np.cos(theta)
    sphere = np.array([1, s*s, s*s*t*t])
    angles = np.arange(2, 5)
    g[angles, angles] = radius**2*sphere
    dg[:2, angles, angles] = 2*radius*dr[:2, None]*sphere
    dg[2, 3, 3] = 2*radius**2*s*c
    dg[2, 4, 4] = 2*radius**2*s*c*t*t
    dg[3, 4, 4] = 2*radius**2*s*s*t*d
    return g, dg, radius, dr, ddr


def connection(point, geometry):
    g, dg, *_ = metric_jets(point, geometry)
    lower = dg.transpose(1, 0, 2) + dg.transpose(1, 2, 0) - dg
    return .5*np.einsum('ad,dbc->abc', np.linalg.inv(g), lower)


def einstein(point, geometry):
    # Differentiate the full five-dimensional connection by a complex step in each coordinate.
    gamma = connection(point, geometry)
    step = 1e-20
    deriv = np.stack([connection(point+1j*step*unit, geometry).imag/step for unit in np.eye(5)])
    ricci = (np.einsum('ccab->ab', deriv) - np.einsum('bcac->ab', deriv)
             + np.einsum('ccd,dab->ab', gamma, gamma) - np.einsum('cbd,dac->ab', gamma, gamma))
    g, _, radius, dr, ddr = metric_jets(point, geometry)
    inv = np.linalg.inv(g)
    scalar = np.sum(inv*ricci)
    return ricci - .5*scalar*g, g, inv, radius, dr, ddr - np.einsum('abc,a->bc', gamma, dr)


def projection(geometry, point, velocity):
    """Project full curvature and independently take the areal-radius Hessian along the observer."""
    G, g, inv, r, dr, hessian = einstein(np.array(point, dtype=float), geometry)
    if geometry == 'double_null':
        scale = math.sqrt(-2*g[0, 1])
        u, v = math.exp(velocity)/scale, math.exp(-velocity)/scale
        obs = np.array([u, v, 0., 0., 0.])
        normal = np.array([-u, v, 0., 0., 0.])
    else:
        x = math.sqrt(-g[0, 0]+velocity**2)
        udot = 1/(x+velocity)
        obs = np.array([udot, velocity, 0., 0., 0.])
        normal = np.array([-udot, x, 0., 0., 0.])
    checks = [normal @ g @ normal-1, obs @ g @ normal, obs @ g @ obs+1]
    mass = r*r*(1+r*r-dr @ inv @ dr)
    pressure = normal @ (G-6*g) @ normal
    hess = obs @ hessian @ obs/r
    rhs = -1-mass/r**4-pressure/3
    checks.append(hess-rhs)
    row = dict(radius=float(r), mass=float(mass), kappa2_Tnn=float(pressure),
               hessian_UU_over_radius=float(hess), predicted_hessian_over_radius=float(rhs))
    if geometry == 'vaidya':
        u = point[0]
        exact = -3*(-.18+.04*u)*udot**2/(2*r**3)
        row['null_dust_kappa2_Tnn'] = float(exact)
        checks.extend([pressure-exact, mass-(.7-.18*u+.02*u*u)])
    if max(map(abs, checks)) > 2e-12:
        raise ArithmeticError(f'{geometry} curvature or frame projection failed')
    return row


def calculate():
    coordinates = [(.11, .57, .81, 1.02, .4), (-.23, .63, 1.18, .74, .3)]
    double_null = [projection('double_null', point, boost)
                   for point in coordinates for boost in (-.4, .27)]
    vaidya = [projection('vaidya', (.4, 2.1, .81, .94, .3), speed) for speed in (-.31, .23)]
    streams = []
    for incident, reflected in [(.4, 0.), (.4, .13), (.4, .4)]:
        # Covariant null momenta in the wall's orthonormal time/outward-normal frame.
        plus, minus = np.array([-1., 1.]), np.array([-1., -1.])
        stress = incident*np.outer(plus, plus) + reflected*np.outer(minus, minus)
        pressure, net = float(stress[1, 1]), float(-stress[0, 1])
        if pressure != incident+reflected or net != incident-reflected:
            raise ArithmeticError('Counterstreaming pressure/flux sign failed')
        streams.append(dict(incident=incident, reflected=reflected, Tnn=pressure, net_flux=net))

    e = math.sqrt(1.25)-1
    pcrit = -(e+2*e**2)/(3*(1+e))
    vacuum, radiation = -e, e/3
    def acceleration(p):
        return (1+e)*(1-2*e-3*p)-1
    radiation_acc = acceleration(radiation)
    if not (vacuum < pcrit < radiation and radiation_acc < 0 and abs(acceleration(pcrit)) < 1e-14):
        raise ArithmeticError('Pressure threshold failed')
    samples = []
    for tau in (-1., 0., .5, 1.):
        a, v = 2*math.cosh(tau/2), math.sinh(tau/2)
        h = v/a
        if max(abs(h*h+a**-2-e*(2+e)), abs(acceleration(vacuum)-.25)) > 1e-14:
            raise ArithmeticError('Closed vacuum-fluid solution failed')
        samples.append(dict(proper_time=tau, radius=a, H=h))
    # Distinct ordinary interiors at rest: b_i^2=1+a^-2-m_i/a^4.
    density = .05
    radial = np.array([.99, 1.11])
    scaled_masses = 1.25-radial**2
    normal_stresses = np.array([.05, .08])  # kappa^2 T_nn / 3
    weights = (1/radial)/sum(1/radial)
    harmonic = 2/sum(1/radial)
    qbar = weights @ (scaled_masses+normal_stresses)
    necessary = -(density+2*density**2+qbar)/(3*(1+density))
    exact = (1-2*density-(1+qbar)/harmonic)/3
    if not (0 <= min(scaled_masses) and exact < necessary and abs(radial.mean()-1-density) < 1e-14):
        raise ArithmeticError('Asymmetric local stress threshold failed')
    exact_acceleration = (1-2*density-3*exact)*harmonic-1-qbar
    if abs(exact_acceleration) > 1e-14:
        raise ArithmeticError('Acceleration at asymmetric threshold differs from zero')
    unequal = dict(density=density, radial_curvatures=list(map(float, radial)),
                   mass_over_radius_fourth=list(map(float, scaled_masses)),
                   kappa2_Tnn_over_three=list(map(float, normal_stresses)),
                   Qbar=float(qbar), exact_critical_pressure=float(exact),
                   necessary_ceiling=float(necessary))
    return dict(units='AdS length=1; surface matter energy density and pressure divided by tuned tension',
                double_null=double_null, vaidya=vaidya, counterstreaming=streams, asymmetric_pressure=unequal,
                pressure_threshold=dict(radius_at_rest=2., extra_wall_density=e, critical_pressure=pcrit),
                vacuum_bounce=dict(pressure=vacuum, acceleration_over_radius=acceleration(vacuum),
                                   samples=samples),
                radiation_at_rest=dict(pressure=radiation, acceleration_over_radius=radiation_acc))
