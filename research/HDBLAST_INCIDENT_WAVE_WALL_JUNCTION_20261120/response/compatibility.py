"""Israel and scalar junction checks for an excursion in a fixed mirrored wall."""
import math

from scipy.integrate import quad


def geometry(y):
    phi = -math.tanh(y)
    log_a = -(y + 2 * math.log(math.cosh(y)) / 3 + phi**2 / 6) / 3
    a_y_over_a = -(1 - phi + phi**3 / 3) / 3
    phi_y = -1 / math.cosh(y)**2
    return math.exp(log_a), a_y_over_a, phi_y


def trajectory(tau, center, distance, width):
    x = tau / width
    bump = math.exp(-x*x / 2)
    return (center + distance*bump, -distance*x*bump / width,
            distance*(x*x - 1)*bump / width**2)


def derivative(fn, x, step):
    return (fn(x-2*step) - 8*fn(x-step) + 8*fn(x+step) - fn(x+2*step)) / (12*step)


def case(center, distance, width):
    def motion(tau):
        return trajectory(tau, center, distance, width)

    def scaled_enthalpy(tau):
        _, speed, accel = motion(tau)
        return 2 * accel / math.hypot(1, speed)

    def time_component(tau):
        y, speed, _ = motion(tau)
        return math.hypot(1, speed) / geometry(y)[0]

    def scaled_density(tau):
        y, speed, _ = motion(tau)
        return -6 * geometry(y)[1] * math.hypot(1, speed)

    # Rederive both K components from the metric, connection and tangent.
    junction_error = 0.
    conservation_error = 0.
    for index in range(-40, 41):
        tau = width * index / 10
        y, speed, _ = motion(tau)
        gamma = math.hypot(1, speed)
        a, ay, phi_y = geometry(y)
        u_time = gamma/a
        step = 3e-4 * width
        metric_y = derivative(lambda yy: geometry(yy)[0]**2, y, 3e-4)
        accel_time = derivative(time_component, tau, step) + metric_y/a**2*u_time*speed
        accel_y = derivative(lambda t: motion(t)[1], tau, step) + metric_y*u_time**2/2
        k_time = -a*speed*accel_time + gamma*accel_y
        k_space = gamma*metric_y/(2*a*a)
        enthalpy = scaled_enthalpy(tau)
        junction_error = max(junction_error, abs(2*(k_time-k_space)-enthalpy))
        energy_rate = derivative(scaled_density, tau, step) + 3*ay*speed*enthalpy
        conservation_error = max(conservation_error, abs(energy_rate-2*phi_y**2*speed*gamma))
    # Fix the null affine normalization by requiring d(tau)/d(lambda) -> 1.
    asymptotic_a = geometry(center)[0]
    def weighted_enthalpy(tau):
        a = geometry(motion(tau)[0])[0]
        return scaled_enthalpy(tau)*asymptotic_a/a

    def anec_after_parts(tau):
        y, speed, _ = motion(tau)
        a, ay, _ = geometry(y)
        return 2*asymptotic_a/a * ay*speed*math.asinh(speed)

    def full_integral(integrand):
        cuts = [-math.inf, -width, width, math.inf]
        return sum(quad(integrand, low, high, epsabs=1e-13, epsrel=1e-12)[0]
                   for low, high in zip(cuts, cuts[1:]))

    anec = full_integral(weighted_enthalpy)
    anec_independent = full_integral(anec_after_parts)
    if abs(anec-anec_independent) > 1e-11 or (distance != 0 and anec >= 0):
        raise ArithmeticError("Affine null energy integration by parts failed")
    budget = 4 * math.asinh(abs(distance)*math.exp(-.5) / width)
    integrals = []
    for sign in (1, -1):
        def integrand(t):
            return max(0., sign*scaled_enthalpy(t))
        integrals.append(full_integral(integrand))
    # Ratio of normal scalar derivative to the tuned static junction value.
    peak_scalar_ratio = math.hypot(1, distance*math.exp(-.5)/width)
    if max(abs(value-budget) for value in integrals) > 1e-11:
        raise ArithmeticError("Energy budgets disagree with the rapidity integral")
    if max(junction_error, conservation_error) > 2e-9:
        raise ArithmeticError("Independent coordinate or surface-energy check failed")
    return dict(y0=center, displacement=distance, width=width,
                center_scaled_rho_plus_p=scaled_enthalpy(0.),
                exact_each_integrated_magnitude=budget,
                numerical_positive_integral=integrals[0],
                numerical_negative_magnitude=integrals[1],
                scaled_affine_null_energy=anec,
                scaled_affine_null_energy_by_parts=anec_independent,
                peak_normal_scalar_ratio=peak_scalar_ratio,
                max_coordinate_junction_error=junction_error,
                max_surface_energy_identity_error=conservation_error)


def calculate():
    return dict(normal="increasing y in the retained copy",
                normalization="rho, p and their integrals are multiplied by kappa_5 squared",
                examples=[case(.3, .1, 1.), case(.3, -.1, 1.),
                          case(.3, .4, .5), case(.3, 0., 1.)])
