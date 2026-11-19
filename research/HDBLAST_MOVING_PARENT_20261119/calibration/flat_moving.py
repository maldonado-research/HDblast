"""Weak, real, Gaussian-switched scalar decay at fixed parent momentum in Minkowski space."""
import json
import math

SQRT_PI = math.sqrt(math.pi)


def gaussian_moments(lower):
    """Integrate z**j exp(-z**2)/sqrt(pi) from lower to infinity, j = 0,1,2,3."""
    edge = math.exp(-lower*lower) / (2*SQRT_PI)
    zeroth = .5*math.erfc(lower)
    return zeroth, edge, .5*zeroth+lower*edge, (lower*lower+1)*edge


def closed_branch(boost, duration, sign):
    energy = math.hypot(1, boost)
    # Stable expression for (E-q)/m even when q/m is large.
    lower = -duration/(energy+boost) if sign == 1 else duration*(energy+boost)
    i0, i1, i2, i3 = gaussian_moments(lower)
    n = i0 + 2*sign*energy/duration*i1 + i2/duration**2
    w = sign*i0 + (1/energy+2*energy)/duration*i1 + 3*sign/duration**2*i2 + i3/(energy*duration**3)
    return n, w


def simpson(func, lo, hi, panels):
    if lo == hi:
        return 0.
    step = (hi-lo)/panels
    edge = func(lo)+func(hi)
    odd = math.fsum(func(lo+i*step) for i in range(1,panels,2))
    even = math.fsum(func(lo+i*step) for i in range(2,panels,2))
    return step*(edge+4*odd+2*even)/3


def direct(boost, duration, panels):
    energy = math.hypot(1, boost)
    gap = duration/(energy+boost)
    high = duration*(energy+boost)

    def ann(z, field):
        above = (z+gap)/duration
        x = boost+above
        measure = above*(above+2*boost)*math.exp(-z*z)/SQRT_PI
        if field == "energy":
            return x/energy*measure
        if field == "work":
            return z/(duration*energy)*measure
        return measure

    def stim(v, field):
        # Shift the threshold, factor out exp(-high**2), and stretch its boundary layer.
        t = v/(1+high)
        above = t/duration
        x = boost+above
        measure = above*(above+2*boost)*math.exp(-t*(t+2*high))/(SQRT_PI*(1+high))
        if field == "energy":
            return x/energy*measure
        return measure

    def annihilation(field):
        return simpson(lambda z: ann(z,field),max(-gap,-14.),14.,panels)

    def scaled_stimulation(field):
        return simpson(lambda v: stim(v,field),0.,64.,panels)

    added = scaled_stimulation("number")
    log_added = (math.log(added)-high**2)/math.log(10)
    nstim = added*math.exp(-high**2)
    dstim = scaled_stimulation("energy")*math.exp(-high**2)
    nann = annihilation("number")
    dann = annihilation("energy")

    # Pair +/- detunings in the overlap and integrate the two positive remainders.
    # The pair gives 4 z^2/duration^2 in the dimensionless EL normalization.
    work_overlap = simpson(lambda z: 4*z*z*math.exp(-z*z)/(duration**2*SQRT_PI),
                           0.,min(gap,14.),panels)
    work_outer = simpson(lambda z: ann(z,"work"),min(gap,14.),14.,panels)
    w_positive = work_overlap + work_outer + nstim + dstim

    # Parseval integrated over the full line would include a negative spacelike
    # polynomial. Restore it over -high < z < -gap by positive quadrature.
    def addback(z):
        below = (-gap-z)/duration
        return below*(2*boost-below)*math.exp(-z*z)/SQRT_PI
    lo = max(-high,-14.)
    add = simpson(addback,lo,-gap,panels) if -gap > lo else 0.
    return dict(annihilation=nann, stimulated=nstim, log10_stimulated=log_added,
                pair_count=nann+nstim, parent_depletion=nann-nstim, daughter_energy=dann+dstim,
                work_by_detuning=w_positive, spacelike_addback=add)


def test(boost, duration):
    e = math.hypot(1,boost)
    ann, da = closed_branch(boost,duration,1)
    stim, ds = closed_branch(boost,duration,-1)
    count, lost = ann+stim, ann-stim
    daughter = da+ds
    work = daughter-lost
    checks = [direct(boost,duration,p) for p in (4096,8192)]
    for check in checks:
        for key, exact in (("annihilation",ann), ("pair_count",count), ("parent_depletion",lost),
                           ("daughter_energy",daughter), ("work_by_detuning",work),
                           ("spacelike_addback",count-1-.5/duration**2)):
            if not math.isclose(exact,check[key],abs_tol=5e-11,rel_tol=2e-8):
                raise ValueError((boost,duration,key,exact,check[key]))
        if max(stim,check["stimulated"])>1.e-13 and not math.isclose(stim,check["stimulated"],rel_tol=2e-8):
            raise ValueError("Stimulated branch")
    if work <= 0 or not math.isclose(work,daughter-lost,rel_tol=2e-14):
        raise ValueError("Work closure")
    return dict(q_over_m=boost, mT=duration, threshold_gap_times_T=duration/(e+boost),
                inclusive_count=count, signed_parent_depletion=lost, daughter_energy=daughter,
                net_switching_work=work,
                pair_longitudinal_momentum_over_EL=boost/e*lost,
                equal_opposite_parent_transfer_over_EL=-boost/e*lost,
                log10_stimulated_count_over_L=checks[-1]["log10_stimulated"],
                positive_work_quadrature_error=checks[-1]["work_by_detuning"]-work,
                timelike_threshold_addback=checks[-1]["spacelike_addback"])


if __name__ == "__main__":
    cases = [test(q,t) for q in (0.,1.,4.,16.) for t in (1.,4.,16.,64.)]
    print(json.dumps(dict(units={"number":"L = y^2 m^2 sqrt(pi) T/(8 pi E)",
                                 "energy_momentum":"E L", "switch":"exp(-t^2/(2 T^2))"},
                          rows=cases), indent=2, allow_nan=False))
