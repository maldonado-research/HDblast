"""Upper bound for vacuum photons with bounded total variation of the gauge coefficient."""
import math


def total_weight(n):
    return n*(n-1)*(n+1)*(n+2)//4 if n>=2 else 0


def budgets(initial_radius, h0_over_mpl, delta=.1):
    def coefficient(n):
        return h0_over_mpl**2*total_weight(n)/(3*math.pi**2*initial_radius**2)
    def bound(n):
        return coefficient(n)*math.sinh(abs(delta))**2
    n=2
    while bound(n)<.5:
        n+=1
    rows=[]
    for cutoff in (.05,.1,.25,.5,1.):
        k=math.floor(initial_radius*cutoff/h0_over_mpl)
        weight=sum(l*(l*l-1) for l in range(2,k+1))
        if weight!=total_weight(k):
            raise ArithmeticError('Closed-form spectral weight differs')
        c=coefficient(k)
        rows.append(dict(initial_frequency_cap_over_Mpl=cutoff,maximum_shell=k,
                         fraction_ceiling=bound(k),minimal_variation_for_half=math.asinh(math.sqrt(.5/c)) if c else None))
    return dict(total_variation_log_I=abs(delta),first_shell_cutoff_not_excluding_half=n,
                its_bounce_frequency_over_Mpl=n*h0_over_mpl/initial_radius,
                below_threshold_ceiling=bound(n-1),at_threshold_ceiling=bound(n),rows=rows)
