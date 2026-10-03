"""Fabricated-only positive and negative tests for candidate verified moments."""
from fractions import Fraction
import json
from math import factorial
from pathlib import Path
import resource
import time

from flint import acb, arb, arb_series, ctx
import verified_moments as v


def require(condition, message):
    if not condition:
        raise ValueError(message)


def constant_panels(value=Fraction(7, 11), remainder=Fraction(0)):
    return tuple(v.PolynomialPanel(Fraction(-9, 2) + Fraction(2*j+1, 128), v.HALF_WIDTH,
                  (v.rational_ball(value),) + tuple(arb(0) for _ in range(v.SOURCE_DEGREE)), remainder)
                 for j in range(v.PANELS))


def full_direct_oracle(k, value):
    # Independent adaptive validated integration of an entire fabricated kernel.
    omega = acb(0, v.rational_ball(2*k))
    c = v.rational_ball(value)
    mexp = acb.integral(lambda t, analytic: c * (omega * (1-t)).exp(), 0, 1,
                       rel_tol=arb((1,-190)), abs_tol=arb((1,-190)))
    if k == 0:
        mdrift = acb(v.rational_ball(value/2))
    else:
        mdrift = acb.integral(lambda t, analytic: c * ((omega*(1-t)).exp()-1)/omega,
                             0, 1, rel_tol=arb((1,-190)), abs_tol=arb((1,-190)))
    return mexp, mdrift


def main():
    started=time.monotonic()
    ctx.prec = v.PRECISION_BITS
    family=v.EntireKernelFamily()
    tests=[]
    def ok(name, condition):
        require(condition, name)
        tests.append({"name":name,"status":"PASS"})
    # Rational polynomial and integrated kernels at zero phase are exact.
    e,q=family.at(0)
    ok("zero_phase_entire_E",all(e[j].contains(v.rational_ball(v.j_monomial(j))) for j in range(25)))
    ok("zero_phase_entire_Q",all(q[j].contains(v.rational_ball(v.mixed_monomial(j,1))) for j in range(25)))
    panels=constant_panels()
    for k in (Fraction(0),Fraction(1,2**60),Fraction(1,256),Fraction(1),Fraction(64),Fraction(256)):
        result=v.forced_response(panels,k,acb(0),acb(0),Fraction(1),family)
        if k==Fraction(1,2**60):
            # Exact entire closed form evaluated at deliberately increased precision.
            with ctx.workprec(512):
                omega=acb(0,v.rational_ball(2*k))
                phase=omega.exp()
                me=v.rational_ball(Fraction(7,11))*(phase-1)/omega
                mu=v.rational_ball(Fraction(7,11))*(phase-1-omega)/(omega*omega)
        else:
            me,mu=full_direct_oracle(k,Fraction(7,11))
        ok(f"fabricated_constant_Mexp_k={k}",result["Mexp"].overlaps(me))
        ok(f"fabricated_constant_Mdrift_k={k}",result["Mdrift"].overlaps(mu))
        ok(f"fabricated_constant_state_k={k}",result["w"].overlaps(-me) and result["u"].overlaps(-mu))
    # Coefficient uncertainty and source-model remainder widen the result.
    uncertain=list(panels)
    uncertain[0]=v.PolynomialPanel(uncertain[0].center,uncertain[0].half_width,
                    (v.inflate_real(v.rational_ball(Fraction(7,11)),Fraction(1,1000000)),)
                    +tuple(arb(0) for _ in range(24)),Fraction(1,100000))
    narrowed=v.forced_response(panels,0,acb(0),acb(0),Fraction(1),family)
    widened=v.forced_response(uncertain,0,acb(0),acb(0),Fraction(1),family)
    ok("coefficient_ball_and_model_remainder_preserved",widened["Mexp"].real.contains(narrowed["Mexp"].real))
    # A fabricated source bias is detected by non-overlap, not fit differences.
    biased=v.forced_response(constant_panels(Fraction(7,11)+Fraction(1,1000)),0,acb(0),acb(0),Fraction(1),family)
    ok("negative_biased_source",not biased["Mexp"].overlaps(narrowed["Mexp"]))
    # Flipped forcing sign gives a disjoint exact synthetic result.
    ok("negative_forcing_sign",not (-narrowed["u"]).overlaps(narrowed["u"]))
    phase_fixture=v.forced_response(panels,1,acb(0),acb(0),Fraction(1),family)
    ok("negative_phase_sign",not phase_fixture["Mexp"].conjugate().overlaps(phase_fixture["Mexp"]))
    # A deliberately short phase polynomial must retain its positive remainder.
    short=v.EntireKernelFamily(degree=0,phase_degree=2)
    z=acb(0,4)
    direct_e=acb.integral(lambda x, analytic:(z*(1-x)).exp(),-1,1,
                         rel_tol=arb((1,-190)),abs_tol=arb((1,-190)))
    correct_e,_=short.at(256)
    ok("short_phase_with_remainder_contains_oracle",correct_e[0].contains(direct_e))
    ok("negative_phase_remainder_omission",not short.exp_polynomials[0](z).overlaps(direct_e))
    previous_cap=ctx.cap
    ctx.cap=27
    try:
        eta=arb_series([arb(-4),v.rational_ball(v.HALF_WIDTH)],prec=27)
        fabricated_high_jet=arb_series([arb(0)]*26+[arb(1)],prec=27)
        g,work=v.forcing_series_coefficients(eta,fabricated_high_jet)
        exact_g24=-Fraction(26*25,1)/v.HALF_WIDTH**2
        ok("fabricated_high_source_jet_keeps_degree24",g[24].contains(v.rational_ball(exact_g24)))
        ok("fabricated_work_source_degree24",work[24].contains(v.rational_ball(exact_g24/Fraction(4))))
        ok("fabricated_high_source_zero_lower_degrees",all(g[j].contains(0) for j in range(24)))
    finally:ctx.cap=previous_cap
    # Reject invalid geometry, negative remainders, and unproved k cap.
    for name,callback in (
        ("negative_remainder",lambda:v.PolynomialPanel(Fraction(0),v.HALF_WIDTH,(arb(1),),Fraction(-1))),
        ("unproved_phase_cap",lambda:family.at(257)),
        ("physical_callback_without_freeze",lambda:v.future_registered_source_panel(Fraction(-4),"positive_B")),
    ):
        try:callback()
        except (ValueError,RuntimeError):ok(name,True)
        else:ok(name,False)
    # Exact rational serialization includes the interval and does not mutate it.
    a=v.inflate_real(v.rational_ball(Fraction(1,3)),Fraction(1,10000))
    endpoints=v.rational_endpoints(a)
    ok("outward_rational_serialization",Fraction(endpoints["lower"])<=Fraction(str(a.lower().fmpq()))
       and Fraction(endpoints["upper"])>=Fraction(str(a.upper().fmpq())))
    receipt={"status":"PASS_FABRICATED_MOMENT_TESTS","physical_source_evaluations":0,
       "checkpoint_array_decodes":0,"tests":tests,"precision_bits":ctx.prec,
       "source_degree":v.SOURCE_DEGREE,"phase_degree":v.PHASE_DEGREE,
       "wall_seconds":time.monotonic()-started,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
       "limits":"Moment prototype only; no contact/pressure/full-ledger certificate or physical run."}
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":main()
