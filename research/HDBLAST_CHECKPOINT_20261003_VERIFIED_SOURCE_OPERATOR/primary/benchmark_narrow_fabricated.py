"""Complete proposed operator-universe fabricated entry point, no real sources."""
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import resource
import time

from flint import acb,arb,ctx
import verified_moments as v

KS=(Fraction(0),Fraction(1,2**40),Fraction(1,2**12),Fraction(1,4),Fraction(1),
    Fraction(16),Fraction(64),Fraction(128),Fraction(256))


def serialized_radius(value):
    if isinstance(value,acb):
        return serialized_radius(value.real)+serialized_radius(value.imag)
    endpoints=v.rational_endpoints(value)
    return (Fraction(endpoints["upper"])-Fraction(endpoints["lower"]))/2


def main():
    started=time.monotonic();ctx.prec=v.PRECISION_BITS
    source_sha=hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest()
    family=v.EntireKernelFamily();cached={k:family.at(k) for k in KS}
    rows=[];whole=[];work=[];largest_panel=Fraction(0);largest_whole=Fraction(0)
    largest_panel_complete=Fraction(0);largest_whole_complete=Fraction(0)
    for fixture_index in (0,1):
        coefficients=tuple(v.rational_ball(Fraction(((-1)**j if fixture_index else 1),
                         (j+1)*128**j)) for j in range(25))
        panels=tuple(v.PolynomialPanel(Fraction(-9,2)+Fraction(2*j+1,128),v.HALF_WIDTH,
                     coefficients,v.SOURCE_REMAINDER) for j in range(64))
        work_panels=tuple(v.PolynomialPanel(p.center,p.half_width,
                         tuple(a/4 for a in p.coefficients),v.SOURCE_WORK_REMAINDER) for p in panels)
        work_total=arb(0)
        for panel_index,panel in enumerate(panels):
            for k in KS:
                moment=v.panel_moments(panel,cached[k])
                encoded={"M0":v.rational_endpoints(moment["M0"]),
                         "Mexp":v.complex_endpoints(moment["Mexp"]),
                         "Mdrift":v.complex_endpoints(moment["Mdrift"]),
                         "complete_serialized_radii":{name:str(serialized_radius(moment[name]))
                          for name in ("M0","Mexp","Mdrift")},
                         "source_model_disk_remainders":{name:str(moment[name]) for name in
                          ("source_M0_remainder","source_Mexp_remainder","source_Mdrift_remainder")},
                         "integrated_kernel_model_bounds":{"E":str(family.e_tail),"Q":str(family.q_tail)}}
                rows.append({"fabricated_fixture":fixture_index,"panel":panel_index,"k":str(k),**encoded})
                largest_panel_complete=max(largest_panel_complete,*(serialized_radius(moment[name])
                          for name in ("M0","Mexp","Mdrift")))
                for value in (moment["M0"],moment["Mexp"].real,moment["Mexp"].imag,
                              moment["Mdrift"].real,moment["Mdrift"].imag):
                    largest_panel=max(largest_panel,Fraction(str(value.rad().upper().fmpq())))
            work_moment=v.panel_moments(work_panels[panel_index],cached[Fraction(0)])["M0"]
            work_total+=work_moment
            work.append({"fabricated_fixture":fixture_index,"panel":panel_index,
                         "Lg_integral":v.rational_endpoints(work_moment)})
        for k in KS:
            moment=v.forced_response(panels,k,acb(0),acb(0),Fraction(1),family)
            whole.append({"fabricated_fixture":fixture_index,"k":str(k),
                          "M0":v.rational_endpoints(moment["M0"]),
                          "Mexp":v.complex_endpoints(moment["Mexp"]),
                          "Mdrift":v.complex_endpoints(moment["Mdrift"]),
                          "Lg_integral":v.rational_endpoints(work_total),
                          "complete_serialized_radii":{name:str(serialized_radius(moment[name]))
                          for name in ("M0","Mexp","Mdrift")},
                          "source_model_M0_Mexp_disk_remainder":str(v.SOURCE_REMAINDER),
                          "source_model_Mdrift_disk_remainder":str(v.SOURCE_REMAINDER/2)})
            largest_whole_complete=max(largest_whole_complete,serialized_radius(work_total),
                          *(serialized_radius(moment[name]) for name in ("M0","Mexp","Mdrift")))
            for value in (moment["M0"],moment["Mexp"].real,moment["Mexp"].imag,
                          moment["Mdrift"].real,moment["Mdrift"].imag,work_total):
                largest_whole=max(largest_whole,Fraction(str(value.rad().upper().fmpq())))
    if len(rows)!=1152 or len(whole)!=18 or len(work)!=128:
        raise ValueError("Incomplete proposed fabricated universe")
    if largest_whole_complete>Fraction(1,10**26):
        raise ValueError("Fabricated operator radius exceeds proposed1e-26 gate")
    if hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest()!=source_sha:
        raise ValueError("Source bytes changed during fabricated entry point")
    data={"panel_rows":rows,"whole_rows":whole,"work_rows":work}
    output=Path("FABRICATED_NARROW_VALUES.json")
    output.write_text(json.dumps(data,sort_keys=True)+"\n")
    receipt={"status":"PASS_COMPLETE_FABRICATED_OPERATOR_UNIVERSE",
             "physical_source_evaluations":0,"checkpoint_array_decodes":0,"uid":os.getuid(),
             "panels":64,"fabricated_sources":2,"k":[str(x) for x in KS],
             "panel_rows":len(rows),"whole_rows":len(whole),"work_rows":len(work),
             "precision_bits":256,"source_degree":24,"phase_degree":96,
             "largest_panel_component_radius":str(largest_panel),
             "largest_whole_component_radius":str(largest_whole),
             "largest_panel_complete_serialized_radius":str(largest_panel_complete),
             "largest_whole_complete_serialized_radius":str(largest_whole_complete),
             "gate":"1e-26 whole unscaled complete complex L1 radius: sum of real/imaginary outward endpoint half-widths",
             "wall_seconds":time.monotonic()-started,
             "peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
             "source_sha256":source_sha,"values_sha256":hashlib.sha256(output.read_bytes()).hexdigest(),
             "limits":"No real bump evaluations or pressure/contact ledger calculation. Internal benchmark; future freeze/official outer resource receipt remains required."}
    Path("FABRICATED_NARROW_RECEIPT.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":main()
