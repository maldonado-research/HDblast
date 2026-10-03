"""Fabricated grid/operator benchmark; no checkpoint inputs or real bump calls."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import resource
import time

from flint import acb,arb,ctx
import verified_moments as v


def fixture_panels():
    # A fabricated degree24 polynomial in local x with conservative model error.
    # Both coefficient values and uncertainty are unrelated to the bump source.
    coefficients=tuple(v.rational_ball(Fraction((-1)**j, (j+1)*128**j))
                       for j in range(v.SOURCE_DEGREE+1))
    return tuple(v.PolynomialPanel(Fraction(-9,2)+Fraction(2*j+1,128),v.HALF_WIDTH,
                 coefficients,v.SOURCE_REMAINDER) for j in range(v.PANELS))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--nodes",type=int,default=8192)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    if args.nodes<=0:raise ValueError("Positive fabricated node count required")
    started=time.monotonic();ctx.prec=v.PRECISION_BITS
    family=v.EntireKernelFamily();panels=fixture_panels()
    setup=time.monotonic()-started
    maximum_radius=arb(0);digest=hashlib.sha256();sample=[]
    # Completely fabricated exact midpoint nodes over [0,256], equal weights.
    for j in range(args.nodes):
        k=Fraction(256*(2*j+1),2*args.nodes)
        result=v.forced_response(panels,k,acb(v.rational_ball(Fraction(1,1000000))),
                                 acb(0,v.rational_ball(Fraction(1,100000))),
                                 Fraction(1,10000),family)
        if not result["u"].is_finite() or not result["w"].is_finite():
            raise ValueError("Nonfinite fabricated enclosure")
        # Outward exact rational output is exercised for every fabricated mode.
        encoded=json.dumps({"u":v.complex_endpoints(result["u"]),
                            "w":v.complex_endpoints(result["w"])},sort_keys=True).encode()
        digest.update(encoded)
        for name in ("u","w"):
            for component in (result[name].real,result[name].imag):
                maximum_radius=maximum_radius.union(component.rad())
        if j in (0,args.nodes//2,args.nodes-1):sample.append({"index":j,"k":str(k),"bytes":len(encoded)})
        if j and j%2048==0:print(json.dumps({"completed":j,"nodes":args.nodes,"elapsed_seconds":time.monotonic()-started}),flush=True)
    receipt={"status":"PASS_FABRICATED_OPERATOR_BENCHMARK","scope":"one fabricated source, uniform64panels, all explicit fabricated nodes, additive global Duhamel; not contact/ledger/full12case driver",
        "physical_source_evaluations":0,"checkpoint_array_decodes":0,"fabricated_nodes":args.nodes,
        "momentum_interval":["0","256"],"panels":64,"source_degree":24,"phase_degree":96,"precision_bits":256,
        "setup_seconds":setup,"wall_seconds":time.monotonic()-started,
        "peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "maximum_radius_hull":v.rational_endpoints(maximum_radius),"serialized_endpoint_digest":digest.hexdigest(),"samples":sample,
        "source_sha256":hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))


if __name__=="__main__":main()
