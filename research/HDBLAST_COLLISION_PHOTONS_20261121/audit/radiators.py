"""Fold two closed massive radiators over each measured collision without an energy lookup table."""
import argparse
import json
from pathlib import Path
import runpy

import numpy as np
from numpy.polynomial.legendre import leggauss


POLE=8300.**2
BANDS=((0.,.001),(.001,.01),(.01,.10),(.10,.25),(.25,1.))


def born(s):
    """Angular average of the spin-summed/divided-by-four Born amplitude, up to a common constant."""
    return (s+2)**2/((s-POLE)**2+.30**2*(s+2)**2*(1.-4./s))


def tensor_numerator(source,mass,x):
    """Numerator of x times hep-ph/0505236 Eq. (5), including continuation to negative x."""
    beta2=1.-4./source
    beta_pair=np.sqrt(1.-4./mass)
    log_pair=np.log(mass/4.)+2.*np.log1p(beta_pair)
    value=((1.+beta2)*(3.-beta2)-2.*x*(3.-beta2)+2.*x*x)*log_pair
    value-=2.*((3.-beta2)*(1.-x)+x*x)*beta_pair
    return value/(np.pi*(3.-beta2))


def energy_spectra(s,q,x):
    """Return x dN/(alpha dx) for initial and final radiation. q=s(1-x) is computed directly at recoil."""
    b2=1.-4./s
    br=np.sqrt(1.-4./q)
    fsr=tensor_numerator(s,q,x)/np.sqrt(b2)
    isr=(q/s)*br/b2*born(q)/born(s)*tensor_numerator(q,s,-s*x/q)
    return np.stack((isr,fsr),axis=-1)


def rates(s,unit,weights):
    """Integrate angle-inclusive count and energy moments over each event's entire allowed spectrum."""
    s=s[:,None]
    endpoint=1.-4./s
    table=np.zeros((len(s),len(BANDS),2,2))
    for i,(lo,hi) in enumerate(BANDS):
        def add(q,x,dx):
            f=energy_spectra(s,q,x)
            table[:,i,0,:]+=np.sum(.5*dx[...,None]*f,axis=1)
            if lo>0.:
                table[:,i,1,:]+=np.sum((dx/x)[...,None]*f,axis=1)
        if hi<=.25:
            if lo==0.:
                x=hi*endpoint*unit
                dx=hi*endpoint*weights
            else:
                span=np.log(hi/lo)
                x=lo*endpoint*np.exp(span*unit)
                dx=span*x*weights
            add(s*(1.-x),x,dx)
        else:
            # Integrate in sqrt(log(s'/4)) to resolve both the massive threshold and the pole.
            high=s*(1.-lo*endpoint)
            breaks=[np.full_like(s,4.)]
            breaks.extend(np.minimum(high,v) for v in (POLE*.5,POLE,POLE*2))
            breaks.append(high)
            for left,right in zip(breaks[:-1],breaks[1:]):
                y0=np.sqrt(np.log(left/4.))
                dy=np.sqrt(np.log(right/4.))-y0
                y=y0+dy*unit
                q=4.*np.exp(y**2)
                x=1.-q/s
                add(q,x,2.*y*(q/s)*dy*weights)
    return table


def measured(path,n):
    with np.load(path) as d:
        s=d["s"]
        scale=float(d["weight"])*d["a"]*d["E_in"].sum(axis=1)/float(d["gas"])
        alpha=3.*float(d["charge_to_mass"])**2*float(d["daughter_mass_L"])**2/(4*np.pi*float(d["action"])*float(d["K"]))
        tracks=2*int(d["pairs"])
    nodes,weights=leggauss(n)
    unit=(nodes+1.)/2.
    weights=weights/2.
    sums=np.zeros((len(BANDS),2,2))
    per_event=np.empty_like(s)
    for start in range(0,len(s),768):
        chunk=slice(start,min(start+768,len(s)))
        table=rates(s[chunk],unit,weights)
        table[:,:,0,:]*=(alpha*scale[chunk])[:,None,None]
        table[:,:,1,:]*=alpha/tracks
        per_event[chunk]=table[:,:,0,:].sum(axis=(1,2))
        sums+=table.sum(axis=0)
    # Cumulative shares refer to incoming encounter mass (before radiation) and to yield.
    mask=s>=np.quantile(s,.99)
    total=per_event.sum()
    bands=[dict(fraction=cuts,isr_energy_over_gas=float(data[0,0]),fsr_energy_over_gas=float(data[0,1]),
                photons_per_particle=float(data[1].sum()) if cuts[0] else None) for cuts,data in zip(BANDS,sums)]
    return dict(tape=path.name,events=len(s),bands=bands,total_energy_over_gas=float(total),
                hard_count_per_particle=float(sums[2:,1].sum()),isr_energy_share=float(sums[:,0,0].sum()/total),
                above_one_percent_energy_share=float(sums[2:,0].sum()/total),
                upper_one_percent_mass_threshold=float(np.sqrt(np.quantile(s,.99))),
                upper_one_percent_mass_share_of_emitted_energy=float(per_event[mask].sum()/total),
                upper_one_percent_mass_share_of_pair_energy=float(scale[mask].sum()/scale.sum()),
                largest_single_event_share_of_expected_energy=float(per_event.max()/total))


def differential_checks(module):
    nodes,weights=leggauss(48)
    unit=(nodes+1)/2
    rows=[]
    for s in (6.25,25.,8300.**2,1.46782e11):
        residual=[0.,0.]
        for end_fraction in (.01,.43,.95,.9999):
            x=end_fraction*(1-4./s)
            closed=energy_spectra(s,s*(1-x),x)
            spin=np.array([x*module["pair_rate"](s,x,unit,weights/2,leg) for leg in (False,True)])
            residual=np.maximum(residual,np.abs(spin/closed-1.))
        rows.append(dict(s=s,maximum_relative_isr_difference=float(residual[0]),maximum_relative_fsr_difference=float(residual[1])))
    if max(max(r["maximum_relative_isr_difference"],r["maximum_relative_fsr_difference"]) for r in rows)>1.e-6:
        raise ValueError("Spinor spectra differ from the published radiator or its crossed expression")
    return rows


def run(opts):
    hard=runpy.run_path(str(opts.hard/"emission.py"))
    with opts.reference.open() as src:
        baseline={row["tape"]:row for row in json.load(src)["collision_histories"]}
    checks=differential_checks(hard)
    with np.load(max(opts.artifacts.glob("t7.05000-p14-*.npz"))) as data:
        largest=float(data["s"].max())
    max_tail={}
    for n in (24,48,96):
        nodes,weights=leggauss(n)
        a=rates(np.array([largest]),(nodes+1)/2,weights/2)[0]
        max_tail[f"n{n}_mean_energy_fraction_over_alpha"]=float(a[:,0].sum())
    with np.errstate(all="raise"):
        direct=sum(hard["band"](largest,lo,hi,40)["mean_fraction_of_pair_energy_over_alpha"] for lo,hi in BANDS)
    max_tail.update(s=largest,spinor_mean_energy_fraction_over_alpha=float(direct),
                    direct_spinor_fractional_discrepancy=float(direct/max_tail["n96_mean_energy_fraction_over_alpha"]-1.))
    if abs(max_tail["direct_spinor_fractional_discrepancy"])>1.e-5:
        raise ValueError("Unresolved highest observed collision mass")
    rows=[]
    for path in sorted(opts.artifacts.glob("*.npz")):
        ref=baseline[path.name]
        row=measured(path,48)
        coarse=measured(path,24)
        differences=[abs((a["isr_energy_over_gas"]+a["fsr_energy_over_gas"])/b["fraction_of_endpoint_gas"]-1.)
                     for a,b in zip(row["bands"],ref["bands"])]
        row["max_band_energy_relative_to_spinor"]=max(differences)
        row["total_energy_relative_to_spinor"]=row["total_energy_over_gas"]/ref["total_energy_over_gas"]-1.
        row["hard_count_relative_to_spinor"]=row["hard_count_per_particle"]/ref["photons_above_one_percent_per_charged_particle"]-1.
        row["half_nodes_energy_relative_change"]=coarse["total_energy_over_gas"]/row["total_energy_over_gas"]-1.
        if max(row["max_band_energy_relative_to_spinor"],abs(row["half_nodes_energy_relative_change"]))>1.e-4:
            raise ValueError(f"Integrated radiator discrepancy in {path.name}: {row}")
        rows.append(row)
    if len(rows)!=8 or len(baseline)!=8:
        raise ValueError("Expected eight matching collision ensembles")
    opts.output.parent.mkdir(exist_ok=True,parents=True)
    opts.output.write_text(json.dumps(dict(differential=checks,largest_mass_fine_late_tape=max_tail,runs=rows),indent=2,sort_keys=True,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ("artifacts","hard","reference","output"):
        parser.add_argument("--"+arg,type=Path,required=True)
    run(parser.parse_args())
