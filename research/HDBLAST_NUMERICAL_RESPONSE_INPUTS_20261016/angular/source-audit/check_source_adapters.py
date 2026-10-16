#!/usr/bin/env python3
"""Portable normal/-O source replay; all mathematics uses the standard library."""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys
import tempfile
sys.path.insert(0,str(Path(__file__).resolve().parent))
import source_adapters as a
import check_source_identities as angular
import check_scalar_gauge_adapters as scalar_gauge
import check_shell_initial_adapters as shell_initial

CONTROLS={
 'angular_tt_deleted':'SOURCE_Y_HESS',
 'qc_time_sign':'SOURCE_QP_IDENTITY',
 'scalar_potential_omitted':'SOURCE_QPHI_IDENTITY',
 'shell_tt_omitted':'SOURCE_SHELL_TENSOR',
 'normal_gauge_sign':'SOURCE_NORMAL_GAUGE_JOIN',
 'initial_mixed_conversion_omitted':'SOURCE_INITIAL_MIXED_VELOCITY',
 'einstein_scalar_stress_halved':'SOURCE_MATTER_YY',
}
def mutate(name):
    if name=='angular_tt_deleted':a.YHESS=(*a.YHESS[:5],0,0)
    elif name=='qc_time_sign':
        original=a.connection_contractions
        def broken(m,b):
            out=original(m,b);out['0']=a.scale(-1,out['0']);return out
        a.connection_contractions=scalar_gauge.connection_contractions=broken
    elif name=='scalar_potential_omitted':
        original=a.scalar_bulk_source
        def broken(m,s,b):return a.add(original(m,s,b),a.scale(b['Vppp']*s['s']**2,a.Y2))
        a.scalar_bulk_source=scalar_gauge.scalar_bulk_source=broken
    elif name in ('shell_tt_omitted','normal_gauge_sign'):
        original=a.shell_sources
        def broken(**kw):
            out=original(**kw)
            if name=='shell_tt_omitted':out['QKab']=(*out['QKab'][:5],0,0)
            else:out['GN']=a.scale(-1,out['GN'])
            return out
        a.shell_sources=broken
    elif name=='initial_mixed_conversion_omitted':
        original=a.initial_frame_data
        def broken(**kw):
            out=original(**kw);coef=kw['rho']*(-48*kw['u']*kw['uy']-192*kw['a']*kw['u']**2)
            out['t_velocities']['y0']=a.add(out['t_velocities']['y0'],a.scale(coef,a.Y2));return out
        a.initial_frame_data=broken
    elif name=='einstein_scalar_stress_halved':
        original=a.einstein_matter_source
        def broken(m,s,b):
            out=original(m,s,b);out['yy']=a.add(out['yy'],a.scale(s['sy']**2,a.Y2));return out
        a.einstein_matter_source=broken
    else:raise ValueError('UNRECOGNIZED_CONTROL')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--run-controls',action='store_true')
    parser.add_argument('--mutation',choices=tuple(CONTROLS))
    args=parser.parse_args()
    output=args.output.resolve();source_root=Path(__file__).resolve().parents[1]
    if output.exists():raise ValueError('OUTPUT_MUST_BE_NEW')
    if output==source_root or source_root in output.parents:raise ValueError('OUTPUT_MUST_BE_EXTERNAL')
    if args.mutation:mutate(args.mutation)
    angular_count=angular.validate();scalar_gauge.validate();shell=shell_initial.validate()
    controls={}
    if args.run_controls:
        with tempfile.TemporaryDirectory(prefix='hdblast-source-controls-') as directory:
            for name,marker in CONTROLS.items():
                dest=Path(directory)/(name+'.json')
                command=[sys.executable,'-I','-S']
                if sys.flags.optimize:command.append('-O')
                command+=['-B',str(Path(__file__).resolve()),'--output',str(dest),'--mutation',name]
                result=subprocess.run(command,text=True,capture_output=True,timeout=60)
                if result.returncode==0 or dest.exists() or marker not in result.stderr:raise ValueError('CONTROL_NOT_REJECTED_'+name)
                controls[name]={'rejected':True,'marker':marker,'output_absent':True}
    paths=[Path(__file__).resolve().parent/name for name in ('source_adapters.py','source_laurent.py','check_source_identities.py','check_scalar_gauge_adapters.py','check_shell_initial_adapters.py','check_source_adapters.py','SOURCE_INVENTORY.md')]
    paths.append(source_root/'angular_exact.py')
    pins={str(path.relative_to(source_root)):{'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size} for path in paths}
    report={'status':'PASS','scope':{'angular_identities':'universal exact polynomial quotient','radial_time_background_identities':'independent Laurent symbols on rank-complete angular witnesses','complete_sources':['QC lower covector','Qphi','Einstein matter terms','QK/Qp physical shell','QC_y/GN normal shell','second initial frame values and velocities'],'remaining':['full geometric Einstein source QE/QEred radial coefficient adapter','actual evolved master and gauge jets']},'angular_identity_count':angular_count,'scalar_angular_witness_rank':3,'vector_angular_witness_rank':2,**shell,'controls':controls,'pins':pins}
    with output.open('x',encoding='utf-8') as handle:json.dump(report,handle,indent=2,sort_keys=True);handle.write('\n')
    print('PASS: source-adapter replay and requested controls.')
if __name__=='__main__':main()
