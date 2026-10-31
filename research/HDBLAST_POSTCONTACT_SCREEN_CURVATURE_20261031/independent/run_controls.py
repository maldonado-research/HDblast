#!/usr/bin/env python3
"""Actual subprocess corruption controls, including optimized Python."""
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

sys.dont_write_bytecode=True
HERE=Path(__file__).resolve().parent
SOURCE=Path('/workspace/hdblast-absolute-publish-31c')


def need(ok,message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    dest=HERE/'controls'
    need(not dest.exists(),'FRESH_CONTROLS_DIRECTORY')
    dest.mkdir()
    checker=(HERE/'check_independent.py').read_text()
    pins=json.loads((HERE/'SOURCE_PINS.json').read_text())
    oprel='research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/inputs/linear_operator.py'
    maprel='research/HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/projection/GEOMETRIC_COEFFICIENT_MAPS.json'
    cases=[
      ('source_bytes_corrupted','SOURCE_AUTHENTICATION',None),
      ('lift_sign_reversed','LIFT_POSITIVE_INWARD_R_DERIVATIVE',('LIFT_SIGN = -1','LIFT_SIGN = 1')),
      ('Kx_lower_image_sign_reversed','KX_ENDPOINT',('KX_SECOND_SIGN = 1','KX_SECOND_SIGN = -1')),
      ('Neumann_image_parity_reversed','KX_ENDPOINT',('REFLECTION_SIGN = 1','REFLECTION_SIGN = -1')),
      ('TT_bulk_cross_coupling','PRIMITIVE_TT_CLOSURE',None),
      ('TT_source_row_deleted','COMPLETE_78_TERM_TT_MAP',None),
    ]
    jobs=[]
    for name,tag,edit in cases:
        work=dest/name; work.mkdir()
        privatepins=dict(pins); text=checker; source=SOURCE
        if name in ('source_bytes_corrupted','TT_bulk_cross_coupling','TT_source_row_deleted'):
            source=work/'source';source.mkdir()
            for rel in pins:
                target=source/rel;target.parent.mkdir(parents=True,exist_ok=True)
                shutil.copyfile(SOURCE/rel,target)
            if name=='source_bytes_corrupted':
                with (source/oprel).open('a') as stream:
                    stream.write('\n# actual corrupt source bytes\n')
            elif name=='TT_bulk_cross_coupling':
                target=source/oprel
                old=target.read_text(); token='    return out\n\ndef assemble'
                need(old.count(token)==1,'UNIQUE_OPERATOR_MUTATION')
                target.write_text(old.replace(token,"    out['Aphi'][21][0]=Q(1)\n    return out\n\ndef assemble"))
                privatepins[oprel]=sha(target)
            else:
                target=source/maprel;d=json.loads(target.read_text())
                d['maps']['metric_forcing']['ab'][6].pop()
                target.write_text(json.dumps(d,sort_keys=True,indent=2)+'\n')
                privatepins[maprel]=sha(target)
        if edit:
            old,new=edit;need(text.count(old)==1,'UNIQUE_CHECKER_MUTATION:'+name)
            text=text.replace(old,new)
        pinbytes=(json.dumps(privatepins,sort_keys=True,indent=2)+'\n').encode()
        (work/'SOURCE_PINS.json').write_bytes(pinbytes)
        # Semantic source fixtures deliberately authorize their changed pins,
        # then exercise the independent mathematical guard after importing.
        original_pinmap=hashlib.sha256((HERE/'SOURCE_PINS.json').read_bytes()).hexdigest()
        text=text.replace(original_pinmap,hashlib.sha256(pinbytes).hexdigest())
        mutant=work/'check_independent.py';mutant.write_text(text)
        for optimized in (False,True):
            report=work/('OPTIMIZED.json' if optimized else 'NORMAL.json')
            cmd=[sys.executable]+(['-O'] if optimized else [])+[str(mutant),'--source-root',str(source),'--output',str(report)]
            jobs.append((name,tag,optimized,cmd,report))
    def execute(job):
        name,tag,optimized,cmd,report=job
        run=subprocess.run(cmd,capture_output=True,text=True,timeout=90)
        message=run.stdout+run.stderr
        need(run.returncode!=0 and not report.exists() and tag in message,
             'NEGATIVE_CONTROL_ACCEPTED_OR_WRONG_FAILURE:'+name+':'+str(optimized)+':'+message[-500:])
        return {'case':name,'optimized':optimized,'exit_code':run.returncode,
                'expected_rejection':tag,'report_created':report.exists(),
                'stderr_last_line':run.stderr.strip().splitlines()[-1]}
    with ThreadPoolExecutor(max_workers=3) as pool:
        rows=list(pool.map(execute,jobs))
    normal=(HERE/'INDEPENDENT_NORMAL.json').read_bytes()
    optimized=(HERE/'INDEPENDENT_OPTIMIZED.json').read_bytes()
    need(normal==optimized,'NORMAL_OPTIMIZED_REPORT_MISMATCH')
    result={'status':'ALL_REAL_NEGATIVE_CONTROLS_REJECTED',
            'cases':rows,'subprocesses':len(rows),
            'positive_normal_optimized_byte_identical':True,
            'positive_report_sha256':hashlib.sha256(normal).hexdigest(),
            'scope':'Each failed subprocess wrote no report. Source-byte corruption uses original pins; semantic source mutants explicitly repin their altered fixture to reach structural guards.'}
    (dest/'CONTROLS.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'])


if __name__=='__main__':
    main()
