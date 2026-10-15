"""Exercise bootstrap rejections in disposable copies; never alter scientific source."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys


def require(condition, reason):
    if not condition:
        raise RuntimeError(reason)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--fixture',type=Path,required=True);parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();src=args.fixture.resolve();out=args.output_dir.resolve()
    require(not out.exists(),'Fresh scratch directory required');out.mkdir(parents=True)
    manifest_pin=hashlib.sha256((src/'MANIFEST.json').read_bytes()).hexdigest();rows=[]
    cases={'bad_external_pin':'External manifest pin differs','changed_source':'Frozen bytes differ:',
           'missing_source':'Real source inside capsule required:','extra_file':'Frozen tree contains missing or extra files',
           'source_symlink':'Real source inside capsule required:','duplicate_record':'Noncanonical, reserved or duplicate source path',
           'stale_output':'New or empty output required'}
    for mode in ['normal','optimized']:
        flags=['-I','-B']+(['-O'] if mode=='optimized' else [])
        for case,reason in cases.items():
            work=out/mode/case;cap=work/src.name;shutil.copytree(src,cap);pin=manifest_pin;replay=work/'replay'
            if case=='bad_external_pin':pin='0'*64
            elif case=='changed_source':(cap/'waveform_intervals.py').write_bytes((cap/'waveform_intervals.py').read_bytes()+b'\n# altered\n')
            elif case=='missing_source':(cap/'waveform_intervals.py').unlink()
            elif case=='extra_file':(cap/'unlisted.txt').write_text('unlisted')
            elif case=='source_symlink':
                (cap/'waveform_intervals.py').unlink();(cap/'waveform_intervals.py').symlink_to(src/'waveform_intervals.py')
            elif case=='duplicate_record':
                d=json.loads((cap/'MANIFEST.json').read_bytes());d['files'].append(d['files'][0]);raw=(json.dumps(d,indent=2)+'\n').encode();(cap/'MANIFEST.json').write_bytes(raw);pin=hashlib.sha256(raw).hexdigest()
            elif case=='stale_output':replay.mkdir();(replay/'stale.json').write_text('{"status":"PASS"}\n')
            p=subprocess.run([sys.executable,*flags,str(cap/'reproduce.py'),'--manifest-sha256',pin,'--output-dir',str(replay)],capture_output=True,text=True,timeout=30)
            require(p.returncode==1 and reason in p.stderr,'Wrong rejection: '+case)
            require((not replay.exists()) if case!='stale_output' else list(replay.iterdir())==[replay/'stale.json'],'Calculation ran after rejected bootstrap')
            (work/'execution.log').write_text(p.stdout+p.stderr)
            rows.append({'mode':mode,'case':case,'exit_code':p.returncode,'required_reason':reason,'rejected_before_calculations':True})
        p=subprocess.run([sys.executable,*flags,str(src/'verify_manifest.py'),'--expected',manifest_pin],capture_output=True,text=True,timeout=30)
        require(p.returncode==0 and json.loads(p.stdout)['status']=='PASS','Valid source verifier failed')
    receipt={'status':'PASS','scope':'Tool-integrity controls on disposable fixture copies, separate from scientific fresh calculations',
             'valid_verifier_modes':2,'rejections':len(rows),'runner_sha256':hashlib.sha256((src/'reproduce.py').read_bytes()).hexdigest(),
             'verifier_sha256':hashlib.sha256((src/'verify_manifest.py').read_bytes()).hexdigest(),'fixture_manifest_sha256':manifest_pin,'checks':rows}
    (out/'INTEGRITY_CONTROLS.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
