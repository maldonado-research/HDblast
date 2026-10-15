"""Assemble a fresh tree from the explicit public-science directory whitelist."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil

ROOT=Path(__file__).resolve().parent.parent
PUB=ROOT/'publication'
TOP=['README.md','INERTIAL_SHELL_WAVEFORM.md','inertial_error_bounds.py',
     'compute_inertial_waveform.py','check_inertial_waveform.py']
COMPUTED=['INERTIAL_WAVEFORM_ENCLOSURE.json','INERTIAL_WAVEFORM_CELLS.csv',
          'CHECK_INERTIAL_WAVEFORM.json','INERTIAL_SHELL_WAVEFORM.png']
PREVIOUS=['UNIFORM_FIRST_VARIATION_WAVEFORM.md','archive/ANALYTIC_SHELL_BOX_DELTA_1E3.json',
          'archive/PHYSICAL_GAIN_AUDIT.md','archive/PHYSICAL_GAIN_ENCLOSURE.json',
          'archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md','sol/UNIFORM_WAVEFORM_CERTIFICATE.md']
LITERATURE=['INERTIAL_WAVEFORM_PRIMARY_LITERATURE.md','INERTIAL_WAVEFORM_SOURCE_LEDGER.json',
            'LITERATURE_PROVENANCE_CHECK.json','BOUNDED_SEARCH_RECEIPT.json']


def require(ok,reason):
    if not ok:raise RuntimeError(reason)


def selected():
    files={name:ROOT/name for name in TOP}
    files.update({'computed/'+n:ROOT/'computed'/n for n in COMPUTED})
    files.update({'inputs/previous/'+n:ROOT/'inputs/previous'/n for n in PREVIOUS})
    files['inputs/PRIOR_SOURCE_PINS.json']=ROOT/'inputs/PRIOR_SOURCE_PINS.json'
    files.update({'literature/'+n:ROOT/'literature'/n for n in LITERATURE})
    safe='millennium/HDBLAST_METHODS_TRANSFER_PUBLIC_SAFE.md';files[safe]=ROOT/safe
    for folder in ['sol','astra','referee','archive','archive/inputs']:
        require((ROOT/folder).is_dir(),'Missing public source directory: '+folder)
        for path in sorted((ROOT/folder).iterdir()):
            require(not path.is_symlink(),'Source symlink rejected')
            if path.is_dir():continue
            require(path.is_file() and path.suffix in {'.py','.md','.json','.txt'},'Unreviewed top-level source type: '+str(path))
            files[folder+'/'+path.name]=path
    for name in ['reproduce.py','verify_manifest.py']:
        files[name]=PUB/'replay-tools'/name
    for name in ['REPLAY_ROSTER.json','requirements.txt','REPRODUCE.md']:
        files[name]=PUB/name
    return files


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=ROOT/'assembled-source')
    parser.add_argument('--receipt',type=Path,default=PUB/'SOURCE_ASSEMBLY_RECEIPT.json')
    parser.add_argument('--inventory-only',action='store_true')
    args=parser.parse_args();out=args.output_dir.absolute();sources=selected();snapshots={}
    for name,path in sorted(sources.items()):
        require(path.is_file() and not any(p.is_symlink() for p in [path,*path.parents]),'Missing or nonregular source: '+name)
        require(not ({p.lower() for p in Path(name).parts}&{'internal','private','__pycache__'}),'Unpublishable source path')
        snapshots[name]=path.read_bytes()
    if args.inventory_only:
        print(json.dumps({'status':'PREPARED_NOT_ASSEMBLED','source_files':len(snapshots),'paths':sorted(snapshots)},indent=2));return
    require(not out.exists() and not any(p.is_symlink() for p in [out,*out.parents]),'Fresh real assembly output required')
    require(not args.receipt.resolve().is_relative_to(out.resolve()),'Assembly receipt must stay outside source')
    out.mkdir(parents=True)
    for name,raw in snapshots.items():
        target=out/name;target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(raw)
        require(target.read_bytes()==raw,'Assembly copy differs: '+name)
    require(all(path.read_bytes()==snapshots[name] for name,path in sources.items()),'Original changed while assembling')
    receipt={'status':'PASS','source_files':len(snapshots),
             'files':[{'path':n,'bytes':len(b),'sha256':sha256(b).hexdigest()} for n,b in snapshots.items()],
             'scope':'Byte-preserving explicit whitelist assembly only. Nested scratch replays, internal paper bodies and private project audits are excluded; the separate publication-scope audit checks exact private-file copies.'}
    args.receipt.write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({'status':'PASS','source_files':len(snapshots),'output':str(out)},indent=2))


if __name__=='__main__':main()
