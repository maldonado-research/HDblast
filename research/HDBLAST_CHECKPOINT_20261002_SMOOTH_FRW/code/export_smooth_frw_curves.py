"""Curate executed integrated curves without recomputing scientific quantities."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def digest(path):
    sha=hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda:stream.read(2**20),b''):
            sha.update(chunk)
    return sha.hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    parser.add_argument('output',type=Path)
    args=parser.parse_args()
    summary=json.loads((args.run/'summary.json').read_text())
    if summary['status'] not in ('PASS','FAIL') or 'static_negative_control' not in summary:
        raise ValueError('Only a completed matrix can be curated')
    args.output.mkdir(parents=True,exist_ok=False)
    manifest=dict(run_status=summary['status'],source_run=str(args.run),
                  provenance=summary['provenance'],curated_files={},retained_raw_artifacts={})
    for filename in ('summary.json','provenance.json'):
        path=args.output/filename
        path.write_bytes((args.run/filename).read_bytes())
        manifest['curated_files'][filename]=dict(sha256=digest(path),bytes=path.stat().st_size)
    for directory in sorted(args.run.iterdir()):
        if not directory.is_dir() or not (directory/'diagnostics.json').exists():
            continue
        outdir=args.output/directory.name
        outdir.mkdir()
        diagnostic=outdir/'diagnostics.json'
        diagnostic.write_bytes((directory/'diagnostics.json').read_bytes())
        with np.load(directory/'arrays.npz',allow_pickle=False) as archived:
            curves={key:archived[key] for key in archived.files
                    if key in ('t','a','H','x','xp') or key.startswith('K')}
        if not all(np.all(np.isfinite(value)) for value in curves.values()):
            raise ValueError('Nonfinite archived integrated curve')
        curves_path=outdir/'integrated_curves.npz'
        np.savez_compressed(curves_path,**curves)
        for path in (diagnostic,curves_path):
            manifest['curated_files'][str(path.relative_to(args.output))]=dict(sha256=digest(path),bytes=path.stat().st_size)
        for filename in ('arrays.npz','physical_modes.npz'):
            source=directory/filename
            manifest['retained_raw_artifacts'][str(source.relative_to(args.run))]=dict(
                sha256=digest(source),bytes=source.stat().st_size,
                disposition='Retained locally; exactly reproducible with registered code and protocol')
    (args.output/'CURATION_MANIFEST.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(output=str(args.output),status=summary['status'],curated_files=len(manifest['curated_files']),
                         retained_raw_files=len(manifest['retained_raw_artifacts'])),allow_nan=False))


if __name__=='__main__':
    main()
