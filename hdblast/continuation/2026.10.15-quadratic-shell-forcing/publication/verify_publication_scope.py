"""Reject private/internal publication inputs without exposing their contents."""
import argparse
import hashlib
import json
from pathlib import Path


def require(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def main():
    task=Path(__file__).resolve().parent.parent
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--private-reference',action='append',type=Path,default=None,
                        help='Explicit private reference directory; repeat for each corpus. Defaults to existing task-local millennium/internal and literature/internal directories.')
    args=parser.parse_args();source=args.source.absolute()
    require(source.is_dir() and not any(p.is_symlink() for p in [source,*source.parents]),'Real source directory required')
    require(not args.output.resolve().is_relative_to(source.resolve()),'Scope receipt must remain outside inspected source')
    forbidden=set();reference_files=0
    reference_directories = args.private_reference
    if reference_directories is None:
        reference_directories = [p for p in [task/'millennium/internal',task/'literature/internal'] if p.exists()]
    for directory in reference_directories:
        require(directory.is_dir() and not directory.is_symlink(),'Explicit private/internal input directory is missing or not real')
        for path in directory.rglob('*'):
            require(not path.is_symlink(),'Symlink in private reference inventory')
            if path.is_file():
                raw=path.read_bytes()
                if raw:forbidden.add(hashlib.sha256(raw).hexdigest())
                reference_files+=1
            else:require(path.is_dir(),'Special private reference file')
    public=[]
    for path in sorted(source.rglob('*')):
        require(not path.is_symlink(),'Symlink in publication source')
        if path.is_dir():continue
        require(path.is_file(),'Special file in publication source')
        rel=path.relative_to(source).as_posix();parts={p.lower() for p in path.relative_to(source).parts}
        require(not (parts & {'internal','private','_internal','_private','__pycache__'}),'Private/internal/cache path in publication source: '+rel)
        require(not rel.startswith('millennium/') or rel=='millennium/HDBLAST_METHODS_TRANSFER_PUBLIC_SAFE.md',
                'Only the public-safe methods note may be published from the millennium audit')
        raw=path.read_bytes();pin=hashlib.sha256(raw).hexdigest()
        require(pin not in forbidden,'An exact private audit or raw-paper body was copied into publication source')
        public.append({'path':rel,'bytes':len(raw),'sha256':pin})
    require(public,'Empty source scope')
    receipt={'status':'PASS','source_files_inspected':len(public),'private_reference_files_checked_locally':reference_files,
             'private_reference_directories_checked':len(reference_directories),
             'private_or_internal_paths_published':0,'exact_private_input_bodies_published':0,
             'scope':'Path and exact-byte exclusion audit. Private input paths, hashes and contents are not included. This does not detect paraphrases or edited excerpts; the public-safe note requires separate review.',
             'source_files':public}
    args.output.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='source_files'},indent=2))


if __name__=='__main__':main()
