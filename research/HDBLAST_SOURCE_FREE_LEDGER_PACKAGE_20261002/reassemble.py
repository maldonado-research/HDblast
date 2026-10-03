#!/usr/bin/env python3
"""Reassemble the complete standalone ZIP without evaluating physical data."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import zipfile

def require(ok, message):
    if not ok: raise ValueError(message)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--parts-manifest',type=Path,required=True)
    parser.add_argument('--expected-sha256',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(re.fullmatch('[0-9a-f]{64}',args.expected_sha256) is not None,'Explicit expected archive SHA required')
    manifest=json.loads(args.parts_manifest.read_text())
    require(manifest['schema_version']==1 and manifest['zip_sha256']==args.expected_sha256,'Wrong archive manifest')
    require(not args.output.exists() and not args.output.is_symlink(),'Fresh output path required')
    require(len({p['path'] for p in manifest['parts']})==len(manifest['parts']),'Duplicate part')
    total=0; archive_hash=hashlib.sha256()
    with args.output.open('xb') as output:
        for entry in manifest['parts']:
            name=entry['path']
            require(Path(name).name==name and name not in ('.','..') and '\\' not in name,'Unsafe part name')
            part=args.parts_manifest.parent/name
            require(part.is_file() and not part.is_symlink(),'Missing or symlink part')
            require(part.stat().st_size==entry['bytes'] and 0<entry['bytes']<=24*1024**2,'Part size differs')
            digest=hashlib.sha256()
            with part.open('rb') as source:
                for block in iter(lambda:source.read(1024**2),b''):
                    output.write(block);digest.update(block);archive_hash.update(block);total+=len(block)
            require(digest.hexdigest()==entry['sha256'],'Part SHA differs: '+name)
    require(total==manifest['zip_bytes'] and archive_hash.hexdigest()==args.expected_sha256,'Complete archive differs')
    with zipfile.ZipFile(args.output) as archive:
        require(archive.testzip() is None,'ZIP CRC failed')
        path=manifest['checkpoint']+'/MANIFEST.json'
        require(hashlib.sha256(archive.read(path)).hexdigest()==manifest['package_manifest_sha256'],
                'Embedded payload manifest differs')
    print(json.dumps({'status':'PASS_REASSEMBLED_PACKAGE','zip_bytes':total,'zip_sha256':archive_hash.hexdigest(),
                      'physical_evaluations':0}))

if __name__=='__main__':main()
