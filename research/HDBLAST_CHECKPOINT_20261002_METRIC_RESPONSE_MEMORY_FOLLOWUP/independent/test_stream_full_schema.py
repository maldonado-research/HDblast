#!/usr/bin/env python3
"""Synthetic full-size631-member NPZ streaming check; no physical arrays read."""
from __future__ import annotations
import argparse
import hashlib
import io
import json
from pathlib import Path
import resource
import zipfile
import numpy as np
from stream_npz import StreamingNpz,file_sha256


def require(condition,message):
    if not condition:
        raise RuntimeError(message)


def synthetic_array(member,index):
    dtype=np.dtype(member['dtype']);shape=tuple(member['shape'])
    rng=np.random.default_rng(index+177)
    if dtype.kind in 'US':
        value=np.full(shape,'rawtest',dtype=dtype)
    elif dtype.kind in 'iu':
        value=rng.integers(-512,512,size=shape,dtype=dtype)
    elif dtype.kind in 'fc':
        value=rng.integers(-512,512,size=shape,dtype=np.int64).astype(dtype)/dtype.type(16)
        if dtype.kind=='c':
            value=value+dtype.type(1j)*rng.integers(-512,512,size=shape,dtype=np.int64).astype(dtype)/dtype.type(32)
    else:
        raise RuntimeError('Unsupported synthetic schema dtype '+str(dtype))
    return np.asfortranarray(value) if member['fortran_order'] else value


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output_dir.exists(),'Fresh synthetic full-schema directory required')
    args.output_dir.mkdir(parents=True)
    here=Path(__file__).resolve().parent
    schema=json.loads((here/'PRESERVED_ARCHIVE_SCHEMA.json').read_text())
    require(schema['member_count']==631,'Expected prior raw member inventory changed')
    archive=args.output_dir/'synthetic_full_schema.npz'
    member_hashes={}
    with StreamingNpz(archive) as writer:
        for index,member in enumerate(schema['members']):
            value=synthetic_array(member,index)
            reference=io.BytesIO()
            np.lib.format.write_array(reference,value,allow_pickle=False)
            member_hashes[member['name']+'.npy']=hashlib.sha256(reference.getbuffer()).hexdigest()
            writer[member['name']]=value
            del value,reference
    with zipfile.ZipFile(archive) as saved:
        require(saved.namelist()==list(member_hashes),'631-member order or coverage differs')
        for member in saved.infolist():
            with saved.open(member) as stream:
                digest=hashlib.sha256()
                for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
            require(digest.hexdigest()==member_hashes[member.filename],'Exact synthetic member bytes changed')
    with np.load(archive,allow_pickle=False) as saved:
        for index,member in enumerate(schema['members']):
            value=saved[member['name']]
            expected=synthetic_array(member,index)
            require(value.dtype==expected.dtype and value.shape==expected.shape and np.array_equal(value,expected),
                    'Synthetic full-schema dtype/shape/value differs: '+member['name'])
            del value,expected
    peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    require(peak<=256*1024,'Synthetic stream alone exceeds unchanged256MiB budget')
    report={'schema_version':1,'status':'PASS','scope':'Synthetic random integer/fraction arrays only, with prior saved-header shapes/dtypes. No physical saved-array values loaded or recomputed.',
            'physical_evaluations':0,'members_exact_npy_bytes_checked':len(member_hashes),'members_dtype_shape_value_checked':len(member_hashes),
            'fine_raw_npy_schema_bytes':schema['uncompressed_npy_bytes'],'synthetic_zip_bytes':archive.stat().st_size,
            'synthetic_peak_rss_kib':peak,'synthetic_archive_sha256':file_sha256(archive),
            'archive_schema_sha256':file_sha256(here/'PRESERVED_ARCHIVE_SCHEMA.json'),'stream_module_sha256':file_sha256(here/'stream_npz.py'),
            'limit':'This validates format and archive-only memory use; actual future physical producer peak remains an execution gate.'}
    (args.output_dir/'FULL_SCHEMA_STREAM_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
