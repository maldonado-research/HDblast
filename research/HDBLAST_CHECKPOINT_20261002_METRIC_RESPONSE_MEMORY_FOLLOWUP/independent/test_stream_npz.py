#!/usr/bin/env python3
"""Synthetic NPZ format/lifetime controls; never call a physical function."""
from __future__ import annotations
import argparse
import gc
import hashlib
import io
import json
from pathlib import Path
import resource
import weakref
import zipfile

import numpy as np
from stream_npz import StreamingNpz,file_sha256


def require(condition,message):
    if not condition:
        raise RuntimeError(message)


def npy_bytes(array):
    stream=io.BytesIO()
    np.lib.format.write_array(stream,array,allow_pickle=False)
    return stream.getvalue()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    require(not args.output_dir.exists(),'Fresh synthetic output directory required')
    args.output_dir.mkdir(parents=True)
    checks=[]
    arrays={'real':np.arange(21,dtype=np.longdouble).reshape(3,7),
            'complex':np.arange(12,dtype=np.clongdouble)*(np.clongdouble(1)+np.clongdouble(2j)),
            'integer':np.array([1,-4,19],dtype=np.int64),'labels':np.array(['q','rho0','current']),
            'scalar':np.array('0.125',dtype=np.longdouble),'empty':np.empty((0,5),dtype=np.longdouble),
            'fortran':np.asfortranarray(np.arange(20,dtype=np.longdouble).reshape(4,5)),
            'strided':np.arange(30,dtype=np.clongdouble)[::3]}
    increment=np.ldexp(np.longdouble(1),-60)
    precision=np.array([np.longdouble(1)+increment,np.longdouble(1)-increment],dtype=np.longdouble)
    require(np.any(precision != precision.astype(np.float64).astype(np.longdouble)),
            'Synthetic precision witness must exceed binary64 precision')
    arrays['extended_precision']=precision
    arrays['extended_complex_precision']=precision.astype(np.clongdouble)+np.clongdouble(1j)*precision[::-1]
    checks.append('synthetic extended precision witness differs from binary64')
    arrays['real'].flags.writeable=False
    regular=args.output_dir/'synthetic_numpy_reference.npz'
    streamed=args.output_dir/'synthetic_streamed.npz'
    np.savez_compressed(regular,**arrays)
    with StreamingNpz(streamed) as writer:
        writer.update(arrays)
    with zipfile.ZipFile(regular) as first,zipfile.ZipFile(streamed) as second:
        require(first.namelist()==second.namelist(),'NPZ member order changed')
        for name in first.namelist():
            require(first.read(name)==second.read(name),'Exact NPY member bytes changed: '+name)
            checks.append('exact synthetic NPY bytes '+name)
    with np.load(streamed,allow_pickle=False) as saved:
        require(saved.files==list(arrays),'numpy.load key order differs')
        for name,value in arrays.items():
            restored=saved[name]
            require(restored.dtype==value.dtype and restored.shape==value.shape and
                    np.array_equal(restored,value),'dtype/shape/value roundtrip differs: '+name)
            checks.append('numpy.load synthetic dtype/shape/value '+name)
    checks.append('exact bounded SHA equality')
    require(file_sha256(regular,37)==hashlib.sha256(regular.read_bytes()).hexdigest(),'Chunked SHA differs')
    lifetime=args.output_dir/'synthetic_lifetime.npz'
    with StreamingNpz(lifetime) as writer:
        value=np.arange(16384,dtype=np.clongdouble)
        reference=weakref.ref(value)
        writer['temporary']=value
        del value
        gc.collect()
        require(reference() is None,'Streaming writer retained an ndarray reference')
        try:writer['temporary']=np.ones(1)
        except ValueError:checks.append('duplicate member rejected')
        else:raise RuntimeError('Duplicate NPZ member accepted')
        try:writer['object']=np.array([object()],dtype=object)
        except ValueError:checks.append('object dtype rejected')
        else:raise RuntimeError('Object dtype accepted')
        try:writer['../path']=np.ones(1)
        except ValueError:checks.append('path member rejected')
        else:raise RuntimeError('Path NPZ member accepted')
    checks.append('writer releases exact synthetic array reference')
    interrupted=args.output_dir/'synthetic_interrupted.npz'
    try:
        with StreamingNpz(interrupted) as writer:
            writer['completed']=arrays['complex']
            raise RuntimeError('synthetic interruption after completed member')
    except RuntimeError as exc:
        require(str(exc)=='synthetic interruption after completed member','Unexpected interruption failure')
    with np.load(interrupted,allow_pickle=False) as saved:
        require(saved.files==['completed'] and np.array_equal(saved['completed'],arrays['complex']),
                'Completed member lost after interrupted stream')
    checks.append('completed synthetic member survives interruption')
    report={'schema_version':1,'status':'PASS','scope':'Synthetic arrays only; no source, mode, stress, quadrature or physical calls.',
            'physical_evaluations':0,'checks':checks,'check_count':len(checks),
            'peak_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'runtime':{'numpy':np.__version__,'longdouble_nmant':np.finfo(np.longdouble).nmant},
            'stream_module_sha256':file_sha256(Path(__file__).with_name('stream_npz.py'))}
    (args.output_dir/'STREAM_FORMAT_CHECKS.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    main()
