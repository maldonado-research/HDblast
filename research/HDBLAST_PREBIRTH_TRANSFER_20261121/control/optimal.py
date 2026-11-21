"""Maximum field stored at birth by waveform shaping within a given incoming envelope."""
import argparse
import itertools
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np


WAVENUMBERS = (2, 5, 10, 20, 50, 103, 258, 515, 1030)


def maximizing_vertex(v):
    """Farthest vertex of the sum of centered line segments in the E,B plane."""
    # Reorienting a centered segment leaves the reachable set unchanged.
    v = v.copy()
    v[v[:, 1] < 0.] *= -1.
    v = v[np.argsort(np.arctan2(v[:, 1], v[:, 0]), kind='stable')]
    start = -np.sum(v, axis=0)
    vertices = start + 2.*np.cumsum(v, axis=0)
    furthest = vertices[np.argmax(np.einsum('ij,ij->i', vertices, vertices))]
    return furthest if np.linalg.norm(furthest) >= np.linalg.norm(start) else start


def check_polygon():
    edges = np.array([[0.,1.],[-.31,.42],[1.,-1.3],[-.8,-.3],[.1,0.],[.18,.27]])
    all_signs = np.array(list(itertools.product((-1.,1.), repeat=len(edges))))
    brute = np.max(np.linalg.norm(all_signs @ edges, axis=1))
    if abs(brute-np.linalg.norm(maximizing_vertex(edges))) > 1.e-12:
        raise ValueError('Wrong reachable polygon')


def prebirth(raw):
    arr = np.memmap(raw, dtype='=f8', mode='r').reshape(-1,10)
    j = int(np.flatnonzero(arr[:,9] > 0.)[0])
    data = arr[:j, [0,1,3,7,4,5]].copy(order='C')
    return np.vstack(([0.,2.,data[0,2],0.,0.,0.], data))


def response(data, exe, temp, k):
    nodes = temp/'nodes'; weights = temp/'weights'
    data[:,:3].copy(order='C').tofile(nodes)
    subprocess.run([str(exe),str(k),str(nodes),str(weights)], check=True)
    kernel = np.fromfile(weights, dtype='=f8').reshape(-1,2)
    incident = .5*(data[:-1,3]+data[1:,3])
    vectors = kernel * incident[:,None]
    coherent = np.sum(vectors,axis=0)
    # A sign chosen independently on each cell also covers every weaker signal.
    extreme = maximizing_vertex(vectors)
    opt = np.linalg.norm(extreme)
    direction = extreme/opt
    support = np.sum(np.abs(vectors @ direction))
    if abs(support/opt-1.) > 2.e-8:
        raise ValueError('Inconsistent support and vertex')
    peak = np.linalg.norm(coherent)
    return dict(fixed_profile_field_EB=coherent.tolist(),
                fixed_profile_spatial_rms=float(peak/math.sqrt(2.)),
                optimal_envelope_spatial_rms=float(opt/math.sqrt(2.)),
                fixed_over_optimal=float(peak/opt))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--transfer', type=Path, required=True,help='Directory of four fresh capture .bin files, or its parent')
    p.add_argument('--output', type=Path, required=True)
    opt=p.parse_args()
    root=opt.transfer
    if (root/'raw').is_dir(): root=root/'raw'
    result={}
    check_polygon()
    with tempfile.TemporaryDirectory(prefix='prebirth-control-') as wd:
        temp=Path(wd); exe=temp/'adjoint'
        subprocess.run(['g++','-std=c++17','-O2',str(Path(__file__).with_name('adjoint.cpp')),'-o',str(exe)],check=True)
        for resolution in ('native','half'):
            data=prebirth(root/f'dedicated_{resolution}.bin')
            shared=prebirth(root/f'shared_{resolution}.bin')
            if not np.array_equal(shared,data):
                raise ValueError('Cases have different prebirth fields')
            rows={str(k):response(data,exe,temp,k) for k in WAVENUMBERS}
            k2=rows['2']['fixed_profile_field_EB']
            error=np.linalg.norm(data[-1,4:6]-k2)/np.linalg.norm(data[-1,4:6])
            if error > 2.e-5: raise ValueError(f'Lowest harmonic disagrees with recorded solution: {error}')
            result[resolution]=dict(birth_time=float(data[-1,0]), coincident_prebirth_histories=True,
                                    lowest_mode_relative_vector_error=float(error), modes=rows)
    for k in WAVENUMBERS:
        for name in ('optimal_envelope_spatial_rms','fixed_profile_spatial_rms'):
            a,b=[result[resolution]['modes'][str(k)][name] for resolution in ('native','half')]
            if abs(a-b)/max(a,b) > .035: raise ValueError(f'Time refinement changes {name} at {k}: {a} vs {b}')
    opt.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
