"""Construct compact laboratory mass protocols from paired proper-time wall samples."""
from bisect import bisect_right
import hashlib
import io
import json
import csv
from math import exp
from pathlib import Path
import sys
import zipfile

PIN='9a7069e44c379e68555511e82fe841e0f01a72b23e77fa5ec69caa1373beaf70'
RHO=7.838285538073243
INDICES=range(42,59)


def samples(archive, name):
    with archive.open(name) as f:
        rows=list(csv.DictReader(io.TextIOWrapper(f,encoding='ascii')))
    s=[0.]
    for past,now in zip(rows,rows[1:]):
        dt=float(now['T'])-float(past['T'])
        s.append(s[-1]+dt*(exp(float(past['B_wall']))+exp(float(now['B_wall'])))/(2*RHO))
    if not all(x<y for x,y in zip(s,s[1:])) or s[-1]<.15:
        raise ValueError('Invalid proper-time history')
    phi=[float(row['phi_wall']) for row in rows]
    values=[]
    for i in INDICES:
        t=i/400
        j=bisect_right(s,t)-1
        values.append(phi[j]+(phi[j+1]-phi[j])*(t-s[j])/(s[j+1]-s[j]))
    return phi[0],values


def main():
    raw=Path(sys.argv[1]).read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PIN:
        raise ValueError('Expected the published wall dataset')
    out={'input_sha256':PIN,'h':'0.0025','centers_in_h':list(INDICES),'coefficient_unit':'H0 squared','protocols':{}}
    with zipfile.ZipFile(io.BytesIO(raw)) as archive:
        for g in ('g4','g8'):
            for amp,folder,plus,minus in [('full','runs','delayed','delayed_negative'),
                                         ('half','half_runs','delayed_positive','delayed_negative')]:
                prefix=f'{folder}/{g}/'
                phi0,q=samples(archive,prefix+'quiet.csv')
                for sign,source in [('positive',plus),('negative',minus)]:
                    initial,p=samples(archive,prefix+source+'.csv')
                    if initial!=phi0:
                        raise ValueError('Wall values differ at the initial cut')
                    coeff=[format(10000*(a-b)*(a+b-1),'.12g') for a,b in zip(p,q)]
                    out['protocols'][f'{g}_{amp}_{sign}']=coeff
                out['phi0']=format(phi0,'.16g')
    print(json.dumps(out,indent=2))


if __name__=='__main__':
    main()
