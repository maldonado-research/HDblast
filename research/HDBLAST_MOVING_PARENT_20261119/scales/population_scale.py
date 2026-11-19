"""Evaluate the trial population's initial stress and explicit junction energy terms."""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import zipfile

WALL = 'research/HDBLAST_SOURCE_WORK_20261118/coupled/WALLS.zip'
SUMMARY = 'research/HDBLAST_PULSE_PARTICLE_RESPONSE_20261118/classical/runs/fine/quiet/rolling_summary.json'
WALL_SHA = '7f3289631a80291f6180e51f5b648a135c6e33c6366d4f97c32a4b90630cab15'
SUMMARY_SHA = 'abec0887ee4de00ba0757e7d96987610f36b7157a85bc4f00696fee8064a6810'


def load(repo, name, digest):
    result = (repo/name).read_bytes()
    if hashlib.sha256(result).hexdigest() != digest:
        raise ValueError(f'Changed input {name}')
    return result


def calculate(repo):
    summary = json.loads(load(repo,SUMMARY,SUMMARY_SHA))
    reference = summary['matter']['events'][0]['A_star']
    with zipfile.ZipFile(io.BytesIO(load(repo,WALL,WALL_SHA))) as walls:
        with walls.open('fine/y1_b1/pre.csv') as data:
            row = next(csv.DictReader(io.TextIOWrapper(data)))
        with walls.open('fine/y1_b1/quiet.csv') as data:
            recorded = list(csv.DictReader(io.TextIOWrapper(data)))
    phi, h = float(row['phi_b']), float(row['H_over_H0'])
    a = math.exp(float(row['A_b'])-reference)
    g, kmax = 100., 8.
    mass = abs(g*(phi-.5))
    p = kmax/a
    e = math.hypot(mass,p)
    at = math.asinh(p/mass)
    n = p**3/(6*math.pi**2)
    # One real field, homogeneous diagonal excess occupation one on the band.
    rho = (p*e*(2*p**2+mass**2)-mass**4*at)/(16*math.pi**2)
    pressure = (p*e*(2*p**2-3*mass**2)+3*mass**4*at)/(48*math.pi**2)
    force = g**2*(phi-.5)*(p*e-mass**2*at)/(4*math.pi**2)
    weights = (rho,pressure,force)
    panels = 8192
    vals = [[],[],[]]
    for i in range(panels):
        k = p*(i+.5)/panels
        freq = math.hypot(mass,k)
        measure = p*k**2/(panels*2*math.pi**2)
        for group,val in zip(vals,(freq,k**2/(3*freq),g**2*(phi-.5)/freq)):
            group.append(measure*val)
    if any(not math.isclose(math.fsum(vec),value,rel_tol=2e-8) for value,vec in zip(weights,vals)):
        raise ValueError('Mode integrals')
    if any(float(row[k])!=0 for k in ('X','R','rad','chi','j_hat','rho_hat')) or summary['params']['Y']!=0:
        raise ValueError('Original slice has nonzero subtracted matter or friction')
    rb = summary['rho_b']
    delta,c,d = (summary['params'][k] for k in ('delta','c','d'))
    def tension(phi):
        return 2*(1-phi+phi**3/3)+delta*(1+c*phi+d*phi**2/2)
    sigma = tension(phi)
    sigma1 = 2*(phi**2-1)+delta*(c+d*phi)
    # Check the stress normalization against the actual outgoing records.
    for sample in recorded[::29]:
        r, x = float(sample['R']),float(sample['X'])
        sg = tension(float(sample['phi_b']))
        rad = rb**2*r*(2*sg+r)/36
        chi = rb**2*x*(2*sg+x+2*r)/36
        if not (math.isclose(chi,float(sample['chi']),rel_tol=3e-13,abs_tol=1e-15)
                and math.isclose(rad,float(sample['rad']),rel_tol=3e-13,abs_tol=1e-15)):
            raise ValueError('Archived expansion terms')
    scale = 6.49e-6/rb
    x, j = scale*rho, scale*force
    rho_linear, rho_squared = rb**2*sigma*x/18, rb**2*x**2/36
    force_linear, force_squared = -rb**2*sigma1*j/24, -rb**2*j**2/48
    return dict(cutoff=kmax,scale_factor=a,mass=mass,comoving_number=a**3*n,number=n,
                energy=rho,pressure=pressure,w=pressure/rho,H2=h**2,scalar_force=force,X=x,J=j,
                sigma=sigma,sigma1=sigma1,rb=rb,
                terms_in_H0sq=dict(matter_linear=rho_linear, matter_quadratic=rho_squared,
                    force_linear=force_linear, force_quadratic=force_squared),
                fractional_shifts=dict(matter=(rho_linear+rho_squared)/h**2,
                    force=(force_linear+force_squared)/h**2,
                    sum=(rho_linear+rho_squared+force_linear+force_squared)/h**2),
                normal_derivative_changes=dict(A=x/6,B=-(2*x+3*scale*pressure)/6,phi=-j/2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('repository',type=Path)
    print(json.dumps(calculate(parser.parse_args().repository),indent=2,allow_nan=False))
