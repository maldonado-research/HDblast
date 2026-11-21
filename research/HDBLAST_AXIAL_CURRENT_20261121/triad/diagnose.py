"""Charge-blind angular and pressure tests for three orthogonal flavor beams."""
import argparse
import json
from pathlib import Path
from itertools import combinations_with_replacement

import numpy as np
from scipy.optimize import brentq
from scipy.special import i0e


# Fixed angular selections measured in the local orthonormal frame of the shell.
CAP_DEGREES = (2., 3., 5., 8., 10., 15., 20., 25., 30.)
# Ordered upper bounds on the six direction cosines of the boost against the beam poles.
POLE_COSINE_UPPER = np.array([1., 2.**-.5, 3.**-.5, 0., 0., -3.**-.5])


def cap_ceiling(angle):
    """Thermal probability ceilings for retaining the faintest 5, 4 or 3 poles."""
    area = .5*(1.-np.cos(angle))
    cosines = np.cos(np.maximum(0., np.arccos(POLE_COSINE_UPPER[1:])-angle))
    # At fixed |p|, dP/dOmega = k exp(k mu)/(4 pi sinh k).
    # Jensen on [2v-1,1] gives k exp(k v)/sinh k <= 1/(1-v) for v>=0.
    density = 1./(1.-np.maximum(0., cosines))
    return np.array([area*density[i:].sum() for i in range(3)])


def fourth_rank_test(n, weights, cubic, cyclic=True):
    """Distance of the measured l=4 harmonic vector from any axisymmetric orbit."""
    # Integrating products of two l=4 harmonics requires degree 8 in z and eight
    # azimuth harmonics. Five Gauss nodes and nine equally spaced azimuths suffice.
    z, wz = np.polynomial.legendre.leggauss(5)
    az = 2*np.pi*np.arange(9)/9
    radial = np.sqrt(1.-z*z)[:,None]
    axes = np.column_stack(((radial*np.cos(az)).ravel(), (radial*np.sin(az)).ravel(), np.repeat(z,9)))
    area = np.repeat(wz/18.,9)
    powers = list(combinations_with_replacement(range(3), 2))
    basis = np.column_stack([n[:,i]*n[:,j] for i,j in powers])
    second = basis.T@weights
    fourth = (basis.T*weights)@basis

    def legendre_moment(u):
        coeff = np.column_stack([(1+int(i!=j))*u[:,i]*u[:,j] for i,j in powers])
        mu2 = coeff@second
        mu4 = ((coeff@fourth)*coeff).sum(axis=1)
        return (35*mu4-30*mu2+3)/8

    f = (sum(legendre_moment(np.roll(axes,j,axis=1)) for j in range(3))/3.
         if cyclic else legendre_moment(axes))
    harmonic = 2.5*((axes**4).sum(axis=1)-.6)
    if abs(area @ (harmonic*(f-(7./12)*cubic*harmonic)))>1.e-11:
        raise ArithmeticError("Cubic harmonic projection")
    residual = 3*np.sqrt(area@(f-(7./12)*cubic*harmonic)**2)
    coefficient_upper = (7./12)*abs(cubic)+residual
    l4_norm_squared = 9*(area@(f*f))
    gap = np.sqrt(max(0., l4_norm_squared-coefficient_upper**2))/np.sqrt(20./7)
    return dict(natural_cyclic_axial_gap=float(gap), l4_norm_squared=float(l4_norm_squared),
                noncubic_supremum_upper=float(residual))


def angular(n, weights, full_moments=True):
    """One flavor includes both charge signs. Every cyclic copy has weight 1/3."""
    weights = weights/weights.sum()
    z = n[:,2]
    transverse = n[:,0]+1j*n[:,1]
    axial = np.dot(weights,z**4+.75*abs(transverse)**4)
    spin4 = np.dot(weights,transverse**4)
    cubic_mean = 2.5*(axial-.6)
    cubic_amplitude = .625*abs(spin4)
    natural = cubic_mean+.625*spin4.real
    low = cubic_mean-cubic_amplitude
    high = cubic_mean+cubic_amplitude
    runs = []
    for angle_deg in CAP_DEGREES:
        a = np.radians(angle_deg)
        c = np.cos(a)
        plus = weights[z >= c].sum()
        minus = weights[z <= -c].sum()
        # Each sign of each axis contains at least its own beam's contribution.
        # Dropping the brightest pole leaves both signs of two axes and one
        # sign of the last; dropping two or three leaves at least this mass.
        observed = np.array([(2*(plus+minus)+min(plus,minus))/3.,
                             ((plus+minus)+2*min(plus,minus))/3., min(plus,minus)])
        thermal = cap_ceiling(a)
        gap = observed-thermal
        i = int(np.argmax(gap))
        # Axisymmetric distributions may put arbitrary mass on any latitude.
        # Keep both caps of the two axes farthest from the symmetry axis.
        ring_ceiling = 4*np.arcsin(np.sqrt(2)*np.sin(a))/np.pi
        ring_floor = (2./3.)*(plus+minus)
        runs.append(dict(axisymmetric_gap=max(0.,float(ring_floor-ring_ceiling)),
                         axisymmetric_observed_floor=float(ring_floor), axisymmetric_target_ceiling=float(ring_ceiling),
                         angle_degrees=angle_deg, own_axis_forward=float(plus),
                         own_axis_backward=float(minus), retained_poles=5-i,
                         observed_floor=float(observed[i]), target_ceiling=float(thermal[i]),
                         any_boost_gap=max(0.,float(gap[i])),
                         stationary_gap=max(0.,float(plus+minus-3*(1-c)))))
    best = max(runs,key=lambda r:r["any_boost_gap"])
    moments=fourth_rank_test(n,weights,natural) if full_moments else {}
    return dict(**moments, cubic_harmonic_at_zero_roll=float(natural),
                cubic_harmonic_all_rolls=[float(low),float(high)],
                any_boost_gap=best["any_boost_gap"], best_angle=best["angle_degrees"],
                axisymmetric_gap=max(r["axisymmetric_gap"] for r in runs),
                stationary_gap=max(r["stationary_gap"] for r in runs), caps=runs)


def cyclic_stress(stress, roll):
    """Cyclic copies of a locally rolled (T00-normalized) energy-momentum tensor."""
    c, s = np.cos(roll), np.sin(roll)
    r = np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
    space = r @ stress[1:,1:] @ r.T
    # Cyclic averaging gives equal diagonal/off-diagonal components and axial flux.
    diag = np.trace(space)/3.
    off = (space[0,1]+space[1,2]+space[2,0])/3.
    perp, along = diag-off, diag+2*off
    flux = np.sum(r @ stress[0,1:])/np.sqrt(3.)
    discr = np.sqrt((1.+along)**2-4*flux**2)
    rho = .5*(1.-along+discr)
    parallel = .5*(along-1.+discr)
    velocity = 2*flux/(1.+along+discr)
    return perp/rho, parallel/rho, abs(velocity)


def stress_tests(p, gamma):
    four = np.column_stack((gamma,p))
    stress = (four.T/gamma)@four/gamma.sum()
    def anisotropy(angle):
        perp, along, _ = cyclic_stress(stress,angle)
        return along-perp
    knots = np.linspace(0.,2*np.pi,65)
    values = [anisotropy(a) for a in knots]
    root = None
    for i in range(len(knots)-1):
        if values[i]*values[i+1]<=0.:
            root = brentq(anisotropy,knots[i],knots[i+1],xtol=1.e-14)
            break
    def stats(angle):
        perp, along, velocity = cyclic_stress(stress,angle)
        return dict(roll_degrees=float(np.degrees(angle)), velocity=float(velocity),
                    principal_pressures_over_rest_energy=sorted(map(float,(perp,perp,along))))
    return dict(natural=stats(0.), fixed_135_degrees=stats(3*np.pi/4.), isotropic_rest_example=stats(root) if root is not None else None)


def controls():
    # Integrate the conditional cap probability directly, without the Jensen bound.
    nodes, weights = np.polynomial.legendre.leggauss(40)
    rng = np.random.default_rng(48218)
    directions = np.vstack([np.eye(3),np.ones((1,3)),np.array([[1.,1.,0.]]),rng.normal(size=(16,3))])
    directions /= np.linalg.norm(directions,axis=1)[:,None]
    max_error=0.
    checks=0
    for deg in CAP_DEGREES:
        a = np.radians(deg)
        c = np.cos(a)
        x = .5*(1-c)*nodes+.5*(1+c)
        for direction in directions:
            cosines=np.sort(np.r_[direction,-direction])[::-1, None][1:]
            slope=np.sqrt(1-cosines**2)*np.sqrt(1-x**2)
            for k in (0., .1, 1., 5., 20., 300.):
                if k==0.:
                    cap=np.full(5,.5*(1-c))
                else:
                    values=np.exp(k*(cosines*x+slope-1))*i0e(k*slope)
                    cap=.5*(1-c)*k/(-np.expm1(-2*k))*(values@weights)
                error=np.max(np.cumsum(cap[::-1])[::-1][:3]-cap_ceiling(a))
                max_error=max(max_error,float(error))
                checks+=3
    # Independent latitude-ring intersection calculation for the four retained caps.
    ring_tests=0
    for deg in CAP_DEGREES:
        a=np.radians(deg)
        for direction in directions:
            d=np.delete(direction,np.argmax(abs(direction)))
            d=np.r_[d,-d]
            for polar in np.linspace(.01,np.pi-.01,61):
                edge=(np.cos(a)-d*np.cos(polar))/(np.sqrt(1.-d*d)*np.sin(polar))
                share=np.arccos(np.clip(edge,-1.,1.)).sum()/np.pi
                if share>4*np.arcsin(np.sqrt(2)*np.sin(a))/np.pi+1e-14:
                    raise ValueError("Latitude-ring bound")
                ring_tests+=1
    # Positive-and-negative pencil beams have H4=1; spherical cubature gives H4=0.
    direction=np.vstack([np.eye(3),-np.eye(3)])
    if abs(np.mean(np.sum(direction**4,axis=1)) - 1.) > 1e-14:
        raise ValueError("Pencil beam control")
    z,w=np.polynomial.legendre.leggauss(12)
    az=np.arange(16)*2*np.pi/16.
    rr=np.sqrt(1-z*z)[:,None]
    grid=np.column_stack([(rr*np.cos(az)).ravel(),(rr*np.sin(az)).ravel(),np.repeat(z,16)])
    weight=np.repeat(w/32.,16)
    if abs(np.dot(weight,2.5*(np.sum(grid**4,axis=1)-.6)))>1e-13 or max_error>2e-12:
        raise ValueError("Angular integration controls")
    # A pure axial pole around (1,1,1) is unchanged by the cyclic completion.
    # The six Cartesian poles instead have a pure cubic harmonic.
    axis=np.ones((1,3))/np.sqrt(3.)
    axial_test=fourth_rank_test(axis,np.ones(1),-2./3.)
    cubic_test=fourth_rank_test(direction,np.ones(6)/6.,1.)
    cubic_exact=np.sqrt((7./12-(7./12)**2)*7./20)
    if axial_test["natural_cyclic_axial_gap"]>1.e-7 or abs(cubic_test["natural_cyclic_axial_gap"]-cubic_exact)>1.e-13:
        raise ArithmeticError("Axial and cubic harmonic controls")
    # Compute stress eigenvalues two ways on a tilted test population, with
    # nonzero charge-blind flux so a rest-frame transformation is required.
    sample=rng.normal(size=(700,3))*[.12,2.,12.]+[.1,.2,3.]
    gamma=np.hypot(1.,np.linalg.norm(sample,axis=1))
    four=np.column_stack((gamma,sample))
    stress=(four.T/gamma)@four/gamma.sum()
    maximum_stress_error=0.
    for roll in np.linspace(.1,5.,8):
        c,ss=np.cos(roll),np.sin(roll)
        rotated=np.column_stack((c*sample[:,0]-ss*sample[:,1],ss*sample[:,0]+c*sample[:,1],sample[:,2]))
        merged=np.vstack((rotated,rotated[:,[1,2,0]],rotated[:,[2,0,1]]))
        gs=np.hypot(1.,np.linalg.norm(merged,axis=1))
        vec=np.column_stack((gs,merged))
        tensor=(vec.T/gs)@vec/gs.sum()
        eig=np.sort(np.linalg.eigvals(tensor*np.array([-1.,1.,1.,1.])))
        if np.max(abs(eig.imag))>1.e-12: raise ValueError("Nonreal full stress eigenvalues")
        pressure=np.sort(eig.real[1:]/-eig.real[0])
        per,parallel,_=cyclic_stress(stress,roll)
        maximum_stress_error=max(maximum_stress_error,float(np.max(abs(pressure-np.sort([per,per,parallel])))))
    if maximum_stress_error>1.e-11: raise ValueError("Explicit three-flavor stress construction failed")
    return dict(conditional_cap_tests=checks,latitude_checks=ring_tests,max_cap_excess=max_error,
                full_tensor_error=maximum_stress_error, six_pole_axial_floor=cubic_test["natural_cyclic_axial_gap"])


def evaluate(path):
    with np.load(path) as data:
        p=data["p"]
    if p.ndim!=2 or p.shape[1]!=3 or not np.all(np.isfinite(p)):
        raise ValueError("Expected equal-weight wall-frame three-momenta")
    modulus=np.linalg.norm(p,axis=1)
    gamma=np.hypot(1.,modulus)
    if np.min(modulus)<=0:
        raise ValueError("Static particle with undefined direction")
    n=p/modulus[:,None]
    quantiles=np.quantile(gamma,[.25,.5,.75])
    if not np.all(np.diff(quantiles)>0):
        raise ValueError("Degenerate sample energy quartiles")
    bands=np.searchsorted(quantiles,gamma)
    groups=[]
    for i in range(4):
        picked=(bands==i)
        s=angular(n[picked],np.ones(np.count_nonzero(picked)),full_moments=False)
        groups.append(dict(count=int(np.sum(picked)), any_boost_gap=s["any_boost_gap"],
                           axisymmetric_gap=s["axisymmetric_gap"], best_angle=s["best_angle"],
                           cubic_harmonic_all_rolls=s["cubic_harmonic_all_rolls"]))
    return dict(snapshot=path.stem,number=int(len(p)),number_weighted=angular(n,np.ones(len(n))),
                energy_weighted=angular(n,gamma), energy_edges=list(map(float,quantiles)),
                by_energy_quartile=groups,three_flavor_stress=stress_tests(p,gamma))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--snapshots",type=Path,required=True)
    parser.add_argument("--shared-snapshots",type=Path)
    parser.add_argument("--shared-later-snapshots",type=Path)
    parser.add_argument("--cross-snapshots",type=Path)
    parser.add_argument("--return-snapshots",type=Path)
    parser.add_argument("--output",type=Path,required=True)
    opt=parser.parse_args()
    files=sorted(opt.snapshots.glob("*.npz"))
    if not files: raise ValueError("No snapshot arrays")
    result=dict(controls=controls(),runs=[evaluate(p) for p in files])
    if opt.cross_snapshots:
        later=sorted(opt.cross_snapshots.glob("*.npz"))
        if not later: raise ValueError("No later crossing snapshot arrays")
        result["later_runs"]=[evaluate(p) for p in later]
    if opt.return_snapshots:
        returned=sorted(opt.return_snapshots.glob("*.npz"))
        if not returned: raise ValueError("No return-field snapshot arrays")
        result["return_runs"]=[evaluate(p) for p in returned]
    shared_files=[]
    for directory in (opt.shared_snapshots,opt.shared_later_snapshots):
        if directory:
            found=list(directory.glob("*/*.npz"))
            if not found: raise ValueError("No shared-flavor snapshot arrays")
            shared_files.extend(found)
    keys=set()
    for p in sorted(shared_files,key=lambda path:(path.parent.name,path.name)):
        key=(p.parent.name,p.name)
        if key in keys: raise ValueError("Duplicate shared-flavor snapshot")
        keys.add(key)
        group="shared_later_runs" if p.name.startswith("t7.") else "shared_runs"
        result.setdefault(group,[]).append(dict(flavor_mode=p.parent.name,**evaluate(p)))
    opt.output.parent.mkdir(parents=True,exist_ok=True)
    opt.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+"\n")
