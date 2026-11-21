"""Analytic checks on the canonical fold strip and its numerical enclosure."""
import numpy as np
from memory_bound import angular, measure, paired_bound


def steady_trace(n, W, a):
    # K=d=1, axial travel W, target a; exact final canonical occupancy is
    # (a/u) times the multiplicity of the two translating birth strips at |u|.
    t=np.linspace(0,W,n+1)
    f=np.linspace(0,1,n+1)
    ones=np.ones(n)
    return np.column_stack((t[:-1],np.diff(t),ones,ones,ones,ones,np.diff(f),
                            f[:-1],f[1:],0*ones,ones,a*ones))


def main():
    for n in (1,32,97,1024):
        W, a = 16., .20
        data=steady_trace(n,W,a)
        res=measure(data,1.)
        # The exact full excess follows by integrating the two axial intervals
        # of length W whose overlapping length at |u| is W-2|u|.
        exact = a-a*a/W
        assert .95*exact < res['combined_lower_fraction'] < exact
        assert np.isclose(res['angular_lower_fraction'],a/2,rtol=1e-12)
        # Independent direct optimization of the single rectangular fold strip.
        ell=np.linspace(0,2*a,100001)
        rect=np.max((W-2*ell)*(12*a*ell-3*ell**2))/(12*a*W)
        if n in (1,32,1024):
            assert abs(res['fold_only_lower_fraction']-rect)<1e-10
        print(f"steady grid {n}: exact {exact:.10f}, lower {res['combined_lower_fraction']:.10f}")
    # Reversal d/T from -2 to 2 and c=0. Compare quadrature bounds to the
    # independent closed integral of the triangular angular deficit.
    one=np.ones(1)
    exact=(.75+.5*np.log(2))/2
    coarse=angular(0*one,-2*one,2*one,one,64)[0]
    fine=angular(0*one,-2*one,2*one,one,4096)[0]
    assert coarse<fine<exact and exact-fine<.001
    # An allowed Z jump can only decrease the reported strip certificate.
    jumps=np.array([0.,.001,.01,.1,1.])
    trial=paired_bound(8.,1.,1.,1.,jumps,2.,True)
    assert np.all(np.diff(trial)<=0) and trial[-1]==0.
    print("reversal and jump checks passed")


if __name__=='__main__':
    main()
