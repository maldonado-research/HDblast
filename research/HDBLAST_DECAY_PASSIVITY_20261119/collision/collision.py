"""Exponentiate a squeezed-oscillator collision in the physical two-mode Fock basis."""
import argparse
import json
import math

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import expm_multiply


CASES = ("physical_vacuum", "opposite_phase", "target_vacuum")
TARGET_N = 0.25


def projected_vacuum(dim, squeeze):
    """Normalized projection of S(r)|0> in the physical oscillator number basis."""
    state = np.zeros(dim)
    state[0] = 1
    for j in range(2, dim, 2):
        state[j] = -math.tanh(squeeze) * math.sqrt((j-1)/j) * state[j-2]
    return state / np.linalg.norm(state)


def collision(dim, eta, check_work):
    r = math.asinh(math.sqrt(TARGET_N))
    theta = math.acos(math.sqrt(eta))
    a = sparse.diags(np.sqrt(np.arange(1, dim)), 1, shape=(dim, dim), format="csr")
    ident = sparse.eye(dim,format="csr")
    b = math.cosh(r)*a + math.sinh(r)*a.T
    # In tensor order (physical system A, physical ancilla C), U=exp(theta G).
    exchange = sparse.kron(b, a.T, format="csr")
    generator = exchange - exchange.T
    physical_zero = projected_vacuum(dim, 0.)
    squeezes = (0., 2*r, r)
    initial = np.column_stack([np.kron(projected_vacuum(dim,v),physical_zero) for v in squeezes])
    final = expm_multiply(theta*generator,initial,traceA=0)
    a_numbers = np.repeat(np.arange(dim),dim)
    c_numbers = np.tile(np.arange(dim),dim)
    def means(psi):
        probs = abs(psi)**2
        n_a = a_numbers @ probs
        n_c = c_numbers @ probs
        n_b = np.sum(abs(sparse.kron(b,ident,format="csr") @ psi)**2, axis=0)
        return n_a, n_c, n_b
    na0,nc0,nb0 = means(initial)
    naf,ncf,nbf = means(final)
    # The dimensionless endpoint Hamiltonians are N_A+1/2 and N_C+1/2.
    # Integrate the external work derivative independently on the finest truncation.
    work_integral = None
    if check_work:
        pair = sparse.kron(a,a,format="csr")
        bare_commutator = 2*math.sinh(r)*(pair+pair.T)
        work_integral = np.zeros(len(CASES))
        points,weights = np.polynomial.legendre.leggauss(8)
        for x,w in zip(points,weights):
            state = expm_multiply((theta*(1+x)/2)*generator,initial,traceA=0)
            velocity = np.sum(state.conj() * (bare_commutator @ state),axis=0).real
            work_integral += w*(theta/2)*velocity
    expected_n_a = (0., math.sinh(2*r)**2, TARGET_N)
    expected_n_b = (TARGET_N, TARGET_N, 0.)
    result = {}
    for j,name in enumerate(CASES):
        system = float(naf[j]-na0[j])
        ancilla = float(ncf[j]-nc0[j])
        work = system+ancilla
        exact_system = (1-eta)*(TARGET_N-expected_n_a[j])
        exact_ancilla = (1-eta)*expected_n_b[j]
        row = dict(nA_initial=na0[j], nB_initial=nb0[j], nB_final=nbf[j],
                   nC_initial=nc0[j], nA_final=naf[j], nC_final=ncf[j],
                   delta_system=system,delta_ancilla=ancilla,work_on_pair=work,
                   max_error_against_infinite_oscillator=max(abs(system-exact_system),
                        abs(ancilla-exact_ancilla),abs(work-exact_system-exact_ancilla)),
                   norm_error=float(abs(np.vdot(final[:,j],final[:,j])-1)))
        if work_integral is not None:
            row["work_commutator_integral"] = float(work_integral[j])
            row["work_integral_minus_endpoint"] = float(work_integral[j]-work)
        result[name] = row
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--eta",type=float,default=0.64)
    parser.add_argument("--dimensions",type=int,nargs="+",default=[16,32,64,96,128])
    args = parser.parse_args()
    if not 0 < args.eta < 1 or any(d<4 for d in args.dimensions):
        parser.error("Choose 0 < eta < 1 and Fock dimensions of at least 4")
    # SciPy may sample columns when estimating the exponential action norm.
    np.random.seed(913)
    print(json.dumps({"units":"one physical oscillator quantum E", "eta":args.eta, "target_nA":TARGET_N,
                      "truncations": {d:collision(d,args.eta,d==max(args.dimensions))
                                      for d in args.dimensions}},allow_nan=False,indent=2))


if __name__ == "__main__":
    main()
