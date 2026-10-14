"""Per-sphere KID charge densities of the second-order constraint sources."""
import json, pickle, sys
import sympy as sp
from engine import X, O
from sphere import sphere_average, parity

case = sys.argv[1]
d = pickle.load(open(f'../results/sources_{case}.pkl', 'rb'))
SH, SM = d['SH'], d['SM']
A, B = d['A'], d['B']
out = {'case': case, 'A': str(A.tolist()), 'B': str(B.tolist())}
out['SH_parity'] = parity(SH)
out['SM_parity'] = [parity(c) for c in SM]
# boost densities: <Omega x_a S_H> (sqrt(gamma) is radial; omitted per sphere)
out['boost_density_zero'] = [sp.simplify(sphere_average(O[0] * X[a] * SH)) == 0 for a in range(4)]
rot = {}
for i in range(4):
    for j in range(i + 1, 4):
        Y = [0] * 4
        Y[i], Y[j] = -X[j], X[i]          # Y = x_i d_j - x_j d_i  (components)
        val = sp.factor(sphere_average(sum(Y[k] * SM[k] for k in range(4))))
        rot[f'J{i+1}{j+1}'] = str(val)
out['rotation_densities'] = rot
f0 = sp.factor(sphere_average(SH))
out['l0_hamiltonian_density'] = str(f0)
out['l0_momentum_radial_density'] = str(sp.factor(sphere_average(sum(X[k] * SM[k] for k in range(4)))))
json.dump(out, open(f'../results/charges_{case}.json', 'w'), indent=2)
print(json.dumps({k: v for k, v in out.items() if k != 'l0_hamiltonian_density'}, indent=2))
print('l0 Hamiltonian density (len):', len(str(f0)))
