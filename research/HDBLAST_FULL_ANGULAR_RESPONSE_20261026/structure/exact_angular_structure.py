#!/usr/bin/env python3
"""New bounded audit of authenticated angular data; no old producer is run."""
import argparse
import copy
import hashlib
import importlib.util
import json
import sys
from fractions import Fraction as Q
from pathlib import Path

sys.dont_write_bytecode = True
ROOTS = {
 'current': ('HDBLAST_COUPLED_MEAN_PILOT_20261025', 'b932157f308c9912fa83ee06990f1132a4e75a498aafd66eb99dabcc30b256de'),
 'I': ('HDBLAST_INPUT_PRECISION_20261023', 'c2e7bd40e94422c5fd98820941231c0da63a8392b14a5a14c4340eabdde37f61'),
 'P': ('HDBLAST_COMMON_GAUGE_FORCING_20261016', '8bf173293b5cecf7f78e138ac5004abd9c7b1821ef29addb7ecfbc2a6409b4f6'),
 'N': ('HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016', '928f3e4e58d222a3c7e1274babf9e94930ec501b74cf6ba04209f97c6bbded44'),
 'G': ('HDBLAST_GEOMETRIC_BULK_SOURCES_20261016', '3ff1f95c151d2da6a1ab5e9b54718c2bb325c10d9c62a7ed9f97ef2bb36f0bb2'),
 'C': ('HDBLAST_FORCED_RESPONSE_BOUNDS_20261023', '72924328a4e79cc77f0da90999b8dd56956bfdebd5e614e15c3c0498c3a8ae76'),
 'T': ('HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023', 'cffb213a69f7f9a27932e990ffdecb7569d1007e5e5baeda708d8dc39d0cce72'),
}

def need(ok, marker):
    if not ok:
        raise ValueError(marker)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def authenticate(repo):
    report = {}
    for key, (name, pin) in ROOTS.items():
        root = repo / 'research' / name
        need(root.is_dir() and not root.is_symlink(), 'ROOT_TYPE_' + key)
        manifest = root / 'MANIFEST.json'
        need(digest(manifest) == pin, 'MANIFEST_' + key)
        seen = {'MANIFEST.json'}
        for item in json.loads(manifest.read_text())['files']:
            rel = Path(item['path'])
            need(not rel.is_absolute() and '..' not in rel.parts and str(rel) not in seen,
                 'MEMBER_PATH_' + key)
            seen.add(str(rel))
            path = root / rel
            need(path.is_file() and not path.is_symlink(), 'MEMBER_TYPE_' + key)
            need(path.stat().st_size == item['bytes'] and digest(path) == item['sha256'],
                 'MEMBER_BYTES_' + key)
        for path in root.rglob('*'):
            need(not path.is_symlink(), 'TREE_SYMLINK_' + key)
            if path.is_file():
                need(str(path.relative_to(root)) in seen, 'UNLISTED_FILE_' + key)
        report[key] = {'manifest_sha256': pin, 'files_including_manifest': len(seen)}
    return report

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

class L:
    """Minimal independent exact Laurent algebra, with monomial division only."""
    def __init__(self, value=0):
        self.terms = ({k: Q(v) for k, v in value.items() if v} if isinstance(value, dict)
                      else ({(): Q(value)} if value else {}))
    @staticmethod
    def coerce(value): return value if isinstance(value, L) else L(value)
    def __add__(self, other):
        result = self.terms.copy()
        for key, value in L.coerce(other).terms.items(): result[key] = result.get(key, 0) + value
        return L(result)
    __radd__ = __add__
    def __neg__(self): return L({key: -value for key, value in self.terms.items()})
    def __sub__(self, other): return self + -L.coerce(other)
    def __rsub__(self, other): return L.coerce(other) + -self
    def __mul__(self, other):
        result = {}
        for a, v in self.terms.items():
            for b, w in L.coerce(other).terms.items():
                powers = dict(a)
                for name, power in b: powers[name] = powers.get(name, 0) + power
                key = tuple(sorted((name, power) for name, power in powers.items() if power))
                result[key] = result.get(key, 0) + v*w
        return L(result)
    __rmul__ = __mul__
    def __truediv__(self, other):
        terms = L.coerce(other).terms
        need(len(terms) == 1, 'LAURENT_MONOMIAL_DIVISOR')
        (key, value), = terms.items()
        return self * L({tuple((name, -power) for name, power in key): 1/value})
    def __rtruediv__(self, other): return L.coerce(other) / self
    def __eq__(self, other): return self.terms == L.coerce(other).terms
    def __bool__(self): return bool(self.terms)
    def json(self):
        return [{'coefficient': str(value), 'powers': dict(key)}
                for key, value in sorted(self.terms.items())]

def symbols(names): return {name: L({((name, 1),): Q(1)}) for name in names.split()}

def reflection(tensor):
    # Pullback under diag(1,1,1,-1): coordinate monomial AND covector indices.
    return {indices: {powers: value * (-1)**(powers[3] + indices.count(3))
                      for powers, value in poly.items()} for indices, poly in tensor.items()}

def connected(indices, matrices):
    # Strong connectivity is only coordinate indecomposability, not a claim
    # that no other similarity transform can further split a block.
    for start in indices:
        reached = {start}
        while True:
            more = {i for matrix in matrices for i in indices
                    if any(matrix[i][j] for j in reached)}
            if more <= reached: break
            reached |= more
        if reached != set(indices): return False
    return True

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--repo-root', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--mutation', choices=['cross_block', 'tt_boundary', 'tt_bulk_parity', 'reflection'])
    args = p.parse_args()
    repo = args.repo_root.resolve()
    need(not args.output.exists(), 'FRESH_OUTPUT_REQUIRED')
    need(repo not in args.output.resolve().parents, 'OUTPUT_INSIDE_FROZEN_REPO')
    before = authenticate(repo)
    roots = {key: repo / 'research' / name for key, (name, pin) in ROOTS.items()}
    op = load('structure_frozen_operator', roots['G'] / 'projection/inputs/linear_operator.py')
    angular = load('structure_frozen_angular', roots['N'] / 'angular/angular_exact.py')
    adapter = load('structure_frozen_adapters', roots['N'] / 'angular/source-audit/source_adapters.py')
    matrices = op.matrices()
    source_path = roots['G'] / 'projection/GEOMETRIC_COEFFICIENT_MAPS.json'
    source = json.loads(source_path.read_text())
    groups = [[op.IDX[f'{name}:0'] for name in ('K','L','N','J','w')]]
    groups += [[op.IDX[f'{name}:{ell}'] for name in ('K','L','N','J','w','v','u','H')]
               for ell in (2,4)]
    groups += [[op.IDX['T:+']], [op.IDX['T:-']]]
    if args.mutation == 'cross_block': matrices['O'][0][1] += 1
    if args.mutation == 'tt_boundary':
        matrices['By_a'][21][21] += 1
        matrices['By_a'][22][22] += 1
    if args.mutation == 'tt_bulk_parity':
        source['maps']['metric_forcing']['ab'][6].append({'coefficient':'1','powers':{}})
    membership = {index: group for group, indices in enumerate(groups) for index in indices}
    need(sorted(membership) == list(range(23)), 'BLOCK_PARTITION')
    need(all(not value or membership[i] == membership[j]
             for matrix in matrices.values() for i, row in enumerate(matrix)
             for j, value in enumerate(row)), 'PRIMITIVE_BLOCK_INVARIANCE')
    need(all(connected(indices, matrices.values()) for indices in groups), 'COORDINATE_BLOCK_CONNECTIVITY')
    blocks = []
    for name, indices in zip(('l0','l2','l4','TTplus','TTminus'), groups):
        blocks.append({'name':name,'indices':indices,'labels':[op.LABELS[i] for i in indices],
                       'dimension':len(indices),
                       'Dirichlet':sum(int(matrices['PD'][i][i]) for i in indices),
                       'Robin':sum(int(matrices['PN'][i][i]) for i in indices)})
    need([b['dimension'] for b in blocks] == [5,8,8,1,1], 'BLOCK_DIMENSIONS')
    need([b['Dirichlet'] for b in blocks] == [1,2,2,0,0], 'BOUNDARY_PARITY_COUNTS')
    need(all(matrix[21][21] == matrix[22][22] for matrix in matrices.values()), 'TT_OPERATOR_EQUALITY')
    need(all(not matrices[key][i][j] for key in ('By_a','By_p','By_sigma1','By_sigma2')
             for i in (21,22) for j in range(23)), 'TT_NEUMANN_LINEAR_PART')

    z = symbols('ell rho a h c Ht p py b d V V1 V2 sigma1 sigma2')
    assembled = op.assemble(**z)
    tt = {key: assembled[key][21][21] for key in ('A','B','D','M')}
    expected_tt = {
      'A': -3*z['ell']*z['Ht'], 'B': -3*z['ell']*z['rho']*z['a'], 'M': L(0),
      'D': z['ell']*z['ell']*(-6/(z['c']*z['c']) + z['rho']*z['rho']*
            (-2*z['a']*z['a']+2*z['h']*z['h']-2*z['b']-8*z['d']+Q(4,3)*z['V']))}
    need(tt == expected_tt, 'INDEPENDENT_TT_SCALAR_EQUATION')
    tt_maps = {}
    for name, fields in source['maps'].items():
        if isinstance(fields, dict) and 'ab' in fields:
            need(len(fields['ab']) == 7 and fields['ab'][5] == fields['ab'][6], 'TT_BULK_SOURCE_PARITY_' + name)
            tt_maps[name] = {'equal':True,'terms_per_TT_slot':len(fields['ab'][5])}
    need('metric_forcing' in tt_maps and 'QEred' in tt_maps, 'COMPLETE_TT_MAPS_PRESENT')
    for name, fields in source['maps'].items():
        if isinstance(fields, dict) and set(fields) == {'yy','y0','00','ya','0a','ab'}:
            need({key:len(value) for key,value in fields.items()} ==
                 {'yy':3,'y0':3,'00':3,'ya':2,'0a':2,'ab':7}, 'SOURCE_SPACE_SHAPE_' + name)
    shell_z = symbols('a p sigmap sigmapp sigmappp h X V L d v')
    shell = adapter.shell_sources(**shell_z)
    ell, rho = z['ell'], z['rho']
    expected_shell = -Q(8,3)*ell*rho*(shell_z['a']*shell_z['h']+shell_z['sigmap']*shell_z['v']/3)*shell_z['d']
    packed_shell = op.boundary_forcing(ell,rho,shell)
    need(packed_shell[21] == packed_shell[22] == expected_shell, 'TT_SHELL_SOURCE_PARITY_AND_NORMALIZATION')
    initial_z = symbols('rho a p u uy ut v vt Fmom Gmom Fmom_y f0 f2 f4 fy0 fy2 fy4')
    initial = adapter.initial_frame_data(**{k:v for k,v in initial_z.items() if not k.startswith('f')},
             f=tuple(initial_z[k] for k in ('f0','f2','f4')),fy=tuple(initial_z[k] for k in ('fy0','fy2','fy4')))
    initial_u, initial_v = op.pack_initial(initial,ell)
    need(initial_u[21] == initial_u[22] == initial_v[21] == initial_v[22] == 0, 'TT_ZERO_INITIAL_DATA')

    basis = angular.TT
    reflected = reflection(basis[0])
    if args.mutation == 'reflection': reflected = copy.deepcopy(basis[0])
    need(angular.equal(reflected,basis[1]) and angular.equal(reflection(basis[1]),basis[0]), 'ANGULAR_REFLECTION_EXCHANGES_TT')
    need(all(angular.equal(reflection(t),t) for t in angular.SCALARS+angular.GRADS+angular.HESS), 'NON_TT_REFLECTION_EVEN')
    for i, tensor in enumerate(basis):
        need(all(not p for p in angular.trace(tensor).values()), 'TT_TRACE')
        need(all(not p for p in angular.div(tensor).values()), 'TT_DIVERGENCE')
        need(angular.equal(angular.lap(tensor),angular.tscale(-6,tensor)), 'TT_ROUGH_LAPLACIAN')
    need(angular.inner(basis[0],basis[0]) == angular.inner(basis[1],basis[1]) == Q(1,6)
         and angular.inner(basis[0],basis[1]) == 0, 'TT_GRAM')
    # E maps 22 active coordinates to all 23 slots. R averages the TT pair.
    embedding = [[Q(i == j or i == 22 and j == 21) for j in range(22)] for i in range(23)]
    restriction = [[Q(i == j) for j in range(23)] for i in range(22)]
    restriction[21][21] = restriction[21][22] = Q(1,2)
    need(op.mm(restriction,embedding) == op.ident(22), 'ACTIVE_LEFT_INVERSE')
    for key, matrix in matrices.items():
        action = op.mm(matrix,embedding)
        need(action == op.mm(embedding,op.mm(restriction,action)), 'ACTIVE_INTERTWINING_' + key)
    gram = op.gram()
    active_gram = [sum(gram[i]*embedding[i][j]*embedding[i][j] for i in range(23)) for j in range(22)]
    need(active_gram[-1] == Q(1,3), 'ACTIVE_TT_GRAM')
    # Verify the separate reflection in the NORMAL coordinate of the shell.
    signs = [1-2*int(matrices['PD'][i][i]) for i in range(23)]
    need([signs[i] for i in (21,22)] == [1,1], 'TT_NORMAL_SHELL_PARITY')
    after = authenticate(repo)
    need(after == before, 'FROZEN_CLOSURES_CHANGED')
    result = {
      'schema':'HDBLAST_ANGULAR_STRUCTURE_EXACT_V1', 'status':'PASS_EXACT_BLOCKS_AND_SELECTED_PARITY',
      'arithmetic':'Python standard-library Fraction and independent Laurent algebra; no angular sampling',
      'frozen_closures':before, 'blocks':blocks, 'primitive_matrices_checked':len(matrices),
      'ambient_dimension':23, 'selected_reflection_even_dimension':22,
      'selected_active_blocks':[5,8,8,1], 'selected_Dirichlet':5,'selected_Robin':17,
      'active_labels':op.LABELS[:21]+['T:sum'],
      'active_embedding_sparse':[[i,j,str(v)] for i,row in enumerate(embedding) for j,v in enumerate(row) if v],
      'active_restriction_sparse':[[i,j,str(v)] for i,row in enumerate(restriction) for j,v in enumerate(row) if v],
      'active_gram_diagonal':list(map(str,active_gram)),
      'active_intertwining_all_primitive_matrices':True,
      'normal_shell_parities':dict(zip(op.LABELS,signs)),
      'TT_primitive_nonzero_diagonals':{key:str(matrix[21][21]) for key,matrix in matrices.items() if matrix[21][21]},
      'TT_equation_coefficients':{key:value.json() for key,value in tt.items()},
      'TT_shell_forcing':expected_shell.json(), 'TT_initial_values_and_velocities':'identically zero',
      'TT_angular_reflection':'diag(1,1,1,-1) pullback exchanges Tplus and Tminus',
      'TT_geometry':{'trace':0,'divergence':0,'rough_laplacian_eigenvalue':-6,'each_norm_squared':'1/6','cross_inner_product':'0','Tsum_norm_squared':'1/3'},
      'TT_source_map_equalities':tt_maps,'source_maps_sha256':digest(source_path),
      'linear_operator_sha256':digest(roots['G']/'projection/inputs/linear_operator.py'),
      'conclusion':'For the registered selected data and unique formal linear response, U_Tplus=U_Tminus; the antisymmetric TT amplitude is zero. General independent second data need not obey this relation.',
      'limits':['No old numerical producer executed','No evolved full response or numerical residual computed by this audit','No nonlinear family or physical instability established','Coordinate block connectivity does not assert absolute algebraic irreducibility']}
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':result['status'],'blocks':result['selected_active_blocks'],'output':str(args.output)}))

if __name__ == '__main__': main()
