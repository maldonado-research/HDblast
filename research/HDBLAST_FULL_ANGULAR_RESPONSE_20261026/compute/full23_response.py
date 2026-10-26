#!/usr/bin/env python3
"""Bounded native full23 first-feedback response and its conditional X error.

The registered rational displacement remains exact in B0. The full23 P1
velocity defines B0, and P1 u,du appear only in the exact ridge feedback.
Every scientific comparison uses Fractions and explicit checks.
"""
import argparse
import hashlib
import importlib.util
import json
import math
import sys
from fractions import Fraction as F
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEPENDENCY_PINS = {
    'initial_producer': ('initial/full23_initial.py', '2e0b3fe803df47060550a4be716095450e710627fd7c5258d2b4e2868f1beda1'),
    'initial_report': ('initial/results/INITIAL_REPORT.json', 'bbde14bab0f52dc4371ee016aba48dd6e91874f610d09c197a942f8c8f2e1fad'),
    'initial_profiles': ('initial/results/INITIAL_PROFILES.json', '555f911f6c336862dfacb6aa56ef5179ac06701812d593b5bfbbf6ea9bd4d5d1'),
    'kernel': ('compute/angular_ridge.py', 'd6e935bd958512cb032e1ef3d9503d6fa159418bcf9b28051de6d2f43a43c02f'),
}
Q = F(3409, 8000000)
SOURCE = F(52162547172519, 2)
TARGET = F(10**29)
BLOCK_NAMES = ('degree0', 'degree2', 'degree4', 'TTplus', 'TTminus')


def need(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def authenticate_local():
    result = {}
    for key, (relative, expected) in DEPENDENCY_PINS.items():
        path = ROOT/relative
        need(path.is_file() and not path.is_symlink(), 'LOCAL_DEPENDENCY_TYPE:'+key)
        need(sha(path) == expected, 'LOCAL_DEPENDENCY_HASH:'+key)
        result[key] = {'path': relative, 'sha256': expected, 'bytes': path.stat().st_size}
    return result


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {key: encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    return value


def sqrt_upper(value):
    need(value >= 0, 'NEGATIVE_NORM_SQUARE')
    if not value:
        return F(0)
    grid = 1 << 112
    return F(math.isqrt((value*grid*grid).__ceil__())+1, grid)


class Response:
    def __init__(self, repo_root):
        self.local_before = authenticate_local()
        initial = load_module('full23_initial_dependency', ROOT/'initial/full23_initial.py')
        # Context verifies seven complete upstream roots before project imports.
        self.context = initial.Context(repo_root)
        self.initial = initial
        self.kernel = load_module('full23_ridge_dependency', HERE/'angular_ridge.py')
        self.profile_data = json.loads((ROOT/'initial/results/INITIAL_PROFILES.json').read_text())
        self.initial_report = json.loads((ROOT/'initial/results/INITIAL_REPORT.json').read_text())
        self.I, self.ell = self.context.I, self.context.ell
        self.labels = self.context.op.LABELS
        self.weights = self.context.weights
        need(tuple(self.labels) == self.kernel.LABELS and tuple(self.weights) == self.kernel.GRAM,
             'KERNEL_OPERATOR_LABEL_GRAM')
        self.primitive = self.context.op.matrices()  # Build exact primitives once.
        self.blocks = [self.context.mean_ids,
                       [self.context.op.IDX[name+':2'] for name in ('K', 'L', 'N', 'J', 'w', 'v', 'u', 'H')],
                       [self.context.op.IDX[name+':4'] for name in ('K', 'L', 'N', 'J', 'w', 'v', 'u', 'H')],
                       [self.context.op.IDX['T:+']], [self.context.op.IDX['T:-']]]
        need(sorted(i for block in self.blocks for i in block) == list(range(23)), 'BLOCK_PARTITION')
        self.block_of = {i: b for b, block in enumerate(self.blocks) for i in block}
        for name, matrix in self.primitive.items():
            need(all(matrix[i][j] == 0 for i in range(23) for j in range(23)
                     if self.block_of[i] != self.block_of[j]), 'PRIMITIVE_BLOCK_CLOSURE:'+name)
        self._assemble_and_freeze()

    def matrix_sum(self, *terms):
        # Skip exact primitive zeros, preserving the exact orthogonal sparsity.
        out = [[self.I(0) for _ in range(23)] for _ in range(23)]
        for coefficient, matrix in terms:
            for i, row in enumerate(matrix):
                for j, value in enumerate(row):
                    if isinstance(value, self.I):
                        if value.lo == value.hi == 0:
                            continue
                    elif value == 0:
                        continue
                    out[i][j] = out[i][j]+coefficient*value
        return out

    def assemble(self, x, r):
        x, r = self.I.cv(x), self.I.cv(r)
        need(0 <= x.lo <= x.hi <= F(3, 2) and 0 <= r.lo <= r.hi <= F(7, 4), 'COEFFICIENT_DOMAIN')
        state = self.context.actual_background(r)
        rho, eta, p, a = [state[key] for key in ('rho', 'eta', 'p', 'a')]
        need(rho.lo > 0 and p.hi < 0, 'ACTUAL_BACKGROUND_SIGNS')
        t = self.ell*x
        ht = self.I(t.lo-t.lo**3/3, t.hi)
        invc = self.I(1-t.hi**2/2, 1)
        need(invc.lo > 0, 'POSITIVE_SECH_ENCLOSURE')
        h = ht/rho
        V, V1, V2 = [self.context.init.poly(self.context.init.UP[j], eta) for j in range(3)]
        d = 1/rho**2-a*a
        b = d+p*p/3
        m = self.primitive
        A = self.matrix_sum((-3*self.ell*ht, m['I']), (-2*self.ell*rho*a, m['O']),
                            (self.ell*rho*p, m['Aphi']))
        B = self.matrix_sum((-3*self.ell*rho*a, m['I']), (self.ell*rho*p, m['Bphi']))
        C = self.matrix_sum((a*a, m['Sa']), (-a*a, m['O2']),
                            (a*h, m['Sah']), (-3*a*h, m['O']), (h*h, m['Sh']),
                            (-2*b, m['Rb']), (-2*d, m['Rd']),
                            (F(4, 3)*V, m['Pmet']), (F(4, 3)*V1, m['Pg']),
                            (-V2, m['Pww']), (-V1, m['W_K']), (3*p*h, m['W_L']))
        D = self.matrix_sum((self.ell**2*invc**2, m['Delta']),
                            (2*self.ell**2*rho*invc*a, m['Ca']),
                            (2*self.ell**2*rho*invc*h, m['Ch']),
                            (self.ell**2*rho**2, C),
                            (self.ell**2*rho*p*invc, m['Ephi']))
        M = self.matrix_sum((-self.ell*rho*a, m['By_a']),
                            (-self.ell*rho*p, m['By_p']),
                            (2*self.ell*rho*p, m['By_sigma1']),
                            (-4*self.ell*rho*(1+eta), m['By_sigma2']))
        return {'A': A, 'B': B, 'D': D, 'M': M}

    def block_norms(self, matrix):
        need(all(matrix[i][j] == 0 for i in range(23) for j in range(23)
                 if self.block_of[i] != self.block_of[j]), 'NORM_MATRIX_BLOCK_CLOSURE')
        return [sqrt_upper(sum((self.weights[i]/self.weights[j]*matrix[i][j]**2
                                for i in block for j in block), F(0))) for block in self.blocks]

    def matrix_norm(self, matrix):
        return max(self.block_norms(matrix))

    def _assemble_and_freeze(self):
        local = self.assemble(self.I(0), self.I(0))
        self.actual_boxes = self.assemble(self.I(0, F(3, 2)), self.I(0, F(7, 4)))
        self.actual_boxes['M'] = local['M']  # Robin coefficients are shell constants.
        self.matrices = {name: [[(box.lo+box.hi)/2 for box in row] for row in matrix]
                         for name, matrix in local.items()}
        old = self.context.frozen_profiles['matrices']
        self.mean_entries_reused = 0
        for name in ('A', 'B', 'D', 'M'):
            for i0, i in enumerate(self.context.mean_ids):
                for j0, j in enumerate(self.context.mean_ids):
                    self.matrices[name][i][j] = F(old[name][i0][j0])
                    self.mean_entries_reused += 1
        self.coefficient_entry_errors = {
            name: [[(box-self.matrices[name][i][j]).mag() for j, box in enumerate(row)]
                   for i, row in enumerate(matrix)]
            for name, matrix in self.actual_boxes.items()}
        self.matrix_caps = {name: self.matrix_norm(matrix) for name, matrix in self.matrices.items()}
        self.coefficient_errors = {name: self.matrix_norm(matrix)
                                   for name, matrix in self.coefficient_entry_errors.items()}
        self.matrix_block_caps = {name: self.block_norms(matrix) for name, matrix in self.matrices.items()}
        self.error_block_caps = {name: self.block_norms(matrix)
                                for name, matrix in self.coefficient_entry_errors.items()}
        parity = self.kernel.PARITY
        need(all(self.matrices['M'][i][j] == 0 for i in range(23) if parity[i] == -1 for j in range(23)),
             'ROBIN_DIRICHLET_ROW_MIXING')
        need(all(self.matrices['M'][i][j] == 0 for j in range(23) if parity[j] == -1 for i in range(23)),
             'ROBIN_DIRICHLET_COLUMN_MIXING')

    def ledger(self):
        nominal_u, nominal_du, nominal_v = map(F, self.initial_report['nominal_H_caps'])
        delta_u, delta_du, delta_v = map(F, self.initial_report['total_H_caps'])
        eps_u, eps_du, eps_v = map(F, self.initial_report['P1_H_error_caps'])
        beta_nominal = max((nominal_u+F(3, 2)*nominal_v)/3, nominal_du+nominal_v)
        X0 = max((nominal_u+F(3, 2)*(nominal_v+eps_v))/3, nominal_du+nominal_v+eps_v)
        efree = max((delta_u+F(3, 2)*delta_v)/3, delta_du+delta_v)
        e0 = efree+eps_v+SOURCE
        need(0 <= Q < 1, 'FULL_H_CONTRACTION')
        Ucap = (X0+e0)/(1-Q)
        ea, eb, ed, em = [self.coefficient_errors[name] for name in ('A', 'B', 'D', 'M')]
        na, nb, nd, nm = [self.matrix_caps[name] for name in ('A', 'B', 'D', 'M')]
        df_coeff = (ea+eb+F(3, 4)*ed)*X0
        df_profile = (na+nb)*eps_du+nd*eps_u
        dg_coeff = F(1, 4)*em*X0
        dg_profile = nm*eps_u
        df, dg = df_coeff+df_profile, dg_coeff+dg_profile
        e1 = F(3, 2)*df+dg
        tail = Q**2*Ucap
        total = tail+(1+Q)*e0+e1
        sharp = sqrt_upper(F(321, 19))
        return {'nominal_u_H': nominal_u, 'nominal_du_H': nominal_du, 'nominal_v_H': nominal_v,
                'delta_u_H': delta_u, 'delta_du_H': delta_du, 'delta_v_H': delta_v,
                'epsilon_u_H': eps_u, 'epsilon_du_H': eps_du, 'epsilon_v_H': eps_v,
                'beta_nominal_free_X_upper': beta_nominal, 'B0_X_upper': X0,
                'retained_full_H_contraction': Q, 'free_input_error_X_upper': efree,
                'omitted_full_H_forcing_X_upper': SOURCE, 'e0': e0,
                'actual_response_X_upper': Ucap,
                'delta_f_coefficient': df_coeff, 'delta_f_profile': df_profile,
                'delta_g_coefficient': dg_coeff, 'delta_g_profile': dg_profile,
                'delta_f': df, 'delta_g': dg, 'e1': e1,
                'two_term_tail': tail, 'total_X_error_upper': total,
                'selected_evaluation_factor_upper': sharp,
                'zero_wedge_pointwise_value_error_upper': F(3, 4)*sharp*total,
                'pointwise_scaled_derivative_error_upper': sharp*total}

    def reflected_nominal(self, r):
        r = F(r)
        values, _ = self.context.profiles(self.I(abs(r)))
        negative = r < 0
        u = [box*(sign if negative else 1) for box, sign in zip(values[0], self.kernel.PARITY)]
        du = [box*(-sign if negative else 1) for box, sign in zip(values[1], self.kernel.PARITY)]
        return u, du

    def free_wave(self, candidate, x, r):
        a, b = r+x, r-x
        ua, da = self.reflected_nominal(a)
        ub, db = self.reflected_nominal(b)
        va, vb = candidate.v.primitive(a), candidate.v.primitive(b)
        pa, pb = candidate.v.primitive(a, 1), candidate.v.primitive(b, 1)
        value = [(ua[i]+ub[i]+pa[i]-pb[i])/2 for i in range(23)]
        dx = [(da[i]-db[i]+va[i]+vb[i])/2 for i in range(23)]
        dr = [(da[i]+db[i]+va[i]-vb[i])/2 for i in range(23)]
        return {'value': value, 'dx': dx, 'dr': dr}

    def run(self, output_dir):
        output_dir = Path(output_dir).resolve()
        need(not output_dir.exists(), 'OUTPUT_DIRECTORY_MUST_BE_NEW')
        need(all(root not in output_dir.parents for root in self.context.roots.values()), 'OUTPUT_INSIDE_FROZEN_ROOT')
        ledger = self.ledger()
        passed = ledger['total_X_error_upper'] < TARGET
        profiles = dict(self.profile_data)
        profiles['schema'] = 'HDBLAST_FULL23_NATIVE_RESPONSE_PROFILES_V1'
        profiles['matrices'] = self.matrices
        profiles['candidate_definition'] = 'B0=free(exact_registered_full23_u,vP1); B1=K(ridge(uP1,duP1,vP1))+L(M PN(uP1+primitive(vP1))).'
        candidate = self.kernel.from_data(profiles)
        TT = [self.context.op.IDX['T:+'], self.context.op.IDX['T:-']]
        need(all(F(row[i]) == 0 for name in ('u', 'du', 'v') for row in profiles['profiles'][name] for i in TT),
             'TT_INITIAL_PROFILE_ZERO')
        need(all(self.matrices[name][i][j] == 0 for name in ('A', 'B', 'D', 'M') for i in TT for j in range(23) if j not in TT),
             'TT_CANDIDATE_BLOCK_INVARIANCE')
        old_evaluated = json.loads((self.context.roots['current']/'results/EVALUATED_CANDIDATE.json').read_text())
        locations = [(F(row['x']), F(row['r'])) for row in old_evaluated['samples']]
        need(len(locations) == 6, 'FROZEN_SIX_LOCATIONS')
        samples, mean_feedback_checks, mean_B0_overlaps = [], 0, 0
        for k, (x, r) in enumerate(locations):
            feedback = candidate.feedback(x, r, self.matrices['M'])
            old_feedback = old_evaluated['samples'][k]['feedback_exact']
            for name in ('value', 'dx', 'dr'):
                for j, i in enumerate(self.context.mean_ids):
                    need(feedback[name][i] == F(old_feedback[name][j]), 'FROZEN_EXACT_MEAN_FEEDBACK_REGRESSION')
                    mean_feedback_checks += 1
            free = self.free_wave(candidate, x, r)
            complete = {name: [free[name][i]+feedback[name][i] for i in range(23)] for name in ('value', 'dx', 'dr')}
            complete['ordinary_t_velocity'] = [box/self.ell for box in complete['dx']]
            enclosure = {}
            for name in ('value', 'dx', 'dr'):
                factor = 3 if name == 'value' else 1
                caps = [factor*ledger['total_X_error_upper']*sqrt_upper(1/weight) for weight in self.weights]
                enclosure[name] = [(box+self.I(-cap, cap)).row() for box, cap in zip(complete[name], caps)]
            need(all(feedback[name][i] == 0 and complete[name][i].lo == complete[name][i].hi == 0
                     for name in ('value', 'dx', 'dr') for i in TT), 'TT_SAMPLE_CANDIDATE_ZERO')
            old_free = self.context.frozen_profiles['free_wave_samples'][k]['B0']
            for name in ('value', 'dx', 'dr'):
                for j, i in enumerate(self.context.mean_ids):
                    lo, hi = map(F, old_free[name][j])
                    need(max(lo, free[name][i].lo) <= min(hi, free[name][i].hi), 'FROZEN_MEAN_FREE_INTERVAL_OVERLAP')
                    mean_B0_overlaps += 1
            samples.append({'x': x, 'r': r,
                            'free_wave_interval': {name: [box.row() for box in row] for name, row in free.items()},
                            'feedback_exact': self.kernel.encoded(feedback),
                            'complete_candidate_interval_evaluation': {name: [box.row() for box in row] for name, row in complete.items()},
                            'conditional_exact_response_component_enclosure': enclosure})
            print(json.dumps({'full23_response_sample': k+1, 'x': str(x), 'r': str(r)}), flush=True)
        block_labels = {name: [self.labels[i] for i in block] for name, block in zip(BLOCK_NAMES, self.blocks)}
        report = {
            'status': 'PASS_BOUNDED_FULL23_RESPONSE_X_ERROR' if passed else 'BOUNDED_FULL23_UPPER_BOUND_ABOVE_TARGET',
            'passes_target': passed, 'target_X_error': TARGET,
            'labels': self.labels, 'gram': self.weights,
            'ell': self.ell, 'fixed_mesh_denominator': 256,
            'candidate_defined_on_whole_cone': True,
            'domain': {'x': ['0', '3/2'], 'r_min': 'max(0,x-5/4)', 'x_plus_r_max': '7/4'},
            'orthogonal_blocks': block_labels,
            'matrix_norm_method': 'Maximum of Gram-weighted Frobenius upper bounds over orthogonal5+8+8+1+1 blocks; exact block closure checked.',
            'constant_matrix_H_operator_upper': self.matrix_caps,
            'coefficient_error_H_operator_upper': self.coefficient_errors,
            'constant_matrix_block_H_operator_upper': {name: dict(zip(BLOCK_NAMES, values)) for name, values in self.matrix_block_caps.items()},
            'coefficient_error_block_H_operator_upper': {name: dict(zip(BLOCK_NAMES, values)) for name, values in self.error_block_caps.items()},
            'constant_matrices': self.matrices,
            'coefficient_entry_error_caps': self.coefficient_entry_errors,
            'actual_coefficient_boxes': {name: [[box.row() for box in row] for row in matrix] for name, matrix in self.actual_boxes.items()},
            'ledger': ledger, 'display_only': {key: float(value) for key, value in ledger.items()},
            'ledger_identity': 'E=q^2*Ucap+(1+q)*e0+e1; Ucap=(X0+e0)/(1-q); e1=3df/2+dg.',
            'candidate_norm': 'X0=max((nomU+3(nomV+epsV)/2)/3,nomUr+nomV+epsV).',
            'full_H_source_and_contraction_scope': 'Authenticated T/analytic/LOCAL_TOTAL_REVIEW.md: q=3409/8000000 and source52162547172519/2 apply in full H; zero-wedge factors3/4 bulk and1/4 Robin.',
            'coefficient_assembly': 'Actual retained Taylor background; tanh(t) in[tlo-tlo^3/3,thi], sech(t) in[1-thi^2/2,1]. Includes Delta,Ca,Ch,Ephi and all curvature/potential/metric-scalar terms. py+4ap=V1; b=d+p^2/3. Robin M evaluated only at actual shell r=0.',
            'regression': {'frozen_mean_constant_entries_reused': self.mean_entries_reused,
                           'exact_mean_feedback_coefficients_equal': mean_feedback_checks,
                           'frozen_mean_B0_interval_overlaps': mean_B0_overlaps},
            'tensor_scope': 'Both candidate TT slots are exactly zero on the whole cone by zero initial tables and autonomous TT blocks. The actual omitted forcing is covered by the full-H source bound; actual TT response is not asserted zero.',
            'pointwise_scope': 'Selected angular evaluation C*=sqrt(321/19), bounded rationally above. C1 X error and causal zero wedge give value<=3C*E/4 and scaled first derivatives<=C*E. No C2 Einstein residual claim.',
            'local_dependency_hashes': self.local_before,
            'upstream_complete_root_authentication': self.context.before,
            'complete_closures_authenticated_before_and_after': True,
            'limits': ['Conditional on retained exact branch, operator, source and formal physical completion premises.',
                       'A full23 C1 response candidate and error bound; no pointwise second-derivative residual.',
                       'Finite rational feedback approximation; remaining Robin defect is covered, not set to zero.',
                       'No nonlinear family, finite-amplitude remainder or physical-origin inference.',
                       'A failed upper target gate is not an actual error lower bound.'],
        }
        evaluated = {
            'status': 'EVALUATED_FULL23_CONTINUOUS_RESPONSE_CANDIDATE',
            'accuracy_target_passed': passed, 'total_X_error_upper': ledger['total_X_error_upper'],
            'sample_count': len(samples), 'labels': self.labels, 'gram': self.weights,
            'candidate_defined_on_whole_cone': True, 'samples': samples,
            'sample_intervals': 'Arithmetic enclosures of the candidate; separate whole-cone conditional X error is the accuracy certificate.',
        }
        self.context.reauthenticate()
        need(authenticate_local() == self.local_before, 'LOCAL_DEPENDENCIES_CHANGED')
        output_dir.mkdir(parents=True)
        for filename, value in (('RESPONSE_REPORT.json', report), ('RESPONSE_PROFILES.json', profiles),
                                ('EVALUATED_CANDIDATE.json', evaluated)):
            with (output_dir/filename).open('x') as handle:
                json.dump(encode(value), handle, indent=2, sort_keys=True)
                handle.write('\n')
        print(json.dumps({'status': report['status'], 'passes_target': passed,
                          'total_X_error_upper_display': float(ledger['total_X_error_upper']),
                          'actual_response_X_upper_display': float(ledger['actual_response_X_upper']),
                          'report_sha256': sha(output_dir/'RESPONSE_REPORT.json'),
                          'profiles_sha256': sha(output_dir/'RESPONSE_PROFILES.json'),
                          'evaluated_sha256': sha(output_dir/'EVALUATED_CANDIDATE.json')}), flush=True)
        return report, profiles, evaluated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo-root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    Response(args.repo_root).run(args.output_dir)


if __name__ == '__main__':
    main()
