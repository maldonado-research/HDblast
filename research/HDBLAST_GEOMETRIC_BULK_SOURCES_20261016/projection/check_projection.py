#!/usr/bin/env python3
"""Exact all-point projection replay with independent coefficient symbols."""
from pathlib import Path
from fractions import Fraction as Q
import argparse
import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
sys.path.insert(0, str(Path(__file__).resolve().parent))
import angular_generic as a
import angular_project as projection
from formal_coefficients import P, symbol, formal_inputs

HERE = Path(__file__).resolve().parent
CONTROLS = {
    'tt_deleted': 'BASIS_tensor_ALL_POINT_RECONSTRUCTION_tensor',
    'one_tt_orientation_deleted': 'BASIS_tensor_ALL_POINT_RECONSTRUCTION_tensor',
    'tensor_factor_two': 'BASIS_tensor_ALL_POINT_RECONSTRUCTION_tensor',
    'witness_only_acceptance': 'OFFSPAN_REJECTION_REQUIRED',
    'tt_packing_swapped': 'REGISTERED23_PACKING',
}


def require(condition, marker):
    if not condition: raise ValueError(marker)


def load_input(name):
    spec = importlib.util.spec_from_file_location('pinned_' + name, HERE / 'inputs' / (name + '.py'))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate(mutation=None):
    parent = load_input('angular_exact')
    operator = load_input('linear_operator')
    for name in ('Y', 'Z', 'F', 'P', 'SCALARS', 'GRADS', 'HESS', 'TT', 'TENSORS'):
        require(getattr(a, name) == getattr(parent, name), 'PARENT_BASIS_' + name)
    if mutation is None:
        require(a.primitive_matrices() == parent.primitive_matrices(), 'GENERIC_BACKEND_PARENT_OPERATIONS')
    require(projection.LABELS == operator.LABELS, 'REGISTERED23_LABELS')
    original_reconstruct = projection.reconstruct
    if mutation in ('tt_deleted', 'one_tt_orientation_deleted', 'tensor_factor_two'):
        def broken(coefficients, kind):
            coefficients = list(coefficients)
            if kind == 'tensor':
                if mutation == 'tt_deleted': coefficients[-2:] = [0, 0]
                elif mutation == 'one_tt_orientation_deleted': coefficients[-1] = 0
                else: coefficients[3] *= 2
            return original_reconstruct(coefficients, kind)
        projection.reconstruct = broken
    if mutation == 'witness_only_acceptance':
        def witness_only(tensor, kind, label='ANGULAR'):
            tensor = projection.canonical(tensor)
            rows, inv = projection.witnesses(kind)
            values = [tensor[component].get(exponent, 0) for component, exponent in rows]
            return tuple(sum((weight*value for weight, value in zip(row, values)), Q(0))
                         for row in inv)
        projection.project_tensor = witness_only
    # Every basis coefficient is an independent symbol. Equality is coefficient
    # equality in the Laurent ring, not a substitution at selected jet values.
    basis_counts = {}
    for kind, basis in projection.BASES.items():
        coefficients = tuple(symbol(kind + str(i)) for i in range(len(basis)))
        tensor = original_reconstruct(coefficients, kind)
        recovered = projection.project_tensor(tensor, kind, 'BASIS_' + kind)
        require(recovered == coefficients, 'FORMAL_BASIS_RECOVERY_' + kind)
        basis_counts[kind] = len(basis)
    # Deliberately off-span odd scalar is invisible to every scalar coefficient
    # witness, but must fail the COMPLETE polynomial remainder test.
    offspan = {(): a.X[0]}
    accepted = False
    try: projection.project_tensor(offspan, 'scalar', 'OFFSPAN')
    except projection.ProjectionError as error:
        require('ALL_POINT_RECONSTRUCTION' in str(error), 'OFFSPAN_REASON')
    else: accepted = True
    require(not accepted, 'OFFSPAN_REJECTION_REQUIRED')
    fields = {name: tuple(symbol(label) for label in labels)
              for name, labels in projection.FIELD_LABELS.items()}
    packed = projection.pack23(fields)
    if mutation == 'tt_packing_swapped': packed[-2:] = packed[-2:][::-1]
    require(packed == operator.pack_fields(fields), 'REGISTERED23_PACKING')
    # Quotient reduction vanishes exactly for all real points on S3.
    sphere_relation = a.add(*(a.mul(x, x) for x in a.X), a.scale(-1, a.ONE))
    require(not a.reduce(sphere_relation), 'SPHERE_QUOTIENT_RELATION')
    m, scalar, background = formal_inputs()
    require(len(m) == 42 and len(scalar) == 6 and len(background) == 11,
            'INDEPENDENT_FORMAL_INPUT_SCHEMA')
    return {'basis_ranks': basis_counts, 'independent_metric_jet_symbols': len(m),
            'independent_scalar_jet_symbols': len(scalar),
            'independent_background_symbols': len(background),
            'full_reconstruction': 'all components, every canonical monomial',
            'offspan_witness_counterexample_rejected': True,
            'parent_primitive_matrices_equal': True,
            'registered_23_packing_equal': True}


def pins():
    paths = sorted(HERE.glob('*.py')) + sorted((HERE / 'inputs').iterdir())
    return {str(path.relative_to(HERE)): {
        'bytes': path.stat().st_size,
        'sha256': hashlib.sha256(path.read_bytes()).hexdigest()} for path in paths}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--run-controls', action='store_true')
    parser.add_argument('--mutation', choices=tuple(CONTROLS))
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists(): raise ValueError('OUTPUT_MUST_BE_NEW')
    if HERE == output or HERE in output.parents: raise ValueError('OUTPUT_MUST_BE_EXTERNAL')
    report = validate(args.mutation)
    controls = {}
    if args.run_controls:
        with tempfile.TemporaryDirectory(prefix='hdblast-projection-controls-') as directory:
            for name, marker in CONTROLS.items():
                destination = Path(directory) / (name + '.json')
                command = [sys.executable, '-I', '-S']
                if sys.flags.optimize: command.append('-O')
                command += ['-B', str(Path(__file__).resolve()), '--output', str(destination), '--mutation', name]
                result = subprocess.run(command, capture_output=True, text=True, timeout=120)
                require(result.returncode != 0 and not destination.exists() and marker in result.stderr,
                        'CONTROL_NOT_REJECTED_' + name)
                controls[name] = {'rejected': True, 'required_stderr': marker, 'output_absent': True}
    report = {'status': 'PASS', **report, 'certificate': projection.certificate(),
              'controls': controls, 'control_count': len(controls), 'pins': pins()}
    with output.open('x', encoding='utf-8') as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print('PASS: exact all-point projection and requested controls.')


if __name__ == '__main__': main()
