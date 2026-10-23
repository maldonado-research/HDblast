"""Exact all-point projection of complete geometric sources into 23 slots.

Candidate coefficients come from a rank-complete rational coefficient witness.
Acceptance always additionally checks every reduced polynomial coefficient in
every ambient tensor component. A finite set of angular values is never used.
"""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
import angular_generic as a

LABELS = ([f'{field}:{degree}' for field in ('K', 'L', 'N', 'J', 'w')
           for degree in (0, 2, 4)]
          + [f'{field}:{degree}' for field in ('v', 'u') for degree in (2, 4)]
          + ['H:2', 'H:4', 'T:+', 'T:-'])
FIELD_LABELS = {
    'yy': ('K:0', 'K:2', 'K:4'), 'y0': ('L:0', 'L:2', 'L:4'),
    '00': ('N:0', 'N:2', 'N:4'), 'ya': ('v:2', 'v:4'),
    '0a': ('u:2', 'u:4'),
    'ab': ('J:0', 'J:2', 'J:4', 'H:2', 'H:4', 'T:+', 'T:-'),
    'scalar': ('w:0', 'w:2', 'w:4'),
}
BASES = {'scalar': a.SCALARS, 'vector': a.GRADS, 'tensor': a.TENSORS}
FIELD_BASES = {'yy': 'scalar', 'y0': 'scalar', '00': 'scalar',
               'ya': 'vector', '0a': 'vector', 'ab': 'tensor', 'scalar': 'scalar'}


class ProjectionError(ValueError):
    pass


def canonical(tensor):
    return {key: a.reduce(poly) for key, poly in tensor.items()}


def inverse(matrix):
    n = len(matrix)
    rows = [[Q(value) for value in row] + [Q(i == j) for j in range(n)]
            for i, row in enumerate(matrix)]
    for col in range(n):
        pivot = next((i for i in range(col, n) if rows[i][col]), None)
        if pivot is None: raise ProjectionError('ANGULAR_BASIS_RANK_DEFICIENT')
        rows[col], rows[pivot] = rows[pivot], rows[col]
        weight = rows[col][col]
        rows[col] = [value / weight for value in rows[col]]
        for i in range(n):
            if i != col and rows[i][col]:
                weight = rows[i][col]
                rows[i] = [x - weight * y for x, y in zip(rows[i], rows[col])]
    return tuple(tuple(row[n:]) for row in rows)


@lru_cache(None)
def witnesses(kind):
    """Select independent MONOMIAL coefficients; compute their exact inverse."""
    basis = BASES[kind]
    keys = sorted({(component, exponent) for tensor in basis
                   for component, poly in tensor.items() for exponent in poly})
    selected, echelon, matrix = [], [], []
    for key in keys:
        component, exponent = key
        row = [tensor[component].get(exponent, Q(0)) for tensor in basis]
        residual = row.copy()
        for pivot, old in echelon:
            weight = residual[pivot]
            if weight:
                residual = [x - weight * y for x, y in zip(residual, old)]
        pivot = next((i for i, value in enumerate(residual) if value), None)
        if pivot is not None:
            weight = residual[pivot]
            echelon.append((pivot, [value / weight for value in residual]))
            selected.append(key)
            matrix.append(row)
            if len(selected) == len(basis): break
    if len(selected) != len(basis):
        raise ProjectionError('ANGULAR_BASIS_RANK_DEFICIENT_' + kind)
    return tuple(selected), inverse(matrix)


def reconstruct(coefficients, kind):
    basis = BASES[kind]
    if len(coefficients) != len(basis):
        raise ProjectionError('ANGULAR_COEFFICIENT_SHAPE_' + kind)
    return a.tadd(*(a.tscale(value, tensor)
                    for value, tensor in zip(coefficients, basis)))


def project_tensor(tensor, kind, label='ANGULAR'):
    """Return exact coefficients ONLY after complete polynomial reconstruction."""
    basis = BASES[kind]
    if set(tensor) != set(basis[0]):
        raise ProjectionError(label + '_COMPONENT_SHAPE')
    tensor = canonical(tensor)
    rows, inv = witnesses(kind)
    values = [tensor[component].get(exponent, 0) for component, exponent in rows]
    coefficients = tuple(sum((weight * value for weight, value in zip(row, values)), Q(0))
                         for row in inv)
    recovered = reconstruct(coefficients, kind)
    for component in tensor:
        remainder = a.reduce(a.add(tensor[component], a.scale(-1, recovered[component])))
        if remainder:
            raise ProjectionError(label + '_ALL_POINT_RECONSTRUCTION_' + kind)
    return coefficients


def project_fields(fields):
    if set(fields) - set(FIELD_BASES):
        raise ProjectionError('UNKNOWN_GEOMETRIC_FIELD')
    return {name: project_tensor(tensor, FIELD_BASES[name], name)
            for name, tensor in fields.items()}


def tensor_fields(tensor):
    """Full covariant 6-slot tensor -> scalar/vector/S3 tensor fields.

    Indices are 0=y, 1=orthonormal time, 2..5=redundant ambient tangent.
    No frame rescaling or raising of a time index is performed here.
    """
    if set(tensor) != set(product(range(6), repeat=2)):
        raise ProjectionError('FULL_GEOMETRIC_TENSOR_SHAPE')
    if any(a.reduce(a.add(tensor[i, j], a.scale(-1, tensor[j, i])))
           for i, j in product(range(6), repeat=2)):
        raise ProjectionError('FULL_GEOMETRIC_TENSOR_SYMMETRY')
    return {'yy': {(): tensor[0, 0]}, 'y0': {(): tensor[0, 1]},
            '00': {(): tensor[1, 1]},
            'ya': {(i,): tensor[0, i + 2] for i in range(4)},
            '0a': {(i,): tensor[1, i + 2] for i in range(4)},
            'ab': {(i, j): tensor[i + 2, j + 2]
                   for i, j in product(range(4), repeat=2)}}


def covector_fields(covector):
    if set(covector) != {(i,) for i in range(6)}:
        raise ProjectionError('FULL_GEOMETRIC_COVECTOR_SHAPE')
    return {'y': project_tensor({(): covector[0,]}, 'scalar', 'QC_y'),
            '0': project_tensor({(): covector[1,]}, 'scalar', 'QC_0'),
            'a': project_tensor({(i,): covector[i + 2,] for i in range(4)},
                                'vector', 'QC_a')}


def project_geometric(sources):
    """Project engine outputs; every accepted geometric sector is exact."""
    tensor_names = ('QRic', 'QE', 'QE_matter', 'QEred', 'metric_forcing', 'QEred_general',
                    'gauge_gradient', 'PC1')
    out = {}
    for name in tensor_names:
        if name in sources: out[name] = project_fields(tensor_fields(sources[name]))
    for name in ('QC', 'QC_actual', 'QC_background'):
        if name in sources: out[name] = covector_fields(sources[name])
    for name in ('Qphi', 'Qphi_general', 'scalar_forcing'):
        if name in sources:
            out[name] = project_tensor({(): sources[name]}, 'scalar', name)
    return out


def pack23(fields):
    if set(fields) - set(FIELD_LABELS):
        raise ProjectionError('UNKNOWN_PACK23_FIELD')
    values = {label: Q(0) for label in LABELS}
    for name, coefficients in fields.items():
        labels = FIELD_LABELS[name]
        if len(coefficients) != len(labels):
            raise ProjectionError('PACK23_FIELD_SHAPE_' + name)
        values.update(zip(labels, coefficients))
    return [values[label] for label in LABELS]


def certificate():
    """Serializable independence and basis convention certificate."""
    return {'proof': 'canonical quotient coefficient witness plus full remainder',
            'quotient': 'n3^2=1-n0^2-n1^2-n2^2',
            'angular_coordinates': ['n0', 'n1', 'n2', 'n3'],
            'fields': {name: list(labels) for name, labels in FIELD_LABELS.items()},
            'labels23': LABELS,
            'basis': {kind: {'rank': len(BASES[kind]),
                             'witness_coefficients': [
                                 {'component': list(component), 'exponent': list(exponent)}
                                 for component, exponent in witnesses(kind)[0]],
                             'witness_inverse': [[str(value) for value in row]
                                                 for row in witnesses(kind)[1]]}
                      for kind in BASES}}
