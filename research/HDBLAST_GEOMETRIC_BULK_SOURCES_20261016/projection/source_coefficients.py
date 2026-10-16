"""Evaluate exported exact geometric source coefficient maps on any scalar ring.

In particular interval users evaluate these algebraically collected maps, rather
than seeking cancellation among dependent interval angular tensor expressions.
"""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent


@lru_cache(None)
def load_maps(path=None):
    path = Path(path) if path is not None else HERE / 'GEOMETRIC_COEFFICIENT_MAPS.json'
    report = json.loads(path.read_text())
    if report['status'] != 'PASS_FORMAL_ALL_POINT_GEOMETRIC_SOURCE_PROJECTION':
        raise ValueError('COEFFICIENT_MAP_STATUS')
    return report


def evaluate_encoded(terms, values, monomials=None):
    monomials = {} if monomials is None else monomials
    out = 0
    for row in terms:
        powers = tuple(sorted(row['powers'].items()))
        if powers not in monomials:
            value = 1
            for name, power in powers: value *= values[name] ** power
            monomials[powers] = value
        out += Q(row['coefficient']) * monomials[powers]
    return out


def evaluate_sources(m, s, b, names=None, path=None):
    """Return projected source fields with Fraction, float, or interval values.

    All derivatives are ordinary physical y/t derivatives in the fixed
    background orthonormal frame; named time sectors remain covariant. The
    QEred/Qphi actual mode assumes C1=0 as an identity on a neighbourhood.
    """
    report = load_maps(path)
    values = dict(m)
    if set(values) & set(s) or set(values) & set(b) or set(s) & set(b):
        raise ValueError('OVERLAPPING_INPUT_SYMBOLS')
    values.update(s); values.update(b)
    values = {name: Q(value) if isinstance(value, int) else value
              for name, value in values.items()}
    missing = set(report['input_symbols']) - set(values)
    if missing: raise ValueError('MISSING_INPUT_SYMBOLS_' + ','.join(sorted(missing)))
    names = list(report['maps']) if names is None else names
    cache = {}
    out = {}
    for name in names:
        mapped = report['maps'][name]
        if isinstance(mapped, dict):
            out[name] = {field: tuple(evaluate_encoded(row, values, cache) for row in coefficients)
                         for field, coefficients in mapped.items()}
        else:
            out[name] = tuple(evaluate_encoded(row, values, cache) for row in mapped)
    return out


def evaluate_bulk_forcing23(m, s, b, path=None):
    """Complete reduced Einstein plus scalar PDE forcing in registered order."""
    from angular_project import pack23
    sources = evaluate_sources(m, s, b, ('metric_forcing', 'scalar_forcing'), path)
    return pack23({**sources['metric_forcing'], 'scalar': sources['scalar_forcing']})
