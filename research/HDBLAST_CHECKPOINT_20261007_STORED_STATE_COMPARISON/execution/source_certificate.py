"""Build/export fixed source models using the three immutable prior files.

There are no source constructions at module import. Actual construction is
reachable only with a new authenticated SourceAuthorization callback.
"""
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
import sys
from registration_guard import SOURCES, require, safe_file, load, rational, sha

def canonical(q):
    require(type(q) is Q, 'Exact real rational coefficient required')
    return f'{q.numerator}/{q.denominator}'

def coefficient_sha(coefficients):
    return sha(json.dumps([canonical(c) for c in coefficients], separators=(',', ':')).encode())

def import_prior(root):
    pins = load(safe_file(root, 'SOURCE_PINS.json').read_bytes())
    require(type(pins) is dict and set(pins['files']) == {
        'source/route.py','source/source_algebra_baseline.py','source/generic_operator_baseline.py',
        'prior/BUDGET.json','prior/primary_DATA.json','prior/independent_DATA.json'}, 'Exact prior inventory required')
    for name, item in pins['files'].items():
        raw = safe_file(root, name).read_bytes()
        require(len(raw) == item['bytes'] and sha(raw) == item['sha256'], 'Immutable prior pin differs: '+name)
    source = root/'source'
    sys.path.insert(0, str(source))
    spec = importlib.util.spec_from_file_location('unchanged_prior_independent_route', source/'route.py')
    require(spec is not None and spec.loader is not None, 'Prior builder import failed')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    require(module.SOURCES == SOURCES and module.SOURCE_DEGREE == 112 and module.BITS == 512 and
            module.EXP_DEGREE == 200 and module.DELTA == Q(1,128), 'Unchanged builder constants differ')
    return module

def geometry_proof(row, source, parent, prior, sid):
    expected = {'source':source, 'parent':parent, 'center':canonical(row['center']),
                'interval':[canonical(row['left']), canonical(row['right'])],
                'half_width':canonical(row['half']), 'analytic_radius':canonical(row['radius']),
                'disk_bound':canonical(row['M'][sid]), 'children':row['children']}
    require(set(prior) == set(expected)|{'chosen_coefficients_sha256','analytic_tail','coefficient_error',
                                      'uniform_source_error','exp_squarings'}, 'Prior source proof schema differs')
    for key, val in expected.items():
        require(type(prior[key]) is type(val) and prior[key] == val, 'Prior source geometry/index differs: '+key)
    require(row['right']-row['left'] == 2*row['half'] and row['center'] == (row['left']+row['right'])/2,
            'Internal source geometry inconsistent')
    require(rational(prior['analytic_tail']) == row['M'][sid]/2**112, 'Prior Cauchy source tail differs')
    for key in ('analytic_tail','coefficient_error','uniform_source_error'):
        require(rational(prior[key]) >= 0, 'Negative prior source error')
    require(type(prior['exp_squarings']) is int and prior['exp_squarings'] >= 0, 'Invalid prior exponential proof')
    return expected

def build_models(root, prior_route, authorization=None, fabricated=False):
    geom = prior_route.panels()
    require(len(geom) == 11, 'Eleven fixed source parents required')
    budget = load(safe_file(root,'prior/BUDGET.json').read_bytes())
    prior_rows = budget['source_model_rows']
    require(budget['fabricated_source_provider'] is False and budget['source_callbacks'] == 22 and
            budget['source_degree'] == 112 and budget['coefficient_bits'] == 512 and len(prior_rows) == 22,
            'Prior independent successful budget identity differs')
    models, certificate = {}, []
    for sid, source in enumerate(SOURCES):
        rows = []
        for parent, row in enumerate(geom):
            prior = prior_rows[sid*11+parent]
            proof = geometry_proof(row, source, parent, prior, sid)
            if fabricated:
                coeff, error, ev = prior_route.fabricated_model(row, sid)
                uniform = error
            else:
                require(authorization is not None and callable(authorization), 'New source authorization required')
                coeff, error, ev = prior_route.physical_model(row, source, authorization)
                require(coefficient_sha(coeff) == prior['chosen_coefficients_sha256'],
                        'Fresh source coefficients differ from prior independent proof')
                require(set(ev) == {'analytic_tail','coefficient_error','uniform_source_error','exp_squarings'},
                        'Fresh source proof schema differs')
                for key in ev:
                    require(type(ev[key]) is type(prior[key]) and ev[key] == prior[key],
                            'Fresh source proof differs: '+key)
                require(canonical(error) == ev['uniform_source_error'], 'Rebuilt model error inconsistent')
                uniform = max(error, rational(prior['uniform_source_error']),
                              rational(prior['analytic_tail'])+rational(prior['coefficient_error']))
            require(type(coeff) is tuple and len(coeff) == 113 and all(type(c) is Q for c in coeff) and
                    type(uniform) is Q and uniform >= 0, 'Exactly 113 certified exact real coefficients required')
            rows.append({'center':row['center'],'half':row['half'],'coefficients':coeff,'uniform_error':uniform})
            certificate.append({**proof, 'coefficients':[canonical(c) for c in coeff],
                                'coefficient_sha256':coefficient_sha(coeff),
                                'prior_coefficient_sha256':prior['chosen_coefficients_sha256'],
                                'rebuilt_proof':ev, 'uniform_error':canonical(uniform),
                                'prior_proof':prior, 'fabricated':fabricated})
        models[source] = rows
    return models, {'schema_version':1,'scope':'FABRICATED_FINITE_POLYNOMIAL' if fabricated else 'CERTIFIED_UNCHANGED_SOURCE_MODELS',
                    'source_count':2,'parent_count_per_source':11,'vector_count':22,'coefficients_per_vector':113,
                    'uniform_error_rule':'max(rebuilt_error,prior_uniform_error,prior_analytic_tail+prior_coefficient_error)',
                    'rows':certificate}

def validate_certificate(certificate, root):
    require(type(certificate) is dict and set(certificate)=={'schema_version','scope','source_count',
            'parent_count_per_source','vector_count','coefficients_per_vector','uniform_error_rule','rows'} and
            type(certificate['schema_version']) is int and certificate['schema_version']==1 and
            type(certificate['source_count']) is int and certificate['source_count']==2 and
            type(certificate['parent_count_per_source']) is int and certificate['parent_count_per_source']==11 and
            type(certificate['vector_count']) is int and certificate['vector_count'] == 22 and
            certificate['coefficients_per_vector'] == 113 and len(certificate['rows']) == 22,
            'Full standalone source certificate required')
    require(certificate['uniform_error_rule']=='max(rebuilt_error,prior_uniform_error,prior_analytic_tail+prior_coefficient_error)',
            'Certificate source error rule differs')
    prior_budget = load(safe_file(root,'prior/BUDGET.json').read_bytes())
    prior_rows = prior_budget['source_model_rows']
    require(type(certificate['rows'][0]['fabricated']) is bool, 'Source certificate mode missing')
    fabricated=certificate['rows'][0]['fabricated']
    require(certificate['scope']==('FABRICATED_FINITE_POLYNOMIAL' if fabricated else 'CERTIFIED_UNCHANGED_SOURCE_MODELS'),
            'Certificate source scope differs')
    for i, row in enumerate(certificate['rows']):
        require(type(row) is dict and set(row)=={'source','parent','center','interval','half_width','analytic_radius',
                'disk_bound','children','coefficients','coefficient_sha256','prior_coefficient_sha256',
                'rebuilt_proof','uniform_error','prior_proof','fabricated'} and row['fabricated'] is fabricated,
                'Consistent source certificate row schema/mode required')
        require(row['source'] == SOURCES[i//11] and row['parent'] == i%11, 'Certificate source order changed')
        require(row['prior_proof'] == prior_rows[i], 'Certificate prior source proof changed')
        for field in ('source','parent','center','interval','half_width','analytic_radius','disk_bound','children'):
            require(type(row[field]) is type(prior_rows[i][field]) and row[field] == prior_rows[i][field],
                    'Certificate source geometry/proof alias differs: '+field)
        co = tuple(rational(c) for c in row['coefficients'])
        require(len(co) == 113 and coefficient_sha(co) == row['coefficient_sha256'], 'Certificate vector changed')
        require(row['prior_coefficient_sha256'] == prior_rows[i]['chosen_coefficients_sha256'], 'Certificate prior coefficient pin changed')
        if not row['fabricated']:
            require(row['coefficient_sha256'] == row['prior_coefficient_sha256'], 'Certificate fresh/prior vector differs')
            require(row['rebuilt_proof'] == {k:prior_rows[i][k] for k in
                    ('analytic_tail','coefficient_error','uniform_source_error','exp_squarings')}, 'Certificate proof changed')
            expected = max(rational(row['rebuilt_proof']['uniform_source_error']),
                           rational(prior_rows[i]['uniform_source_error']),
                           rational(prior_rows[i]['analytic_tail'])+rational(prior_rows[i]['coefficient_error']))
            require(rational(row['uniform_error']) == expected, 'Certificate uniform error weakened')
        else:
            require(rational(row['uniform_error']) == 0, 'Fabricated finite model error differs')
    return {'vectors_verified':22,'exact_real_coefficients_verified':2486,'physical_source_constructed':False}
