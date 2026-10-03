"""Strict, versioned comparison for complete active-source diagnostics.

This module never runs checkpoint code, reads an array, or contacts a service.
Only declared elapsed time, RSS, and the ZIP-membership manifest receipt are
variable. Complete raw diagnostics and these observations remain in artifacts.
The frozen normal/-O scientific validators still check every numerical schema.
"""
import copy
import hashlib
import json
import math
import re

POLICY = 'HDBLAST_ACTIVE_SOURCE_SCIENTIFIC_PROJECTION_V1'
PRIMARY_TOP = ('schema_version', 'route', 'status', 'records', 'old_metric_status',
    'fixed_cases', 'interval', 'registration_sha256', 'freeze_commit',
    'input_manifest_sha256', 'classification', 'failures', 'fatal_failures',
    'attribution_cases', 'precision_gap_max_exact_rational',
    'quadrature_control_gap_max_exact_rational', 'producer_sha256',
    'source_unchanged_during_execution', 'resources', 'arithmetic')
INDEPENDENT_TOP = ('schema_version', 'old_metric_status', 'fixed_cases',
    'input_manifest_sha256', 'resource_limits', 'route', 'status', 'provenance',
    'interval', 'midpoint', 'controls', 'canonical_control',
    'reduction_decimal_digits', 'source_phase_arithmetic', 'quadrature_constants',
    'high_precision_scope', 'units', 'momentum_target', 'rows',
    'serialization_max_error_by_dps', 'native_profile_serialization_max_error',
    'native_profile_serialization_scope', 'limits', 'resources', 'producer_sha256')
PROVENANCE = ('freeze_commit', 'registration_sha256', 'registered_files_checked',
    'input_manifest_sha256', 'full_frozen_verification', 'post_run_frozen_verification')
FROZEN_VERIFICATION = ('public_freeze_commit', 'registration_sha256',
    'freeze_receipt_sha256', 'frozen_files', 'package_manifest_sha256',
    'input_capsule_verification')
CAPSULE_VERIFICATION = ('capsules_verified', 'members_verified',
    'array_payload_values_interpreted', 'input_manifest_sha256')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def exact_keys(value, fields, label):
    require(type(value) is dict and set(value) == set(fields),
            'Unexpected or missing projection fields: ' + label)


def hexadecimal(value, length):
    return isinstance(value, str) and re.fullmatch('[0-9a-f]{%d}' % length, value) is not None


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False, allow_nan=False).encode('utf-8')


def project(result, route):
    require(route in {'primary', 'independent'}, 'Unknown projection route')
    exact_keys(result, PRIMARY_TOP if route == 'primary' else INDEPENDENT_TOP, route)
    require(type(result['schema_version']) is int and result['schema_version'] == 1,
            'Unsupported scientific diagnostic schema')
    require(type(result['fixed_cases']) is int and result['fixed_cases'] == 12,
            'Wrong complete diagnostic case count')
    require(result['old_metric_status'] == 'FAIL', 'Historical metric status changed')
    require(hexadecimal(result['producer_sha256'], 64)
            and hexadecimal(result['input_manifest_sha256'], 64),
            'Missing producer or input-manifest identity')
    elapsed_key = 'seconds' if route == 'primary' else 'elapsed_seconds'
    resources = result['resources']
    exact_keys(resources, (elapsed_key, 'peak_rss_kib', 'scope', 'limit_seconds',
                          'limit_peak_rss_kib'), route + ' resources')
    seconds = resources[elapsed_key]
    require(type(seconds) in (int, float) and math.isfinite(seconds) and 0 < seconds <= 900,
            'Diagnostic resource time failed')
    require(type(resources['peak_rss_kib']) is int and 0 < resources['peak_rss_kib'] <= 262144,
            'Diagnostic resource memory failed')
    require(type(resources['limit_seconds']) is int and resources['limit_seconds'] == 900
            and type(resources['limit_peak_rss_kib']) is int
            and resources['limit_peak_rss_kib'] == 262144,
            'Diagnostic declared resource limits changed')
    require(type(resources['scope']) is str and resources['scope'], 'Missing resource scope')
    observations = {'excluded_resource_values': {
        elapsed_key: seconds, 'peak_rss_kib': resources['peak_rss_kib']}}
    if route == 'primary':
        require(result['status'] == 'COMPLETED_SAVED_DATA_DIAGNOSTIC'
                and result['source_unchanged_during_execution'] is True,
                'Primary route is incomplete or changed')
        require(hexadecimal(result['registration_sha256'], 64)
                and hexadecimal(result['freeze_commit'], 40), 'Missing primary freeze identity')
    else:
        require(result['status'] == 'COMPLETED_DIAGNOSTIC_WITHOUT_RECLASSIFYING_OLD_FAIL',
                'Independent route is incomplete')
        exact_keys(result['resource_limits'], ('wall_seconds', 'peak_rss_kib'),
                   'independent resource_limits')
        require(type(result['resource_limits']['wall_seconds']) is int
                and result['resource_limits']['wall_seconds'] == 900
                and type(result['resource_limits']['peak_rss_kib']) is int
                and result['resource_limits']['peak_rss_kib'] == 262144,
                'Independent declared limits changed')
        provenance = result['provenance']
        exact_keys(provenance, PROVENANCE, 'independent provenance')
        frozen = provenance['full_frozen_verification']
        exact_keys(frozen, FROZEN_VERIFICATION, 'full frozen verification')
        capsule = frozen['input_capsule_verification']
        exact_keys(capsule, CAPSULE_VERIFICATION, 'input capsule verification')
        require(hexadecimal(provenance['freeze_commit'], 40)
                and frozen['public_freeze_commit'] == provenance['freeze_commit']
                and hexadecimal(provenance['registration_sha256'], 64)
                and frozen['registration_sha256'] == provenance['registration_sha256']
                and hexadecimal(frozen['freeze_receipt_sha256'], 64),
                'Independent freeze identity is absent or inconsistent')
        require(provenance['post_run_frozen_verification'] == 'PASS',
                'Independent final authentication failed')
        require(type(frozen['frozen_files']) is dict and frozen['frozen_files']
                and all(type(name) is str and hexadecimal(digest, 64)
                        for name, digest in frozen['frozen_files'].items()),
                'Malformed frozen inventory')
        require(type(provenance['registered_files_checked']) is int
                and provenance['registered_files_checked'] == len(frozen['frozen_files']),
                'Independent checked inventory count differs')
        require(type(capsule['capsules_verified']) is int and capsule['capsules_verified'] == 4
                and type(capsule['members_verified']) is int and capsule['members_verified'] == 76
                and capsule['array_payload_values_interpreted'] is False,
                'Incomplete authenticated capsule inventory')
        require(capsule['input_manifest_sha256'] == result['input_manifest_sha256']
                == provenance['input_manifest_sha256']
                == frozen['frozen_files'].get('inputs/INPUT_MANIFEST.json'),
                'Independent input identity differs')
        package_pin = frozen['package_manifest_sha256']
        require(package_pin is None or hexadecimal(package_pin, 64),
                'Malformed package manifest receipt')
        observations['excluded_package_manifest_sha256'] = package_pin
    retained = copy.deepcopy(result)
    del retained['resources'][elapsed_key]
    del retained['resources']['peak_rss_kib']
    if route == 'independent':
        del retained['provenance']['full_frozen_verification']['package_manifest_sha256']
    # No recursive dropping of fields, normalization of numeric strings,
    # reordering of scientific lists, or tolerance-based comparison is allowed.
    payload = {'scientific_hash_policy': POLICY, 'route': route, 'diagnostic': retained}
    return hashlib.sha256(canonical_bytes(payload)).hexdigest(), observations
