#!/usr/bin/env bash
# Complete registered replay; source and frozen outputs remain unchanged.
set -Eeuo pipefail
source_root=${1:?Usage: replay_smooth_frw.sh CHECKPOINT_ROOT FRESH_WORK_DIRECTORY}
work=${2:?Usage: replay_smooth_frw.sh CHECKPOINT_ROOT FRESH_WORK_DIRECTORY}
source_root=$(cd "$source_root" && pwd)
python=${HDBLAST_PYTHON:-python3}
work=$("$python" -c 'from pathlib import Path; import sys; print(Path(sys.argv[1]).resolve())' "$work")
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONOPTIMIZE=0 MPLBACKEND=Agg
if [[ -e "$work" ]]; then
    printf 'Replay directory already exists: %s\n' "$work" >&2
    exit 2
fi
case "$work/" in "$source_root/"*) printf 'Replay directory must be outside frozen checkpoint\n' >&2; exit 2;; esac
mkdir -p "$work/logs" "$work/runtime"
export MPLCONFIGDIR="$work/matplotlib-cache"
first_failure=0
recorded_exit=0
record_command() {
    local label=$1
    shift
    local rc
    printf 'Running %s\n' "$label"
    if "$@" > "$work/logs/$label.log" 2>&1; then rc=0; else rc=$?; fi
    printf '%s\t%s\n' "$label" "$rc" >> "$work/exit-codes.tsv"
    recorded_exit=$rc
    printf '%s exit: %s; log: %s\n' "$label" "$rc" "$work/logs/$label.log"
}
run_check() {
    record_command "$@"
    if [[ "$recorded_exit" -ne 0 && "$first_failure" -eq 0 ]]; then first_failure=$recorded_exit; fi
}

# Verify frozen bytes BEFORE copying or running any scientific output writer.
run_check frozen-package "$python" "$source_root/code/verify_package.py" "$source_root/$(basename "$source_root").zip"
if [[ "$first_failure" -ne 0 ]]; then exit "$first_failure"; fi
mkdir "$work/source"
source_copy="$work/source/$(basename "$source_root")"
cp -a "$source_root" "$source_copy"
code="$source_copy/code"
protocol="$source_copy/REGISTRATION.md"
protocol_hash=37fda5b30332d0079ec44db708e266788fdf5859303707bb709750fdcf59d898
run_check package-copy "$python" "$code/verify_package.py" "$source_copy/$(basename "$source_root").zip"
run_check dependencies "$python" -m pip check
run_check runtime-versions "$python" - <<'PY'
import json
import platform
import sys
import numpy
import scipy
import matplotlib
import sympy
import mpmath
assert sys.flags.optimize == 0
assert sys.version_info[:2] == (3, 12)
versions = {name: module.__version__ for name, module in [('numpy', numpy), ('scipy', scipy), ('matplotlib', matplotlib), ('sympy', sympy), ('mpmath', mpmath)]}
assert versions == dict(numpy='2.2.6', scipy='1.15.3', matplotlib='3.10.1', sympy='1.14.0', mpmath='1.3.0'), versions
assert numpy.finfo(numpy.longdouble).eps < 1e-18
print(json.dumps(dict(python=platform.python_version(), packages=versions, longdouble_epsilon=float(numpy.finfo(numpy.longdouble).eps))))
PY
run_check derive-local-trace "$python" "$code/derive_local_trace.py"
run_check verify-local-identities "$python" "$code/verify_local_identities.py"
run_check verify-trace-integrand "$python" "$code/verify_trace_integrand.py"
run_check numeric-preflight "$python" "$code/verify_numeric_counterterms.py" --output "$work/runtime/preflight.json"
if [[ "$first_failure" -ne 0 ]]; then exit "$first_failure"; fi

# Preserve the actual original command exit separately. Only the exact historical
# scientific failure qualifies as prerequisite evidence for the new registration.
# The original summary remains FAIL; every unexpected failure is still fatal.
record_command original-matrix "$python" -u "$code/smooth_frw_control.py" --execute --protocol "$protocol" --protocol-sha256 "$protocol_hash" --output "$work/runtime/primary"
original_exit=$recorded_exit
printf 'Original matrix actual exit: %s; scientific outcome preserved in summary.json\n' "$original_exit"
run_check independent-radau "$python" -u "$code/independent_radau.py" --protocol "$protocol" --protocol-sha256 "$protocol_hash" --output "$work/runtime/independent.json"
run_check independent-comparison "$python" -u "$code/compare_independent_modes.py" --protocol "$protocol" --protocol-sha256 "$protocol_hash" --primary "$code/smooth_frw_control.py" --independent "$work/runtime/independent.json" --output "$work/runtime/comparison.json"
run_check plotting "$python" "$code/plot_smooth_frw.py" "$work/runtime/primary"
run_check original-outcome-validation "$python" "$code/validate_original_outcome.py" --run-root "$work/runtime/primary" --matrix-exit-code "$original_exit" --output "$work/runtime/original-outcome.json"
if [[ "$recorded_exit" -ne 0 && "$original_exit" -ne 0 ]]; then
    printf 'Original outcome was unexpected; retaining original command exit %s\n' "$original_exit" >&2
    exit "$original_exit"
fi
if [[ "$first_failure" -ne 0 ]]; then exit "$first_failure"; fi

followup_hash=b067afda194895255081553e7bd9a3833b781f5277d3011b619829b4b3033d51
# The new experiment consumes the prospectively pinned historical baseline,
# whose bytes remain distinct from a fresh original replay's metadata.
run_check registered-followup "$python" -u "$code/smooth_frw_followup.py" --execute --core "$code/smooth_frw_control.py" --baseline "$source_copy/followup_baseline" --protocol "$source_copy/FOLLOWUP_REGISTRATION.md" --protocol-sha256 "$followup_hash" --output "$work/runtime/followup"
run_check followup-plotting "$python" "$code/plot_smooth_frw_followup.py" "$work/runtime/followup"
run_check nonempty-terminal-gates "$python" - "$work" <<'PY'
import json
from pathlib import Path
import sys

work = Path(sys.argv[1])
def read(path):
    return json.loads(path.read_text(), parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
summary = read(work / 'runtime/primary/summary.json')
assert summary['status'] == 'FAIL', summary['status']
assert set(summary['amplitudes']) == {'0.0', '0.2'}
assert summary['amplitudes']['0.0']['status'] == 'PASS' and summary['amplitudes']['0.2']['status'] == 'FAIL'
assert summary['static_negative_control']['passed'] is True
outcome = read(work / 'runtime/original-outcome.json')
assert outcome['scientific_status'] == 'FAIL' and outcome['original_command_exit_code'] == 1
assert outcome['validation_status'] == 'EXPECTED_SCIENTIFIC_FAILURE_CONFIRMED' and outcome['all_other_acceptance_gates_passed'] is True
followup = read(work / 'runtime/followup/summary.json')
assert followup['status'] == 'PASS' and followup['original_matrix_status'] == 'FAIL'
assert followup['cutoff'] == 384 and followup['amplitude'] == .2
expected_followup_gates = {'cutoff', 'solver', 'quadrature', 'subset_continuity', 'time', 'wronskian', 'integrated_exchange', 'local_exchange', 'independent_ODE_work', 'direct_trace', 'sampled_trace', 'future_spectral_energy', 'static_negative_control'}
assert set(followup['gates']) == expected_followup_gates
assert all(value is True for value in followup['gates'].values())
assert followup['provenance']['followup_protocol_sha256'] == 'b067afda194895255081553e7bd9a3833b781f5277d3011b619829b4b3033d51'
assert followup['provenance']['core_sha256'] == '32f6632bb84cc10dd6cef22c2311941acdbe223d846eb7157a2f1e9459fb160a'
assert len(list((work / 'runtime/followup').glob('*/diagnostics.json'))) == 4
independent = read(work / 'runtime/independent.json')
assert len(independent['cases']) == 12
assert all(len(case['samples']) == 10 for case in independent['cases'])
assert len(independent['gates']) == 6 and all(value is True for value in independent['gates'].values())
comparison = read(work / 'runtime/comparison.json')
assert comparison['status'] == 'PASS' and len(comparison['cases']) == 12
assert all(len(case['samples']) == 10 for case in comparison['cases'])
assert sum(len(row['gates']) for case in comparison['cases'] for row in case['samples']) == 600
assert all(gate['passed'] is True for case in comparison['cases'] for row in case['samples'] for gate in row['gates'].values())
for name in ('verify-local-identities', 'verify-trace-integrand'):
    assert read(work / ('logs/' + name + '.log'))['status'] == 'PASS'
assert read(work / 'runtime/preflight.json')['status'] == 'PASS'
assert (work / 'logs/derive-local-trace.log').read_text().count('PASS:') == 3
print('PASS: original scientific FAIL independently confirmed and preserved; separate K384 follow-up has 13 passing gates; 12 Radau cases, 600 comparisons, three symbolic scripts, preflight')
PY
printf 'Replay evidence: %s; final exit: %s\n' "$work" "$first_failure"
exit "$first_failure"
