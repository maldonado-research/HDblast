"""Deterministic receipt and source-integrity support for portable checkers."""
from pathlib import Path
import hashlib
import json


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_result(result, expected_checks, output, script):
    script = script.resolve()
    here = script.parent
    bundle_root = here.parent
    manifest_path = here / 'SOURCE_INPUTS.json'
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('schema') != 1:
        raise ValueError('Unsupported source manifest schema')
    entries = manifest.get('files')
    if not isinstance(entries, list) or not entries:
        raise ValueError('Source manifest must contain a nonempty files list')
    integrity = {}
    for entry in entries:
        relative = Path(entry['file'])
        if relative.is_absolute() or '..' in relative.parts:
            raise ValueError('Source manifest entries must stay within the bundle')
        target = bundle_root / relative
        actual = sha256(target) if target.is_file() else None
        integrity[relative.as_posix()] = {
            'expected_sha256': entry['sha256'],
            'actual_sha256': actual,
            'match': actual == entry['sha256'],
        }
    checks = result['checks']
    algebra_pass = (len(checks) == expected_checks
                    and all(value is True for value in checks.values()))
    source_pass = all(entry['match'] for entry in integrity.values())
    result.update({
        'expected_checks': expected_checks,
        'algebra_all_pass': algebra_pass,
        'source_integrity_all_pass': source_pass,
        'source_integrity': integrity,
        'source_manifest_sha256': sha256(manifest_path),
        'script_sha256': sha256(script),
        'support_sha256': sha256(Path(__file__)),
        'all_pass': algebra_pass and source_pass,
    })
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n')
    print(f"{result['n_pass']}/{result['n_checks']} exact checks passed; "
          f"source integrity {'passed' if source_pass else 'FAILED'}", flush=True)
    return result['all_pass']
