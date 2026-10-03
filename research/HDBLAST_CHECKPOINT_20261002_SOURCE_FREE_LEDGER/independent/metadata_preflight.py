"""Read only raw hashes, CRCs and NPY headers. Never calls np.load or a profile."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
from diagnostic_independent import input_metadata, runtime_check, sha
from exact_binary import require


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--inputs', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), 'Metadata receipt must be fresh')
    runtime = runtime_check()
    manifest, capsules = input_metadata(args.inputs)
    receipt = {'status': 'PASS_METADATA_ONLY', 'physical_array_payloads_decoded': False,
               'physical_source_or_mode_functions_called': False,
               'input_manifest_sha256': sha(args.inputs / 'INPUT_MANIFEST.json'),
               'capsule_files': [{'path': path.name, 'bytes': path.stat().st_size,
                                  'sha256': sha(path)} for path in capsules.values()],
               'raw_NPY_hash_and_header_verified': sum(len(c['members']) for c in manifest['capsules']),
               'runtime': runtime}
    args.output.write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps(receipt, indent=2))
