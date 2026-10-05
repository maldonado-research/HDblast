"""Check the archive inventory, extract safely and replay the sealed exact proofs."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import stat
import subprocess
import sys
import zipfile


def require(value, message):
    if not value:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--packet', type=Path, required=True)
    parser.add_argument('--output-directory', type=Path, required=True)
    args = parser.parse_args()
    packet, output = args.packet.resolve(), args.output_directory.resolve()
    require(not output.exists(), 'Fresh archive replay directory required')
    require(output != packet and packet not in output.parents and output not in packet.parents,
            'Output must be outside and disjoint from the packet')
    inventory = json.loads((packet / 'ARCHIVE.json').read_text())
    name = inventory['checkpoint_name']
    require('/' not in name and name not in ('', '.', '..'), 'Unsafe checkpoint name')
    archive_name = inventory['archive_filename']
    require(PurePosixPath(archive_name).name == archive_name and '\\' not in archive_name, 'Unsafe archive name')
    archive = packet / archive_name
    require(not archive.is_symlink() and archive.is_file(), 'Archive must be a regular file')
    raw = archive.read_bytes()
    require(len(raw) == inventory['bytes'] and hashlib.sha256(raw).hexdigest() == inventory['sha256'],
            'Archive byte size or SHA256 differs')
    with zipfile.ZipFile(archive) as package:
        members = package.infolist()
        require(len(members) == inventory['members'] and len({item.filename for item in members}) == len(members),
                'Archive member count or uniqueness differs')
        require(len(members) <= 256 and sum(item.file_size for item in members) <= 20*1024*1024,
                'Archive exceeds declared extraction safety limit')
        for item in members:
            path = PurePosixPath(item.filename)
            require(not path.is_absolute() and path.parts[0] == name
                    and all(part not in ('', '.', '..') for part in item.filename.split('/'))
                    and '\\' not in item.filename and stat.S_ISREG(item.external_attr >> 16)
                    and not item.flag_bits & 1, 'Unsafe or encrypted archive member')
        require(package.testzip() is None, 'Archive CRC differs')
        output.mkdir(parents=True)
        package.extractall(output / 'extracted')
    checkpoint = output / 'extracted' / name
    replay = output / 'proof-runs'
    process = subprocess.run([sys.executable, '-B', str(checkpoint / 'execution/replay_checkpoint.py'),
                              '--checkpoint', str(checkpoint),
                              '--expected-manifest-sha256', inventory['manifest_sha256'],
                              '--output-directory', str(replay)], capture_output=True, timeout=300)
    (output / 'replay.stdout.log').write_bytes(process.stdout)
    (output / 'replay.stderr.log').write_bytes(process.stderr)
    require(process.returncode == 0, 'Extracted checkpoint replay failed; inspect retained log')
    result = json.loads((replay / 'REPLAY_RECEIPT.json').read_text())
    require(result['status'] == 'PASS_FRESH_EXACT_FINITE_BAND_REPLAY', 'Replay result is not PASS')
    receipt = {'status': 'PASS_CLEAN_ARCHIVE_FINITE_BAND_REPLAY', 'archive_sha256': inventory['sha256'],
               'manifest_sha256': inventory['manifest_sha256'], 'archive_members': len(members),
               'all_payload_bytes_unchanged': result['all_payload_bytes_unchanged'],
               'proof_replay': result, 'raw_archive_bytes': len(raw)}
    (output / 'FRESH_REPRODUCTION.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n')
    print(json.dumps({'status': receipt['status'], 'archive_members': len(members),
                      'verifier_executions': len(result['routes'])}))


if __name__ == '__main__':
    main()
