"""Install an authenticated new capsule and dedicated CI while preserving all old workflows."""
import argparse
from hashlib import sha256
import importlib.util
import io
import json
from pathlib import Path
import re
import shutil
import subprocess
import zipfile

CAPSULE = 'HDBLAST_COORDINATE_SHELL_SOURCES_20261015'
OLD_CAPSULE = 'HDBLAST_QUADRATIC_SHELL_FORCING_20261015'
WORKFLOW = '.github/workflows/hdblast_coordinate_shell_sources.yml'


def require(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def main():
    task = Path(__file__).resolve().parent.parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worktree', type=Path, required=True)
    parser.add_argument('--frozen-dir', type=Path, default=task/'frozen-source')
    parser.add_argument('--runner-config', type=Path, default=task/'publication/RUNNER_CONFIG.json')
    parser.add_argument('--expected-existing-workflows', type=int, default=44)
    parser.add_argument('--output', type=Path, default=task/'publication/CI_INSTALL_PREFLIGHT.json')
    args = parser.parse_args()
    tree, frozen = args.worktree.absolute(), args.frozen_dir.absolute()
    require(tree.is_dir() and not any(p.is_symlink() for p in [tree, *tree.parents]), 'Real publication worktree required')
    receipt = json.loads((frozen/'SOURCE_FREEZE_RECEIPT.json').read_bytes())
    require(receipt['status'] == 'PASS' and receipt['capsule'] == CAPSULE, 'Matching completed source freeze required')
    source, target = frozen/CAPSULE, tree/'research'/CAPSULE
    archive = frozen/(CAPSULE+'.zip')
    for path in (target, target.with_suffix('.zip'), tree/WORKFLOW, args.output.absolute()):
        require(not any(p.is_symlink() for p in [path, *path.parents]), 'Real installation and receipt paths required')
    require(not args.output.exists() and not args.output.resolve().is_relative_to(tree.resolve())
            and not args.output.resolve().is_relative_to(frozen.resolve()),
            'New CI receipt must remain outside the public worktree and frozen source')
    require(not target.exists() and not target.with_suffix('.zip').exists() and not (tree/WORKFLOW).exists(),
            'New capsule/archive/workflow paths required')
    manifest_raw = (source/'MANIFEST.json').read_bytes()
    archive_raw = archive.read_bytes()
    require(len(manifest_raw) == receipt['manifest']['bytes'] and sha256(manifest_raw).hexdigest() == receipt['manifest']['sha256']
            and len(archive_raw) == receipt['archive']['bytes'] and sha256(archive_raw).hexdigest() == receipt['archive']['sha256'],
            'External freeze bytes differ')
    spec = importlib.util.spec_from_file_location('source_verifier', task/'publication/replay-tools/verify_manifest.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    # Use the audited publication helper, without importing/executing capsule code.
    verified = module.verify(source, receipt['manifest']['sha256'])
    require(verified['files_verified'] == receipt['source_files'], 'Frozen source count differs')
    require((source/'verify_manifest.py').read_bytes() == (task/'publication/replay-tools/verify_manifest.py').read_bytes(),
            'Frozen verifier must equal the reviewed audited support tool')
    expected_zip = {CAPSULE+'/'+p.relative_to(source).as_posix(): p.read_bytes()
                    for p in source.rglob('*') if p.is_file()}
    with zipfile.ZipFile(io.BytesIO(archive_raw)) as zipped:
        names = zipped.namelist()
        require(len(names) == len(set(names)) == receipt['archive_members']
                and set(names) == set(expected_zip), 'Complete archive inventory differs')
        for name, raw in expected_zip.items():
            require(zipped.read(name) == raw, 'Archive/source bytes differ')
    def git(*args):
        return subprocess.check_output(['git', '-C', str(tree), *args])
    old_names = git('ls-tree', '-r', '--name-only', 'HEAD', '--', '.github/workflows').decode().splitlines()
    require(len(old_names) == args.expected_existing_workflows, 'Existing workflow count differs; inspect baseline')
    local_names = {(p.relative_to(tree).as_posix()) for p in (tree/'.github/workflows').rglob('*') if p.is_file()}
    require(local_names == set(old_names), 'Unreviewed existing workflow paths require inspection')
    for name in old_names:
        require((tree/name).read_bytes() == git('show', 'HEAD:'+name), 'Existing workflow changed: '+name)
    config_raw = args.runner_config.read_bytes()
    config = json.loads(config_raw)
    runner = (source/'reproduce.py').read_bytes()
    require(config['capsule'] == CAPSULE and config['status'] == 'APPROVED_FOR_FROZEN_REPLAY'
            and ("CONFIG_SHA256 = '"+sha256(config_raw).hexdigest()+"'").encode() in runner,
            'CI Python pin must come from embedded approved runner config')
    template_raw = (task/'publication/workflow_template.yml').read_bytes()
    provenance = json.loads((task/'publication/TOOL_REUSE_PROVENANCE.json').read_bytes())
    template_pin = [r for r in provenance['files'] if r['path'] == 'workflow_template.yml']
    require(len(template_pin) == 1 and sha256(template_raw).hexdigest() == template_pin[0]['prior_sha256'],
            'Audited workflow template bytes changed')
    text = template_raw.decode().replace(OLD_CAPSULE, CAPSULE).replace(
        'Quadratic material shell forcing replay', 'Coordinate shell source coefficients replay').replace(
        'hdblast_quadratic_shell_forcing.yml', Path(WORKFLOW).name).replace(
        'quadratic-shell-forcing', 'coordinate-shell-sources')
    text, n = re.subn(r"(python-version: ')[^']+(')", lambda m: m[1]+config['python_version']+m[2], text)
    require(n == 1, 'Unique exact workflow Python version required')
    pins = {'MANIFEST_SHA256': receipt['manifest']['sha256'], 'ARCHIVE_SHA256': receipt['archive']['sha256'],
            'RUNNER_SHA256': sha256(runner).hexdigest(),
            'VERIFIER_SHA256': sha256((source/'verify_manifest.py').read_bytes()).hexdigest()}
    for key, value in pins.items():
        text, n = re.subn(r"("+key+r": ')[0-9a-f]{64}(')", lambda m: m[1]+value+m[2], text)
        require(n == 1, 'Unique workflow pin required: '+key)
    shutil.copytree(source, target)
    target.with_suffix('.zip').write_bytes(archive_raw)
    (tree/WORKFLOW).write_text(text)
    for name in old_names:
        require((tree/name).read_bytes() == git('show', 'HEAD:'+name), 'Existing workflow changed during installation')
    for path in source.rglob('*'):
        if path.is_file():
            require(path.read_bytes() == (target/path.relative_to(source)).read_bytes(), 'Installed source bytes differ')
    result = {'status': 'PASS', 'capsule': CAPSULE, 'source_files': receipt['source_files'],
              'archive_members': receipt['archive_members'], 'old_workflows_preserved': len(old_names),
              'workflow_path': WORKFLOW, 'workflow_sha256': sha256(text.encode()).hexdigest(), 'pins': pins,
              'template_sha256': sha256(template_raw).hexdigest(), 'python_version': config['python_version'],
              'scope': 'Local capsule/CI installation and preservation only; no Git action, remote CI or publication.'}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
