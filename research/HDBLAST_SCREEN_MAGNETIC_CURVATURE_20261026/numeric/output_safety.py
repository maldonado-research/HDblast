"""Fresh external outputs and deterministic source-integrity guards."""
import hashlib
from pathlib import Path


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def packaged_ancestors(path):
    """Containing research/<capsule> roots, using resolved real paths."""
    path = Path(path).resolve()
    return [parent for parent in (path,*path.parents) if parent.parent.name == 'research']


def guard_output(output, *, source_file, protected_roots=(), input_paths=()):
    """Reject overwrites and all writes to source/frozen capsule trees."""
    output = Path(output)
    need(not output.is_symlink(),'OUTPUT_SYMLINK')
    need(not output.exists(),'OUTPUT_EXISTS')
    target = output.resolve()
    source = Path(source_file).resolve()
    roots = [('SOURCE_TREE',source.parent)]
    roots.extend(('PACKAGED_CAPSULE',root) for root in packaged_ancestors(source))
    roots.extend(protected_roots)
    for input_path in input_paths:
        if input_path is not None:
            roots.extend(('INPUT_CAPSULE',root) for root in packaged_ancestors(input_path))
    for reason,root in roots:
        root = Path(root).resolve()
        need(target != root and root not in target.parents,'OUTPUT_INSIDE_'+reason)
    need(target.parent.is_dir(),'OUTPUT_PARENT_MISSING')
    return output


def source_snapshot(paths):
    return {str(Path(path).resolve()):sha(path) for path in paths}


def verify_source_snapshot(before):
    for path,pin in before.items():
        need(sha(path) == pin,'SOURCE_CHANGED:'+Path(path).name)


def source_labels(before):
    """Stable canonical report labels; absolute paths remain local only."""
    labels = {Path(path).name:pin for path,pin in before.items()}
    need(len(labels) == len(before),'SOURCE_LABEL_COLLISION')
    return labels


def require_hash(path, expected, reason):
    need(sha(path) == expected,reason)
