"""Pack/extract the exact compressed source tables with deterministic ZIP metadata."""
import argparse
from pathlib import Path
import zipfile


def pack(directory, output):
    directory, output = Path(directory), Path(output)
    paths = [directory/'SOURCES.json', *sorted((directory/'tables').glob('*.npz'))]
    if len(paths) != 33:
        raise ValueError('Expected one manifest and 32 source tables')
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_STORED) as archive:
        for path in paths:
            entry = zipfile.ZipInfo(path.relative_to(directory).as_posix(),(1980,1,1,0,0,0))
            entry.external_attr = 0o100644 << 16
            entry.create_system = 3
            archive.writestr(entry,path.read_bytes())


def extract(archive_path, output):
    output = Path(output)
    with zipfile.ZipFile(archive_path) as archive:
        for name in archive.namelist():
            path = Path(name)
            if path.is_absolute() or '..' in path.parts:
                raise ValueError('Unsafe archive member')
        archive.extractall(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command',required=True)
    packing = commands.add_parser('pack')
    packing.add_argument('--directory',type=Path,required=True)
    packing.add_argument('--output',type=Path,required=True)
    extracting = commands.add_parser('extract')
    extracting.add_argument('--archive',type=Path,required=True)
    extracting.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    if args.command == 'pack':
        pack(args.directory,args.output)
    else:
        extract(args.archive,args.output)
