"""Bounded-memory standard NPZ output; no scientific values are defined here."""
from __future__ import annotations

import hashlib
from pathlib import Path
import zipfile

import numpy as np


def file_sha256(path, chunk_bytes=1024*1024):
    if chunk_bytes <= 0:
        raise ValueError("A positive bounded hash chunk is required")
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(chunk_bytes), b""):
            digest.update(chunk)
    return digest.hexdigest()


class StreamingNpz:
    """Write each exact ndarray immediately without retaining its reference.

    Each member is the same .npy stream used by numpy.savez_compressed. The
    archive remains a standard ZIP_DEFLATED/ZIP64 NPZ, read by numpy.load.
    Only member names and small ZIP directory metadata remain in memory.
    """

    def __init__(self, path):
        self.path = Path(path)
        self._zip = zipfile.ZipFile(self.path, mode="x", compression=zipfile.ZIP_DEFLATED,
                                    allowZip64=True)
        self._names = set()
        self._closed = False

    def __setitem__(self, key, value):
        if self._closed:
            raise RuntimeError("Cannot append to a closed NPZ archive")
        if not isinstance(key, str) or not key or key in self._names or "/" in key or "\\" in key:
            raise ValueError("NPZ member must have a unique nonempty flat name")
        array = np.asanyarray(value)
        if array.dtype.hasobject:
            raise ValueError("Object arrays are excluded from physical raw archives")
        with self._zip.open(key+".npy", mode="w", force_zip64=True) as member:
            np.lib.format.write_array(member, array, allow_pickle=False)
        self._names.add(key)

    def update(self, values):
        for key, value in values.items():
            self[key] = value

    def close(self):
        if not self._closed:
            self._zip.close()
            self._closed = True

    def __enter__(self):
        return self

    def __exit__(self, exception_type, exception, traceback):
        self.close()
        return False
