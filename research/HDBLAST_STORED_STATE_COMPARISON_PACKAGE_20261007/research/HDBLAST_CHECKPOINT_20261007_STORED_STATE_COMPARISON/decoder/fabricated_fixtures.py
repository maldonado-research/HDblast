"""Fabricated byte fixtures ONLY. Never imports/opens a historical input."""
from dataclasses import replace
import hashlib
import io
import struct
import zipfile
import zlib
from exact_binary80 import InputSpec, Limits, MemberSpec, STORAGE_LAYOUT


LIMITS = Limits(8 << 20, 1 << 20, 8 << 20, 32, 4096, 3, 65536, 64, 65536, 65536)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def slot(significand=1 << 63, exponent=16383, negative=False, padding=b"\x00" * 6):
    if len(padding) != 6:
        raise ValueError("fabricated padding length")
    return struct.pack("<QH", significand, exponent | (0x8000 if negative else 0)) + padding


def npy(payload, descr="<f16", shape=(1,), *, literal=None, newline=b"\n", version=b"\x01\x00"):
    text = (literal if literal is not None else repr({"descr": descr, "fortran_order": False,
                                                    "shape": shape})).encode("ascii")
    header = text + b" " * ((-10 - len(text) - len(newline)) % 64) + newline
    return b"\x93NUMPY" + version + struct.pack("<H", len(header)) + header + payload


def member(name, data, descr="<f16", shape=(1,), decode=True, compression=None):
    header_bytes = 10 + struct.unpack("<H", data[8:10])[0]
    return MemberSpec(name, sha(data), len(data), sha(data[:header_bytes]), header_bytes,
                      descr, shape, (1, 0), False, decode, compression,
                      None if compression is None else zlib.crc32(data))


def input_spec(data, members, kind):
    return InputSpec(kind, sha(data), len(data), STORAGE_LAYOUT, sha(b"fabricated producer"),
                     sha(b"fabricated manifest"), sha(b"fabricated producer layout evidence"), tuple(members))


def standalone(data, descr="<f16", shape=(1,), decode=True):
    return input_spec(data, (member("array.npy", data, descr, shape, decode),), "npy")


def archive(entries, *, compression=zipfile.ZIP_STORED, force_zip64=False):
    output = io.BytesIO()
    members = []
    with zipfile.ZipFile(output, "w", compression=compression) as z:
        for name, data, descr, shape, decode in entries:
            if force_zip64:
                info = zipfile.ZipInfo(name)
                info.compress_type = compression
                with z.open(info, "w", force_zip64=True) as destination:
                    destination.write(data)
            else:
                z.writestr(name, data)
            members.append(member(name, data, descr, shape, decode, compression))
    data = output.getvalue()
    return data, input_spec(data, members, "npz")


def repin_file(spec, data):
    return replace(spec, file_sha256=sha(data), file_bytes=len(data))
