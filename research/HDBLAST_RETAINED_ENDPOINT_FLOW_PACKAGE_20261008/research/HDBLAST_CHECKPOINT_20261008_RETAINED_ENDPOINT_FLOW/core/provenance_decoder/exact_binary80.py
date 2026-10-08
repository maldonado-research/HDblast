"""Strict, pinned, stdlib-only NPY/NPZ reader for one known producer layout.

No NumPy, float conversion, pickle, source callbacks, or implicit dtype inference.
Verification snapshots bytes and verifies the file and ALL members before exposing
an iterator. Only canonical finite x87 binary80 values are accepted. See README.
"""
from __future__ import annotations

import ast
from dataclasses import asdict, dataclass
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import stat
import struct
import tempfile
import threading
from typing import BinaryIO, Iterator
import zipfile
import zlib

STORAGE_LAYOUT = "x87-binary80-little-endian-padded16-v1"
_SHA = re.compile(r"[0-9a-f]{64}\Z")
_NAME = re.compile(r"[A-Za-z0-9_]+\.npy\Z")
_ITEMSIZE = {"<f16": 16, "<c32": 32, "<i8": 8, "<U8": 32}


class CodecError(ValueError):
    """Rejected format, pin, limit, or binary80 encoding."""


@dataclass(frozen=True)
class Limits:
    max_file_bytes: int
    max_member_bytes: int
    max_total_member_bytes: int
    max_members: int
    max_header_bytes: int
    max_dimensions: int
    max_elements: int
    max_decompression_ratio: int
    chunk_bytes: int
    spool_memory_bytes: int


@dataclass(frozen=True)
class MemberSpec:
    name: str
    sha256: str
    nbytes: int
    header_sha256: str
    header_bytes: int
    descr: str
    shape: tuple[int, ...]
    npy_version: tuple[int, int]
    fortran_order: bool
    decode: bool
    compression: int | None
    crc32: int | None


@dataclass(frozen=True)
class InputSpec:
    kind: str
    file_sha256: str
    file_bytes: int
    storage_layout: str
    producer_sha256: str
    manifest_sha256: str
    layout_evidence_sha256: str
    members: tuple[MemberSpec, ...]


@dataclass(frozen=True)
class ExactReal:
    value: Fraction
    negative_zero: bool


@dataclass(frozen=True)
class ExactComplex:
    real: ExactReal
    imag: ExactReal


@dataclass(frozen=True)
class CodecCheckpoint:
    input_spec_sha256: str
    file_sha256: str
    member_name: str
    member_sha256: str
    next_index: int


def _integer(value: object, minimum: int, maximum: int, label: str) -> int:
    if type(value) is not int or not minimum <= value <= maximum:
        raise CodecError(f"invalid {label}")
    return value


def _digest(value: object, label: str) -> None:
    if type(value) is not str or not _SHA.fullmatch(value):
        raise CodecError(f"invalid {label} SHA256")


def _validate_spec(spec: InputSpec, limits: Limits) -> str:
    if type(spec) is not InputSpec or type(limits) is not Limits:
        raise CodecError("explicit InputSpec and Limits required")
    # Hard ceilings bound parser/snapshot resource use even for mistaken caps.
    ceilings = {"max_file_bytes": 128 << 20, "max_member_bytes": 32 << 20,
                "max_total_member_bytes": 128 << 20, "max_members": 256,
                "max_header_bytes": 16384, "max_dimensions": 8,
                "max_elements": 1 << 20, "max_decompression_ratio": 1024,
                "chunk_bytes": 1 << 20, "spool_memory_bytes": 8 << 20}
    for name, ceiling in ceilings.items():
        _integer(getattr(limits, name), 1, ceiling, name)
    if limits.chunk_bytes < 32:
        raise CodecError("chunk_bytes must hold a complex slot")
    if (type(spec.kind) is not str or type(spec.storage_layout) is not str
            or spec.kind not in ("npy", "npz") or spec.storage_layout != STORAGE_LAYOUT):
        raise CodecError("unsupported kind or explicitly pinned storage layout")
    for label in ("file_sha256", "producer_sha256", "manifest_sha256", "layout_evidence_sha256"):
        _digest(getattr(spec, label), label)
    _integer(spec.file_bytes, 1, limits.max_file_bytes, "file_bytes")
    if type(spec.members) is not tuple or not 1 <= len(spec.members) <= limits.max_members:
        raise CodecError("invalid member roster")
    names: set[str] = set()
    total = 0
    for member in spec.members:
        if (type(member) is not MemberSpec or type(member.name) is not str
                or len(member.name) > 255 or not _NAME.fullmatch(member.name)):
            raise CodecError("invalid member name/spec")
        if member.name in names:
            raise CodecError("duplicate spec member")
        names.add(member.name)
        _digest(member.sha256, "member")
        _digest(member.header_sha256, "header")
        _integer(member.nbytes, 1, limits.max_member_bytes, "member bytes")
        _integer(member.header_bytes, 10, limits.max_header_bytes, "header bytes")
        if member.header_bytes % 16 or type(member.descr) is not str or member.descr not in _ITEMSIZE:
            raise CodecError("unsupported descriptor or header alignment")
        if type(member.shape) is not tuple or len(member.shape) > limits.max_dimensions:
            raise CodecError("invalid shape")
        count = 1
        for dim in member.shape:
            _integer(dim, 0, limits.max_elements, "dimension")
            count *= dim
            if count > limits.max_elements:
                raise CodecError("element cap exceeded")
        if (type(member.npy_version) is not tuple or len(member.npy_version) != 2
                or any(type(part) is not int for part in member.npy_version)
                or member.npy_version != (1, 0)):
            raise CodecError("only NPY 1.0 is supported")
        if member.fortran_order is not False or type(member.decode) is not bool:
            raise CodecError("only C-order arrays and explicit decode flags supported")
        if member.decode and member.descr not in ("<f16", "<c32"):
            raise CodecError("only binary80 real/complex members may be decoded")
        if member.nbytes != member.header_bytes + count * _ITEMSIZE[member.descr]:
            raise CodecError("spec byte count disagrees with shape")
        if spec.kind == "npz":
            if type(member.compression) is not int or member.compression not in (0, 8):
                raise CodecError("only ZIP_STORED/ZIP_DEFLATED supported")
            _integer(member.crc32, 0, (1 << 32) - 1, "CRC32")
        elif member.compression is not None or member.crc32 is not None:
            raise CodecError("standalone NPY must have no ZIP fields")
        total += member.nbytes
    if total > limits.max_total_member_bytes:
        raise CodecError("total member cap exceeded")
    if spec.kind == "npy" and (len(spec.members) != 1 or spec.members[0].nbytes != spec.file_bytes
                              or spec.members[0].sha256 != spec.file_sha256):
        raise CodecError("standalone NPY member must pin the complete file")
    return hashlib.sha256(json.dumps(asdict(spec), sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True).encode("ascii")).hexdigest()


def _spool(limits: Limits, temporary_directory: str | Path | None) -> BinaryIO:
    return tempfile.SpooledTemporaryFile(max_size=limits.spool_memory_bytes, mode="w+b",
                                       dir=temporary_directory)


def _snapshot(path: str | Path, spec: InputSpec, limits: Limits,
              temporary_directory: str | Path | None) -> BinaryIO:
    snapshot = _spool(limits, temporary_directory)
    try:
        flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        with os.fdopen(os.open(path, flags), "rb") as source:
            info = os.fstat(source.fileno())
            if not stat.S_ISREG(info.st_mode) or info.st_size != spec.file_bytes:
                raise CodecError("source must be regular and match pinned size")
            digest = hashlib.sha256()
            count = 0
            while True:
                block = source.read(min(limits.chunk_bytes, spec.file_bytes - count + 1))
                if not block:
                    break
                count += len(block)
                if count > spec.file_bytes:
                    raise CodecError("source grew beyond pin")
                snapshot.write(block)
                digest.update(block)
            if count != spec.file_bytes or digest.hexdigest() != spec.file_sha256:
                raise CodecError("complete-file SHA256/length mismatch")
        snapshot.seek(0)
        return snapshot
    except BaseException:
        snapshot.close()
        raise


def _read_exact(stream: BinaryIO, count: int) -> bytes:
    data = stream.read(count)
    if len(data) != count:
        raise CodecError("truncated input")
    return data


def _zip_preflight(snapshot: BinaryIO, spec: InputSpec, limits: Limits) -> None:
    # Bound the roster BEFORE ZipFile allocates/parses central-directory data.
    if spec.file_bytes < 22:
        raise CodecError("truncated ZIP trailer")
    snapshot.seek(spec.file_bytes - 22)
    trailer = struct.unpack("<4s4H2IH", _read_exact(snapshot, 22))
    magic, disk, cd_disk, n_disk, n_total, cd_size, cd_offset, comment = trailer
    if (magic != b"PK\x05\x06" or disk or cd_disk or comment or n_disk != n_total
            or n_total != len(spec.members) or cd_offset + cd_size != spec.file_bytes - 22
            or cd_size > limits.max_members * (46 + 255 + limits.max_header_bytes)):
        raise CodecError("unsupported/bounded ZIP central directory or trailer")
    position = cd_offset
    for member in spec.members:
        snapshot.seek(position)
        record = struct.unpack("<4s6H3I5H2I", _read_exact(snapshot, 46))
        if (record[0] != b"PK\x01\x02" or record[10] != len(member.name)
                or record[11] > limits.max_header_bytes or record[12] or record[13]):
            raise CodecError("invalid/bounded central-directory record")
        if _read_exact(snapshot, record[10]) != member.name.encode("ascii"):
            raise CodecError("central-directory roster mismatch")
        position += 46 + record[10] + record[11] + record[12]
        if position > cd_offset + cd_size:
            raise CodecError("central directory exceeded declared region")
    if position != cd_offset + cd_size:
        raise CodecError("extra central-directory records/bytes")
    snapshot.seek(0)


def _zip_layout(snapshot: BinaryIO, archive: zipfile.ZipFile, spec: InputSpec,
                limits: Limits) -> list[tuple[zipfile.ZipInfo, int]]:
    # Deliberately narrow ZIP container contract. Member-local ZIP64 size fields
    # are supported (the known capsule uses these); ZIP64 EOCD is unnecessary.
    snapshot.seek(spec.file_bytes - 22)
    eocd = _read_exact(snapshot, 22)
    magic, disk, cd_disk, n_disk, n_total, cd_size, cd_offset, comment = struct.unpack("<4s4H2IH", eocd)
    if (magic != b"PK\x05\x06" or disk or cd_disk or comment or n_disk != n_total
            or n_total != len(spec.members) or cd_offset + cd_size != spec.file_bytes - 22
            or archive.start_dir != cd_offset or archive.comment):
        raise CodecError("unsupported ZIP trailer, disks, comment, or trailing bytes")
    infos = archive.infolist()
    if len(infos) > limits.max_members or [i.filename for i in infos] != [m.name for m in spec.members]:
        raise CodecError("ZIP roster/order mismatch or duplicate member")
    position = 0
    regions = []
    for info, member in zip(infos, spec.members):
        if (info.filename != info.orig_filename or info.header_offset != position or info.comment
                or info.flag_bits & ~0x800 or stat.S_IFMT(info.external_attr >> 16) not in (0, stat.S_IFREG)
                or info.compress_type != member.compression
                or info.CRC != member.crc32 or info.file_size != member.nbytes
                or info.file_size > limits.max_member_bytes or info.compress_size < 1
                or info.compress_size > limits.max_file_bytes
                or info.file_size > info.compress_size * limits.max_decompression_ratio):
            raise CodecError("ZIP member metadata/pin/expansion limit mismatch")
        snapshot.seek(position)
        local = struct.unpack("<4s5H3I2H", _read_exact(snapshot, 30))
        sig, version, flags, method, _, _, crc, compressed, uncompressed, name_len, extra_len = local
        if (sig != b"PK\x03\x04" or version not in (20, 45) or flags != info.flag_bits
                or method != info.compress_type or crc != info.CRC
                or name_len != len(member.name) or extra_len > limits.max_header_bytes):
            raise CodecError("invalid ZIP local header")
        if _read_exact(snapshot, name_len) != member.name.encode("ascii"):
            raise CodecError("ZIP local filename mismatch")
        extra = _read_exact(snapshot, extra_len)
        if compressed == 0xffffffff or uncompressed == 0xffffffff:
            if (version != 45 or compressed != 0xffffffff or uncompressed != 0xffffffff
                    or len(extra) != 20 or extra[:4] != b"\x01\x00\x10\x00"):
                raise CodecError("unsupported local ZIP64 extra")
            uncompressed, compressed = struct.unpack("<QQ", extra[4:])
        elif extra:
            raise CodecError("unsupported ZIP local extra field")
        if compressed != info.compress_size or uncompressed != info.file_size:
            raise CodecError("ZIP local sizes mismatch")
        regions.append((info, position + 30 + name_len + extra_len))
        position += 30 + name_len + extra_len + compressed
    if position != cd_offset:
        raise CodecError("ZIP member regions do not exactly precede central directory")
    return regions


def _copy_zip_member(snapshot: BinaryIO, info: zipfile.ZipInfo, data_offset: int,
                     member: MemberSpec, limits: Limits, stream: BinaryIO) -> None:
    # ZipExtFile trusts/truncates to declared file_size. Read the complete raw
    # compressed region ourselves so forged sizes cannot hide an inflated tail.
    snapshot.seek(data_offset)
    remaining = info.compress_size
    digest = hashlib.sha256()
    crc = 0
    count = 0
    inflater = zlib.decompressobj(-15) if info.compress_type == 8 else None

    def record(block: bytes) -> None:
        nonlocal count, crc
        count += len(block)
        if count > member.nbytes:
            raise CodecError("actual decompression exceeded complete-member length/cap")
        digest.update(block)
        crc = zlib.crc32(block, crc)
        stream.write(block)

    while remaining:
        block = _read_exact(snapshot, min(limits.chunk_bytes, remaining))
        remaining -= len(block)
        if inflater is None:
            record(block)
            continue
        pending = block
        while pending:
            previous_length = len(pending)
            output = inflater.decompress(pending, min(limits.chunk_bytes, member.nbytes - count + 1))
            record(output)
            pending = inflater.unconsumed_tail
            if inflater.unused_data or (inflater.eof and (pending or remaining)):
                raise CodecError("deflate stream has trailing compressed bytes")
            if pending and len(pending) == previous_length and not inflater.eof:
                # Output may be emitted without consuming more compressed bytes;
                # when that happens count advances and its cap still guarantees
                # termination. A zero-output/no-consumption state is malformed.
                if not output:
                    raise CodecError("deflate decompressor made no progress")
        if inflater.eof and remaining:
            raise CodecError("deflate stream ended before its compressed region")
    if inflater is not None:
        while not inflater.eof:
            output = inflater.decompress(b"", min(limits.chunk_bytes, member.nbytes - count + 1))
            record(output)
            if not output:
                raise CodecError("deflate stream lacks end-of-stream marker")
        if inflater.unused_data or inflater.unconsumed_tail:
            raise CodecError("deflate stream did not end exactly at its region boundary")
    if count != member.nbytes or digest.hexdigest() != member.sha256 or crc != member.crc32:
        raise CodecError("complete-member SHA256/CRC32/actual length mismatch")


def _header(stream: BinaryIO, member: MemberSpec, limits: Limits) -> None:
    stream.seek(0)
    preamble = _read_exact(stream, 10)
    if preamble[:8] != b"\x93NUMPY\x01\x00":
        raise CodecError("only NPY 1.0 magic/version accepted")
    header_length = struct.unpack("<H", preamble[8:])[0]
    if header_length + 10 != member.header_bytes or header_length + 10 > limits.max_header_bytes:
        raise CodecError("header length disagrees with pin/cap")
    raw = _read_exact(stream, header_length)
    if not raw.endswith(b"\n") or b"\n" in raw[:-1] or b"\r" in raw or b"\t" in raw:
        raise CodecError("header must be a single newline-terminated line")
    if hashlib.sha256(preamble + raw).hexdigest() != member.header_sha256:
        raise CodecError("complete-header SHA256 mismatch")
    try:
        tree = ast.parse(raw[:-1].rstrip(b" ").decode("ascii"), mode="eval")
    except (ValueError, SyntaxError, UnicodeError, RecursionError) as error:
        raise CodecError("invalid bounded literal header") from error
    node = tree.body
    if type(node) is not ast.Dict or len(node.keys) != 3:
        raise CodecError("header must have exactly three literal dictionary keys")
    fields = {}
    for key, value in zip(node.keys, node.values):
        if type(key) is not ast.Constant or type(key.value) is not str or key.value in fields:
            raise CodecError("nonliteral/duplicate header key")
        fields[key.value] = value
    if set(fields) != {"descr", "fortran_order", "shape"}:
        raise CodecError("unsupported header keys")
    descr, order, shape = fields["descr"], fields["fortran_order"], fields["shape"]
    if type(descr) is not ast.Constant or type(descr.value) is not str or descr.value != member.descr:
        raise CodecError("descriptor mismatch; objects/structured/foreign layouts unsupported")
    if type(order) is not ast.Constant or order.value is not False:
        raise CodecError("fortran_order must be literal False")
    if type(shape) is not ast.Tuple or len(shape.elts) > limits.max_dimensions:
        raise CodecError("shape must be a bounded literal tuple")
    dimensions = []
    for dim in shape.elts:
        if type(dim) is not ast.Constant or type(dim.value) is not int:
            raise CodecError("shape dimensions must be literal integers")
        dimensions.append(_integer(dim.value, 0, limits.max_elements, "header dimension"))
    if tuple(dimensions) != member.shape:
        raise CodecError("shape mismatch")


def _decode_slot(slot: bytes) -> ExactReal:
    if len(slot) != 16:
        raise CodecError("binary80 slot requires exactly 16 bytes")
    significand = int.from_bytes(slot[:8], "little")
    sign_exponent = int.from_bytes(slot[8:10], "little")
    negative = bool(sign_exponent & 0x8000)
    exponent = sign_exponent & 0x7fff
    integer_bit = significand >> 63
    if exponent == 0x7fff:
        raise CodecError("NaN/infinity/pseudo-special binary80 unsupported")
    if exponent == 0:
        if integer_bit:
            raise CodecError("pseudo-denormal binary80 unsupported")
        power = -16445
    else:
        if not integer_bit:
            raise CodecError("unnormal binary80 unsupported")
        power = exponent - 16446
    if significand == 0:
        return ExactReal(Fraction(0), negative)
    numerator = -significand if negative else significand
    value = Fraction(numerator << power, 1) if power >= 0 else Fraction(numerator, 1 << -power)
    return ExactReal(value, False)


class VerifiedInputs:
    """Owns verified immutable byte snapshots; close after consuming iterators."""
    def __init__(self, spec: InputSpec, limits: Limits, spec_digest: str,
                 streams: dict[str, BinaryIO]):
        self.spec = spec
        self.limits = limits
        self.input_spec_sha256 = spec_digest
        self._streams = streams
        self._members = {m.name: m for m in spec.members}
        self._locks = {m.name: threading.RLock() for m in spec.members}
        self._closed = False

    def __enter__(self) -> VerifiedInputs:
        self._check_open()
        return self

    def __exit__(self, *_: object) -> None:
        self.close()

    def _check_open(self) -> None:
        if self._closed:
            raise CodecError("verified input is closed")

    def close(self) -> None:
        for stream in self._streams.values():
            stream.close()
        self._closed = True

    def checkpoint(self, member_name: str, next_index: int) -> CodecCheckpoint:
        self._check_open()
        member = self._members.get(member_name)
        if member is None or not member.decode:
            raise CodecError("member is not approved for decoding")
        count = (member.nbytes - member.header_bytes) // _ITEMSIZE[member.descr]
        _integer(next_index, 0, count, "checkpoint index")
        return CodecCheckpoint(self.input_spec_sha256, self.spec.file_sha256, member_name,
                               member.sha256, next_index)

    def iter_values(self, member_name: str, *, start: int = 0, stop: int | None = None,
                    checkpoint: CodecCheckpoint | None = None) -> Iterator[ExactReal | ExactComplex]:
        self._check_open()
        member = self._members.get(member_name)
        if member is None or not member.decode:
            raise CodecError("member is not approved for decoding")
        if checkpoint is not None:
            if type(checkpoint) is not CodecCheckpoint or start != 0:
                raise CodecError("explicit checkpoint with default start required")
            expected = self.checkpoint(member_name, checkpoint.next_index)
            if checkpoint != expected:
                raise CodecError("checkpoint belongs to different bytes/spec/member")
            start = checkpoint.next_index
        width = _ITEMSIZE[member.descr]
        count = (member.nbytes - member.header_bytes) // width
        _integer(start, 0, count, "start")
        end = count if stop is None else _integer(stop, start, count, "stop")
        per_chunk = max(1, self.limits.chunk_bytes // width)
        stream = self._streams[member_name]
        for base in range(start, end, per_chunk):
            self._check_open()
            with self._locks[member_name]:
                stream.seek(member.header_bytes + base * width)
                block = _read_exact(stream, min(per_chunk, end - base) * width)
            for offset in range(0, len(block), width):
                real = _decode_slot(block[offset:offset + 16])
                yield real if width == 16 else ExactComplex(real, _decode_slot(block[offset + 16:offset + 32]))


def verify_inputs(path: str | Path, input_spec: InputSpec, limits: Limits, *,
                  temporary_directory: str | Path | None = None) -> VerifiedInputs:
    """Verify an opaque snapshot completely, then permit explicitly selected decode.

    The call never decodes numeric payloads. All hashes include binary80 padding.
    Errors close every owned temporary stream. No original pathname is reopened.
    """
    spec_digest = _validate_spec(input_spec, limits)
    streams: dict[str, BinaryIO] = {}
    snapshot = _snapshot(path, input_spec, limits, temporary_directory)
    try:
        if input_spec.kind == "npy":
            streams[input_spec.members[0].name] = snapshot
        else:
            _zip_preflight(snapshot, input_spec, limits)
            with zipfile.ZipFile(snapshot, "r") as archive:
                regions = _zip_layout(snapshot, archive, input_spec, limits)
                for (info, data_offset), member in zip(regions, input_spec.members):
                    stream = _spool(limits, temporary_directory)
                    streams[member.name] = stream
                    _copy_zip_member(snapshot, info, data_offset, member, limits, stream)
            snapshot.close()
        # Every complete member hash has passed before parsing any NPY header.
        for member in input_spec.members:
            _header(streams[member.name], member, limits)
        return VerifiedInputs(input_spec, limits, spec_digest, streams)
    except BaseException as error:
        snapshot.close()
        for stream in streams.values():
            stream.close()
        if isinstance(error, (zipfile.BadZipFile, NotImplementedError, EOFError, RuntimeError, struct.error, zlib.error)):
            raise CodecError("invalid/unsupported ZIP container") from error
        raise
