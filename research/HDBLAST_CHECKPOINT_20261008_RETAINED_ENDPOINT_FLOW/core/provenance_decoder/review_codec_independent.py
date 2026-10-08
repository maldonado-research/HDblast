"""Independent stdlib-only fabricated-fixture checks; never opens real payloads."""
from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from io import BytesIO
from pathlib import Path
import hashlib
import struct
import tempfile
import unittest
from unittest.mock import patch
import zipfile
import zlib

import exact_binary80 as c

ROOT = Path(__file__).resolve().parent
LIMITS = c.Limits(1 << 20, 1 << 18, 1 << 19, 16, 4096, 8, 4096, 1024, 32, 32)
EVIDENCE = hashlib.sha256(b"fabricated fixture layout evidence").hexdigest()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def slot(significand, exponent, negative=False, padding=b"\xa5" * 6):
    return struct.pack("<QH", significand, exponent | (0x8000 if negative else 0)) + padding


def npy(descr, shape, payload, text=None):
    if text is None:
        text = repr({"descr": descr, "fortran_order": False, "shape": shape})
    raw = text.encode("ascii")
    raw += b" " * ((-(10 + len(raw) + 1)) % 16) + b"\n"
    return b"\x93NUMPY\x01\x00" + struct.pack("<H", len(raw)) + raw + payload


def member(name, data, descr, shape, decode=True, compression=None, crc=None):
    header_bytes = 10 + struct.unpack("<H", data[8:10])[0]
    return c.MemberSpec(name, digest(data), len(data), digest(data[:header_bytes]),
                        header_bytes, descr, shape, (1, 0), False, decode,
                        compression, crc)


def spec(data, members, kind="npy"):
    return c.InputSpec(kind, digest(data), len(data), c.STORAGE_LAYOUT,
                       EVIDENCE, EVIDENCE, EVIDENCE, tuple(members))


def archive(entries, compression=zipfile.ZIP_STORED, force_zip64=False):
    out = BytesIO()
    with zipfile.ZipFile(out, "w", compression=compression) as z:
        for name, data in entries:
            if force_zip64:
                with z.open(name, "w", force_zip64=True) as stream:
                    stream.write(data)
            else:
                z.writestr(name, data)
    blob = out.getvalue()
    with zipfile.ZipFile(BytesIO(blob)) as z:
        infos = z.infolist()
    return blob, infos


class IndependentCodecReview(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory(prefix=".codec-review-", dir=ROOT)
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / "fabricated.bin"

    def verify(self, data, specification, limits=LIMITS):
        self.path.write_bytes(data)
        return c.verify_inputs(self.path, specification, limits,
                               temporary_directory=self.directory.name)

    def assert_verification_rejects_without_decode(self, data, specification, limits=LIMITS):
        with patch.object(c, "_decode_slot", side_effect=AssertionError("decoded before verification")) as decoder:
            with self.assertRaises(c.CodecError):
                with self.verify(data, specification, limits):
                    pass
            self.assertEqual(decoder.call_count, 0)

    def test_canonical_binary80_boundaries_and_signed_zero(self):
        cases = [
            (0, 0, False, Fraction(0), False),
            (0, 0, True, Fraction(0), True),
            (1, 0, False, Fraction(1, 1 << 16445), False),
            ((1 << 63) - 1, 0, False, Fraction((1 << 63) - 1, 1 << 16445), False),
            (1 << 63, 1, False, Fraction(1, 1 << 16382), False),
            (1 << 63, 16383, False, Fraction(1), False),
            ((1 << 63) + 1, 16383, False, Fraction((1 << 63) + 1, 1 << 63), False),
            (1 << 63, 16383, True, Fraction(-1), False),
            ((1 << 64) - 1, 32766, False, Fraction(((1 << 64) - 1) << 16320), False),
        ]
        for significand, exponent, negative, value, negative_zero in cases:
            with self.subTest(significand=significand, exponent=exponent, negative=negative):
                result = c._decode_slot(slot(significand, exponent, negative))
                self.assertEqual(result, c.ExactReal(value, negative_zero))

    def test_padding_is_ignored_numerically(self):
        a = c._decode_slot(slot((1 << 63) + 11, 16387, padding=b"\x00" * 6))
        b = c._decode_slot(slot((1 << 63) + 11, 16387, padding=b"\xff" * 6))
        self.assertEqual(a, b)

    def test_noncanonical_and_nonfinite_slots_rejected(self):
        for significand, exponent in [(1 << 63, 0), (0, 1), (1, 17),
                                      (1 << 63, 32767), ((1 << 64) - 1, 32767),
                                      (0, 32767), (1, 32767)]:
            for negative in (False, True):
                with self.subTest(significand=significand, exponent=exponent, negative=negative):
                    with self.assertRaises(c.CodecError):
                        c._decode_slot(slot(significand, exponent, negative))
        for length in (0, 10, 15, 17, 32):
            with self.subTest(length=length), self.assertRaises(c.CodecError):
                c._decode_slot(bytes(length))

    def test_end_to_end_complex_and_checkpoint(self):
        payload = (slot(0, 0, True) + slot(1 << 63, 16383)
                   + slot((1 << 63) + 1, 16383) + slot(0, 0, True))
        data = npy("<c32", (2,), payload)
        m = member("array.npy", data, "<c32", (2,))
        with self.verify(data, spec(data, [m])) as verified:
            values = list(verified.iter_values(m.name))
            self.assertEqual(values, [c.ExactComplex(c.ExactReal(Fraction(0), True),
                                                    c.ExactReal(Fraction(1), False)),
                                      c.ExactComplex(c.ExactReal(Fraction((1 << 63) + 1, 1 << 63), False),
                                                     c.ExactReal(Fraction(0), True))])
            checkpoint = verified.checkpoint(m.name, 1)
            self.assertEqual(list(verified.iter_values(m.name, checkpoint=checkpoint)), values[1:])
            with self.assertRaises(c.CodecError):
                list(verified.iter_values(m.name, checkpoint=replace(checkpoint, member_sha256="0" * 64)))
            with self.assertRaises(c.CodecError):
                list(verified.iter_values(m.name, start=True))

    def test_no_scalar_decode_during_successful_verification(self):
        data = npy("<f16", (1,), slot(1 << 63, 16383))
        m = member("array.npy", data, "<f16", (1,))
        with patch.object(c, "_decode_slot", side_effect=AssertionError("early decode")) as decoder:
            with self.verify(data, spec(data, [m])):
                pass
            self.assertEqual(decoder.call_count, 0)

    def test_full_file_mismatch_and_padding_change_rejected_before_decode(self):
        data = npy("<f16", (1,), slot(1 << 63, 16383))
        m = member("array.npy", data, "<f16", (1,))
        changed = data[:-1] + bytes([data[-1] ^ 1])
        self.assert_verification_rejects_without_decode(changed, spec(data, [m]))
        self.assert_verification_rejects_without_decode(data + b"x", spec(data, [m]))

    def test_malformed_literal_headers_rejected_before_decode(self):
        texts = [
            "{'descr': '<f16', 'fortran_order': False, 'shape': (True,)}",
            "{'descr': '<f16', 'fortran_order': 0, 'shape': (1,)}",
            "{'descr': '<f16', 'fortran_order': False, 'shape': [1]}",
            "{'descr': '<f16', 'fortran_order': False, 'shape': (1+0,)}",
            "{'descr': '<f16', 'fortran_order': False, 'shape': (1,), 'shape': (1,)}",
            "{'descr': '<f16', 'fortran_order': False, 'shape': (1,), 'unknown': 0}",
            "{'descr': '<f16', 'fortran_order': False, 'shape': (1,)} or {}",
            "dict(descr='<f16', fortran_order=False, shape=(1,))",
            "{'descr': '<f16', 'fortran_order': False, 'shape': (-1,)}",
            "{'descr': '<f16', 'fortran_order': False, 'shape': (4097,)}",
            "{'descr': '<f16', 'fortran_order': False, 'shape': (1.0,)}",
        ]
        for text in texts:
            with self.subTest(text=text):
                data = npy("<f16", (1,), slot(1 << 63, 16383), text)
                m = member("array.npy", data, "<f16", (1,))
                self.assert_verification_rejects_without_decode(data, spec(data, [m]))

    def test_header_pin_and_version_mismatch_rejected_before_decode(self):
        good = npy("<f16", (1,), slot(1 << 63, 16383))
        m = member("array.npy", good, "<f16", (1,))
        self.assert_verification_rejects_without_decode(good, spec(good, [replace(m, header_sha256="0" * 64)]))
        for position, replacement in [(6, b"\x02"), (7, b"\x01"), (0, b"x")]:
            bad = good[:position] + replacement + good[position + 1:]
            bad_m = member("array.npy", bad, "<f16", (1,))
            self.assert_verification_rejects_without_decode(bad, spec(bad, [bad_m]))

    def test_empty_and_scalar_shapes(self):
        for shape, payload, expected in [((), slot(1 << 63, 16383), [c.ExactReal(Fraction(1), False)]),
                                         ((0,), b"", [])]:
            with self.subTest(shape=shape):
                data = npy("<f16", shape, payload)
                m = member("array.npy", data, "<f16", shape)
                with self.verify(data, spec(data, [m])) as verified:
                    self.assertEqual(list(verified.iter_values(m.name)), expected)

    def test_stored_deflated_and_local_zip64(self):
        data = npy("<f16", (2,), slot(1 << 63, 16383) + slot(0, 0, True))
        for compression in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED):
            for force_zip64 in (False, True):
                with self.subTest(compression=compression, force_zip64=force_zip64):
                    blob, infos = archive([("array.npy", data)], compression, force_zip64)
                    m = member("array.npy", data, "<f16", (2,), compression=infos[0].compress_type,
                               crc=infos[0].CRC)
                    with self.verify(blob, spec(blob, [m], "npz")) as verified:
                        self.assertEqual(list(verified.iter_values(m.name)),
                                         [c.ExactReal(Fraction(1), False), c.ExactReal(Fraction(0), True)])

    def test_later_member_hash_failure_never_decodes_first_member(self):
        first = npy("<f16", (1,), slot(1 << 63, 16383))
        second = npy("<f16", (1,), slot(0, 0))
        altered = second[:-1] + bytes([second[-1] ^ 1])
        blob, infos = archive([("first.npy", first), ("second.npy", altered)])
        members = [member("first.npy", first, "<f16", (1,), compression=0, crc=infos[0].CRC),
                   member("second.npy", altered, "<f16", (1,), compression=0, crc=infos[1].CRC)]
        members[1] = replace(members[1], sha256=digest(second))
        self.assert_verification_rejects_without_decode(blob, spec(blob, members, "npz"))

    def test_actual_output_tail_rejected_before_decode(self):
        data = npy("<f16", (1,), slot(1 << 63, 16383))
        for compression in (zipfile.ZIP_STORED, zipfile.ZIP_DEFLATED):
            with self.subTest(compression=compression):
                blob, _ = archive([("array.npy", data + b"UNPINNED-DECOMPRESSED-TAIL")], compression)
                raw = bytearray(blob)
                central = struct.unpack_from("<I", raw, len(raw) - 6)[0]
                crc = zlib.crc32(data)
                struct.pack_into("<I", raw, 14, crc)
                struct.pack_into("<I", raw, 22, len(data))
                struct.pack_into("<I", raw, central + 16, crc)
                struct.pack_into("<I", raw, central + 24, len(data))
                blob = bytes(raw)
                m = member("array.npy", data, "<f16", (1,), compression=compression, crc=crc)
                self.assert_verification_rejects_without_decode(blob, spec(blob, [m], "npz"))

    def test_deflate_compressed_junk_and_incomplete_stream_rejected(self):
        data = npy("<f16", (1,), slot(1 << 63, 16383))
        original, _ = archive([("array.npy", data)], zipfile.ZIP_DEFLATED)
        oldcentral = struct.unpack_from("<I", original, len(original) - 6)[0]
        oldsize = struct.unpack_from("<I", original, 18)[0]
        start = 30 + len("array.npy")
        inflater = zlib.decompressobj(-15)
        self.assertEqual(inflater.decompress(original[start:start + oldsize - 1]), data)
        self.assertFalse(inflater.eof)
        junk = b"IGNORED-COMPRESSED-JUNK"
        variants = [("unused compressed bytes", original[:oldcentral] + junk + original[oldcentral:], len(junk)),
                    ("incomplete raw deflate", original[:oldcentral - 1] + original[oldcentral:], -1)]
        for label, initial, delta in variants:
            with self.subTest(label=label):
                raw = bytearray(initial)
                central = oldcentral + delta
                struct.pack_into("<I", raw, 18, oldsize + delta)
                struct.pack_into("<I", raw, central + 20, oldsize + delta)
                struct.pack_into("<I", raw, len(raw) - 6, central)
                blob = bytes(raw)
                m = member("array.npy", data, "<f16", (1,), compression=8, crc=zlib.crc32(data))
                self.assert_verification_rejects_without_decode(blob, spec(blob, [m], "npz"))

    def test_extra_central_records_rejected_before_zipfile_allocation(self):
        data = npy("<f16", (1,), slot(1 << 63, 16383))
        original, infos = archive([("array.npy", data)])
        central = struct.unpack_from("<I", original, len(original) - 6)[0]
        record = original[central:-22]
        raw = bytearray(original[:central] + record * 64 + original[-22:])
        struct.pack_into("<I", raw, len(raw) - 10, len(record) * 64)
        blob = bytes(raw)
        m = member("array.npy", data, "<f16", (1,), compression=0, crc=infos[0].CRC)
        with patch.object(c.zipfile, "ZipFile", side_effect=AssertionError("unbounded central-directory parse")) as parser:
            self.assert_verification_rejects_without_decode(blob, spec(blob, [m], "npz"))
            self.assertEqual(parser.call_count, 0)

    def test_later_header_failure_never_decodes_first_member(self):
        first = npy("<f16", (1,), slot(1 << 63, 16383))
        second = npy("<f16", (1,), slot(0, 0), "{'descr':'<f16','fortran_order':True,'shape':(1,)}")
        blob, infos = archive([("first.npy", first), ("second.npy", second)])
        members = [member("first.npy", first, "<f16", (1,), compression=0, crc=infos[0].CRC),
                   member("second.npy", second, "<f16", (1,), compression=0, crc=infos[1].CRC)]
        self.assert_verification_rejects_without_decode(blob, spec(blob, members, "npz"))

    def test_verify_only_metadata_members_not_decodable(self):
        data = npy("<i8", (1,), struct.pack("<q", 13))
        m = member("metadata.npy", data, "<i8", (1,), decode=False)
        with self.verify(data, spec(data, [m])) as verified:
            with self.assertRaises(c.CodecError):
                list(verified.iter_values(m.name))

    def test_zip_roster_crc_trailing_and_cap_rejections(self):
        data = npy("<f16", (1,), slot(1 << 63, 16383))
        blob, infos = archive([("array.npy", data)])
        m = member("array.npy", data, "<f16", (1,), compression=0, crc=infos[0].CRC)
        self.assert_verification_rejects_without_decode(blob, spec(blob, [replace(m, crc32=m.crc32 ^ 1)], "npz"))
        trailing = blob + b"X"
        self.assert_verification_rejects_without_decode(trailing, spec(trailing, [m], "npz"))
        duplicate, _ = archive([("array.npy", data), ("array.npy", data)])
        self.assert_verification_rejects_without_decode(duplicate, spec(duplicate, [m], "npz"))
        extra, _ = archive([("array.npy", data), ("extra.npy", data)])
        self.assert_verification_rejects_without_decode(extra, spec(extra, [m], "npz"))
        self.assert_verification_rejects_without_decode(blob, spec(blob, [m], "npz"),
                                                         replace(LIMITS, max_member_bytes=len(data) - 1))

    def test_unpinned_source_changes_do_not_change_snapshot(self):
        data = npy("<f16", (1,), slot(1 << 63, 16383))
        m = member("array.npy", data, "<f16", (1,))
        with self.verify(data, spec(data, [m])) as verified:
            self.path.write_bytes(npy("<f16", (1,), slot(1 << 63, 16384)))
            self.assertEqual(list(verified.iter_values(m.name)), [c.ExactReal(Fraction(1), False)])

    def test_npy_version_types_strict(self):
        data = npy("<f16", (1,), slot(1 << 63, 16383))
        m = member("array.npy", data, "<f16", (1,))
        for version in ((True, False), (1.0, 0.0), (1, False), (True, 0)):
            with self.subTest(version=version):
                self.assert_verification_rejects_without_decode(data, spec(data, [replace(m, npy_version=version)]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
