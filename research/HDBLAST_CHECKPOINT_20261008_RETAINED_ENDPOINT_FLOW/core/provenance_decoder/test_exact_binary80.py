"""Meaningful canonical/malformed/pinning tests using fabricated bytes only."""
from dataclasses import replace
from fractions import Fraction
import io
import os
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch
import warnings
import zipfile
from exact_binary80 import CodecError, CodecCheckpoint, ExactComplex, ExactReal, verify_inputs
from fabricated_fixtures import LIMITS, archive, input_spec, member, npy, repin_file, sha, slot, standalone


class FabricatedCodecTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="binary80-fabricated-")
        self.addCleanup(self.tmp.cleanup)
        self.path = Path(self.tmp.name) / "fabricated-input.bin"

    def read(self, data, spec=None, **kwargs):
        self.path.write_bytes(data)
        return verify_inputs(self.path, spec or standalone(data), kwargs.pop("limits", LIMITS),
                             temporary_directory=self.tmp.name, **kwargs)

    def reject(self, data, spec=None, limits=LIMITS):
        with patch("exact_binary80._decode_slot", side_effect=RuntimeError("must not decode")) as decode:
            with self.assertRaises(CodecError):
                self.read(data, spec, limits=limits)
            self.assertEqual(decode.call_count, 0)

    def test_normal_sign_endian_and_ulp(self):
        cases = [(slot(), Fraction(1)), (slot(3 << 62), Fraction(3, 2)),
                 (slot((1 << 63) + 1), 1 + Fraction(1, 1 << 63)),
                 (slot(1 << 63, 16381, True), Fraction(-1, 4))]
        for raw, expected in cases:
            with self.subTest(expected=str(expected)):
                with self.read(npy(raw)) as reader:
                    self.assertEqual(list(reader.iter_values("array.npy")), [ExactReal(expected, False)])

    def test_signed_zero_and_arbitrary_padding(self):
        values = [slot(0, 0, False, b"ABCDEF"), slot(0, 0, True, b"\xff" * 6),
                  slot(padding=b"\xde\xad\xbe\xef\x10\x20")]
        data = npy(b"".join(values), shape=(3,))
        with self.read(data, standalone(data, shape=(3,))) as reader:
            self.assertEqual(list(reader.iter_values("array.npy")),
                             [ExactReal(Fraction(0), False), ExactReal(Fraction(0), True), ExactReal(Fraction(1), False)])
        changed = data[:-1] + bytes([data[-1] ^ 1])
        self.reject(changed, standalone(data, shape=(3,)))

    def test_subnormal_and_extreme_boundaries(self):
        raw = [slot(1, 0), slot((1 << 63) - 1, 0), slot(1 << 63, 1),
               slot((1 << 64) - 1, 0x7ffe), slot(1, 0, True)]
        expected = [Fraction(1, 1 << 16445), Fraction((1 << 63) - 1, 1 << 16445),
                    Fraction(1, 1 << 16382), Fraction(((1 << 64) - 1) << 16320),
                    -Fraction(1, 1 << 16445)]
        data = npy(b"".join(raw), shape=(5,))
        with self.read(data, standalone(data, shape=(5,))) as reader:
            self.assertEqual([v.value for v in reader.iter_values("array.npy")], expected)

    def test_complex_component_order_and_signed_zero(self):
        data = npy(slot(3 << 62) + slot(0, 0, True, b"abcdef"), "<c32")
        with self.read(data, standalone(data, "<c32")) as reader:
            self.assertEqual(list(reader.iter_values("array.npy")),
                             [ExactComplex(ExactReal(Fraction(3, 2), False), ExactReal(Fraction(0), True))])

    def test_noncanonical_and_special_are_rejected(self):
        for raw in [slot(1 << 63, 0), slot(1, 1), slot(0, 16383),
                    slot(1 << 63, 0x7fff), slot((1 << 63) + 1, 0x7fff), slot(0, 0x7fff)]:
            with self.subTest(raw=raw.hex()):
                with self.read(npy(raw)) as reader:
                    with self.assertRaises(CodecError):
                        list(reader.iter_values("array.npy"))

    def test_no_value_decode_during_complete_verification(self):
        data, spec = archive([("first.npy", npy(slot()), "<f16", (1,), True),
                              ("last.npy", npy(slot()), "<f16", (1,), True)])
        with patch("exact_binary80._decode_slot", side_effect=RuntimeError("must not decode")) as decode:
            with self.read(data, spec):
                self.assertEqual(decode.call_count, 0)
        bad = replace(spec, members=spec.members[:-1] + (replace(spec.members[-1], sha256=sha(b"wrong last")),))
        self.reject(data, bad)

    def test_zip_stored_deflated_and_local_zip64(self):
        for compression in (0, 8):
            for force in (False, True):
                data, spec = archive([("x.npy", npy(slot(padding=b"123456")), "<f16", (1,), True)],
                                     compression=compression, force_zip64=force)
                with self.read(data, spec) as reader:
                    self.assertEqual(list(reader.iter_values("x.npy")), [ExactReal(Fraction(1), False)])

    def test_opaque_pinned_members_verified_but_not_decoded(self):
        data, spec = archive([("x.npy", npy(slot()), "<f16", (1,), True),
                              ("i.npy", npy(struct.pack("<q", 19), "<i8"), "<i8", (1,), False),
                              ("text.npy", npy(b"\x00" * 32, "<U8"), "<U8", (1,), False)])
        with self.read(data, spec) as reader:
            with self.assertRaises(CodecError):
                list(reader.iter_values("i.npy"))
            self.assertEqual(len(list(reader.iter_values("x.npy"))), 1)

    def test_strict_headers(self):
        literals = ["{'descr':'<f16','fortran_order':False,'shape':(1,), 'descr':'<f16'}",
                    "{'descr':'<f16','fortran_order':0,'shape':(1,)}",
                    "{'descr':'<f16','fortran_order':False,'shape':(True,)}",
                    "{'descr':'<f16','fortran_order':False,'shape':[1]}",
                    "{'descr':'<f16','fortran_order':False,'shape':(65537,)}",
                    "{'descr':'<f16','fortran_order':False,'shape':(1,), 'extra':None}",
                    "dict(descr='<f16',fortran_order=False,shape=(1,))",
                    "{'descr':'|O','fortran_order':False,'shape':(1,)}",
                    "{'descr':'>f16','fortran_order':False,'shape':(1,)}",
                    "{'descr':'<f16','fortran_order':False,'shape':(1+0,)}"]
        for literal in literals:
            with self.subTest(literal=literal):
                self.reject(npy(slot(), literal=literal))
        self.reject(npy(slot(), newline=b" "))
        self.reject(npy(slot(), version=b"\x02\x00"))
        data = npy(slot())
        self.reject(data[:8] + b"\xff\xff" + data[10:], standalone(data))

    def test_payload_lengths_and_trailing_bytes(self):
        data = npy(slot())
        self.reject(data + b"extra", repin_file(standalone(data), data + b"extra"))
        self.reject(data[:-1], repin_file(standalone(data), data[:-1]))
        self.reject(data, replace(standalone(data), members=(replace(standalone(data).members[0], shape=(2,)),)))

    def test_zip_roster_duplicate_unlisted_trailing(self):
        entries = [("x.npy", npy(slot()), "<f16", (1,), True)]
        data, spec = archive(entries)
        self.reject(data + b"trailing", repin_file(spec, data + b"trailing"))
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            duplicate, _ = archive(entries + entries)
        self.reject(duplicate, repin_file(spec, duplicate))
        other, _ = archive(entries + [("other.npy", npy(slot()), "<f16", (1,), True)])
        self.reject(other, repin_file(spec, other))
        self.reject(data, replace(spec, members=(replace(spec.members[0], name="../x.npy"),)))

    def test_zip_crc_and_local_metadata_mismatch(self):
        data, spec = archive([("x.npy", npy(slot()), "<f16", (1,), True)])
        self.reject(data, replace(spec, members=(replace(spec.members[0], crc32=0),)))
        bad = bytearray(data)
        bad[14] ^= 1  # local CRC; outer pin deliberately updated for malformed test
        self.reject(bytes(bad), repin_file(spec, bytes(bad)))
        bad = bytearray(data)
        payload_position = 30 + len("x.npy") + 128
        bad[payload_position] ^= 1
        self.reject(bytes(bad), repin_file(spec, bytes(bad)))

    def test_hash_and_resource_caps(self):
        data = npy(slot())
        spec = standalone(data)
        self.reject(data, replace(spec, file_sha256=sha(b"wrong")))
        self.reject(data, replace(spec, members=(replace(spec.members[0], header_sha256=sha(b"wrong")),)))
        self.reject(data, spec, replace(LIMITS, max_file_bytes=len(data) - 1))
        self.reject(data, spec, replace(LIMITS, max_member_bytes=len(data) - 1))
        self.reject(data, spec, replace(LIMITS, max_header_bytes=64))
        self.reject(data, spec, replace(LIMITS, max_elements=True))
        for version in ((True, False), (1.0, 0.0)):
            self.reject(data, replace(spec, members=(replace(spec.members[0], npy_version=version),)))
        compressed, cspec = archive([("x.npy", npy(slot(0, 0) * 4096, shape=(4096,)),
                                     "<f16", (4096,), True)], compression=8)
        self.reject(compressed, cspec, replace(LIMITS, max_decompression_ratio=2))

    def test_empty_scalar_multidimensional_and_checkpoints(self):
        for shape, payload, count in [((), slot(), 1), ((0,), b"", 0), ((2, 2), slot() * 4, 4)]:
            data = npy(payload, shape=shape)
            with self.read(data, standalone(data, shape=shape)) as reader:
                self.assertEqual(len(list(reader.iter_values("array.npy"))), count)
        data = npy(b"".join(slot((1 << 63) + i) for i in range(5)), shape=(5,))
        with self.read(data, standalone(data, shape=(5,))) as reader:
            all_values = list(reader.iter_values("array.npy"))
            saved = reader.checkpoint("array.npy", 2)
            self.assertEqual(list(reader.iter_values("array.npy", checkpoint=saved)), all_values[2:])
            self.assertEqual(list(reader.iter_values("array.npy", start=1, stop=3)), all_values[1:3])
            with self.assertRaises(CodecError):
                list(reader.iter_values("array.npy", checkpoint=replace(saved, member_sha256=sha(b"foreign"))))
        with self.assertRaises(CodecError):
            list(reader.iter_values("array.npy"))

    def test_snapshot_prevents_path_replacement_and_symlink(self):
        data = npy(slot())
        with self.read(data) as reader:
            self.path.write_bytes(npy(slot(3 << 62)))
            self.assertEqual(list(reader.iter_values("array.npy")), [ExactReal(Fraction(1), False)])
        linked = Path(self.tmp.name) / "symlink.bin"
        linked.symlink_to(self.path)
        with self.assertRaises(OSError):
            verify_inputs(linked, standalone(npy(slot(3 << 62))), LIMITS)


if __name__ == "__main__":
    unittest.main()
