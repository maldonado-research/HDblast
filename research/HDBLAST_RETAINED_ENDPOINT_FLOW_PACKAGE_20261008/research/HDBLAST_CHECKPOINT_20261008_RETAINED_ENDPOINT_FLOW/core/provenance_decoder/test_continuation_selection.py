"""Integration of the nine selected streams using manufactured bytes only."""
from dataclasses import asdict
from fractions import Fraction
from pathlib import Path
import struct
import tempfile
import unittest
from unittest.mock import patch

import exact_binary80 as codec
from fabricated_fixtures import LIMITS, archive, npy, slot
from input_spec_adapter import SELECTED_MEMBERS, hydrate_reader_spec, load_metadata, metadata_receipt


class ContinuationSelectionTests(unittest.TestCase):
    def test_manufactured_full_roster_and_exact_nine_selected_streams(self):
        roster = [
            "k.npy", "momentum_weights.npy", "u_1.npy", "w_1.npy", "u_2.npy", "w_2.npy",
            "u_3.npy", "w_3.npy", "observation_eta.npy", "history_eta.npy", "history_values.npy",
            "history_ledger_integrand.npy", "history_baseline_contact.npy", "history_source_jet.npy",
            "history_forcing_jet.npy", "history_geometry.npy", "history_cutoffs.npy",
            "history_quantity_names.npy", "history_geometry_names.npy",
        ]
        entries = []
        expected = {}
        for index, name in enumerate(roster):
            descr = "<f16"
            payload = slot((1 << 63) + index, padding=b"ABCDEF")
            if name.startswith(("u_", "w_")):
                descr = "<c32"
                payload += slot(0, 0, True, padding=b"UVWXYZ")
            elif name == "history_cutoffs.npy":
                descr, payload = "<i8", struct.pack("<q", 64)
            elif name in ("history_quantity_names.npy", "history_geometry_names.npy"):
                descr, payload = "<U8", bytes(32)
            # Noncanonical binary80 in an opaque history is never interpreted.
            elif name.startswith("history_"):
                payload = slot(1, 1)
            entries.append((name, npy(payload, descr), descr, (1,), name in SELECTED_MEMBERS))
            expected[name] = Fraction((1 << 63) + index, 1 << 63)
        data, original = archive(entries, force_zip64=True)
        spec = hydrate_reader_spec(asdict(original))
        with tempfile.TemporaryDirectory(prefix="continuation-manufactured-") as directory:
            path = Path(directory) / "manufactured-only.npz"
            path.write_bytes(data)
            with patch.object(codec, "_decode_slot", wraps=codec._decode_slot) as decoder:
                with codec.verify_inputs(path, spec, LIMITS, temporary_directory=directory) as verified:
                    self.assertEqual(decoder.call_count, 0)
                    for name in SELECTED_MEMBERS:
                        values = list(verified.iter_values(name))
                        self.assertEqual(len(values), 1)
                        value = values[0]
                        real = value.real if isinstance(value, codec.ExactComplex) else value
                        self.assertEqual(real.value, expected[name])
                        if isinstance(value, codec.ExactComplex):
                            self.assertEqual(value.imag.value, 0)
                            self.assertTrue(value.imag.negative_zero)
                    self.assertEqual(decoder.call_count, 15)
                    for name in roster:
                        if name not in SELECTED_MEMBERS:
                            with self.assertRaises(codec.CodecError):
                                list(verified.iter_values(name))
                    self.assertEqual(decoder.call_count, 15)

    def test_actual_candidate_json_hydrates_without_payload_access(self):
        # This reads the locally generated JSON pins only; archive opening and
        # numeric decoding are replaced with failures while hydration runs.
        document = load_metadata()
        with patch.object(codec, "_snapshot", side_effect=AssertionError("array opened")), \
                patch.object(codec, "_decode_slot", side_effect=AssertionError("retained decode")):
            receipt = metadata_receipt(document)
        self.assertEqual(receipt["status"], "PASS_ALL_FOUR_CONTINUATION_METADATA_HYDRATIONS")
        self.assertEqual(len(receipt["capsules"]), 4)
        self.assertTrue(all(c["selected_decode_members"] == 9 for c in receipt["capsules"]))
        self.assertEqual(receipt["retained_array_paths_opened"], 0)


if __name__ == "__main__":
    unittest.main()
