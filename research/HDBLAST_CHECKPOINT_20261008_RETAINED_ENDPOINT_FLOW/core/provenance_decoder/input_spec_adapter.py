"""Metadata-only hydration for the continuation candidate; no array I/O.

This module does not provide a retained-input decoding entry point. A caller that
later uses the copied codec remains responsible for its registered decode gate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

try:
    from . import exact_binary80 as codec
except ImportError:
    import exact_binary80 as codec


SELECTED_MEMBERS = (
    "k.npy", "momentum_weights.npy", "observation_eta.npy",
    "u_1.npy", "w_1.npy", "u_2.npy", "w_2.npy", "u_3.npy", "w_3.npy",
)


def hydrate_reader_spec(data: dict) -> codec.InputSpec:
    """Convert frozen JSON metadata to typed specs without opening any arrays."""
    try:
        members = []
        for raw in data["members"]:
            clean = {key: value for key, value in raw.items() if key != "payload_bytes"}
            clean["shape"] = tuple(clean["shape"])
            clean["npy_version"] = tuple(clean["npy_version"])
            members.append(codec.MemberSpec(**clean))
        fields = {key: value for key, value in data.items() if key != "members"}
        return codec.InputSpec(**fields, members=tuple(members))
    except (KeyError, TypeError, ValueError) as error:
        raise codec.CodecError("invalid candidate reader metadata") from error


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise codec.CodecError("duplicate candidate JSON key")
        result[key] = value
    return result


def load_metadata(path: str | Path | None = None) -> dict:
    """Read only the candidate JSON document, never its referenced archives."""
    selected_path = Path(path) if path is not None else Path(__file__).with_name("INPUT_SPEC.json")
    return json.loads(selected_path.read_text(encoding="utf-8"), object_pairs_hook=_unique_object)


def metadata_receipt(document: dict) -> dict:
    """Validate all four typed member rosters and resource bounds in memory."""
    limits = codec.Limits(**document["explicit_limits"])
    expected_pairs = [
        ("positive_B", "coarse"), ("positive_B", "fine"),
        ("signed_uB", "coarse"), ("signed_uB", "fine"),
    ]
    if [(c["source"], c["grid"]) for c in document["capsules"]] != expected_pairs:
        raise codec.CodecError("candidate capsule universe/order differs")
    receipts = []
    for capsule in document["capsules"]:
        spec = hydrate_reader_spec(capsule["reader_inputspec"])
        spec_digest = codec._validate_spec(spec, limits)
        selected = {member.name for member in spec.members if member.decode}
        if (len(spec.members) != 19 or selected != set(SELECTED_MEMBERS)
                or capsule["selected_decode_members"] != list(SELECTED_MEMBERS)):
            raise codec.CodecError("candidate selected or complete member roster differs")
        node_count = capsule["full_node_count"]
        for member in spec.members:
            if member.name in selected:
                expected_shape = (6,) if member.name == "observation_eta.npy" else (node_count,)
                expected_descr = "<c32" if member.name.startswith(("u_", "w_")) else "<f16"
                if member.shape != expected_shape or member.descr != expected_descr:
                    raise codec.CodecError("candidate continuation member alignment/type differs")
        receipts.append({
            "source": capsule["source"], "grid": capsule["grid"],
            "complete_members": len(spec.members), "selected_decode_members": len(selected),
            "canonical_reader_inputspec_sha256": spec_digest,
        })
    return {
        "status": "PASS_ALL_FOUR_CONTINUATION_METADATA_HYDRATIONS",
        "capsules": receipts, "retained_array_paths_opened": 0,
        "retained_values_decoded": 0, "physical_source_callbacks": 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-spec", type=Path)
    args = parser.parse_args()
    print(json.dumps(metadata_receipt(load_metadata(args.input_spec)), sort_keys=True, indent=2))


if __name__ == "__main__":
    main()
