#!/usr/bin/env python3
"""Read source ASTs only; never import or execute either physical producer."""
from __future__ import annotations

import argparse
import ast
import copy
import hashlib
import json
from pathlib import Path
import sys


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def dump(node):
    return ast.dump(node, include_attributes=False)


def statement(source):
    return ast.parse(source).body[0]


def check(original_source, revised_source):
    first, second = ast.parse(original_source), ast.parse(revised_source)
    original = {n.name: n for n in first.body if isinstance(n, ast.FunctionDef)}
    revised = {n.name: n for n in second.body if isinstance(n, ast.FunctionDef)}
    require(set(original) == set(revised), "Producer function inventory changed")
    unchanged = []
    for name in original:
        if name not in ("sha", "provenance_gate", "evolve"):
            require(ast.get_source_segment(original_source, original[name]) ==
                    ast.get_source_segment(revised_source, revised[name]),
                    "Scientific function bytes changed: " + name)
            unchanged.append(name)

    added_imports = [n for n in second.body if isinstance(n, ast.ImportFrom)
                     and n.module == "stream_npz"]
    require(len(added_imports) == 1 and dump(added_imports[0]) ==
            dump(statement("from stream_npz import StreamingNpz, file_sha256")),
            "Unexpected archive helper import")
    old_top = [n for n in first.body if not isinstance(n, ast.FunctionDef)]
    new_top = [n for n in second.body if not isinstance(n, ast.FunctionDef)
               and n is not added_imports[0]]
    require(dump(ast.Module(body=old_top, type_ignores=[])) ==
            dump(ast.Module(body=new_top, type_ignores=[])),
            "Top-level declarations, gates or entrypoint changed")

    provenance = copy.deepcopy(revised["provenance_gate"])
    removed = 0
    for node in ast.walk(provenance):
        if (isinstance(node, ast.Set) and
                all(isinstance(x, ast.Constant) for x in node.elts) and
                {x.value for x in node.elts} ==
                {"forced_metric.py", "metric_wkb.py", "stable_baselines.py", "stream_npz.py"}):
            node.elts = [x for x in node.elts if x.value != "stream_npz.py"]
            removed += 1
    require(removed == 1 and dump(provenance) == dump(original["provenance_gate"]),
            "Provenance change exceeds the additional hash-pinned archive module")
    expected_sha = statement("def sha(path):\n return file_sha256(path)")
    require(dump(revised["sha"]) == dump(expected_sha),
            "SHA helper exceeds the bounded file hash substitution")

    old_evolve, new_evolve = copy.deepcopy(original["evolve"]), copy.deepcopy(revised["evolve"])
    contexts = [n for n in new_evolve.body if isinstance(n, ast.With)]
    require(len(contexts) == 1, "Exactly one archive close-on-exit wrapper is required")
    context = contexts[0]
    expected_context = statement("with StreamingNpz(path) as archive:\n pass")
    require(len(context.items) == 1 and
            dump(context.items[0]) == dump(expected_context.items[0]),
            "Unexpected archive context expression")
    flat = []
    for node in new_evolve.body:
        flat.extend(node.body if node is context else [node])
    new_evolve.body = flat
    removals = {"original": [], "revised": []}
    path_assignment = dump(statement('path = output_dir/f"metric_modes_{source_id}_{name}.npz"'))
    old_archive = dump(statement('archive = {"k": k, "momentum_weights": weights}'))
    permitted_expressions = {
        "original": {
            dump(statement("np.savez_compressed(path,**archive)")): "deferred archive write",
            dump(statement('active_arrays.update({"archived_"+key: value for key,value in archive.items()})')):
                "duplicate archival lifetime references",
        },
        "revised": {
            dump(statement('archive.update({"k": k, "momentum_weights": weights})')):
                "immediate initial standard-NPY member writes",
            dump(statement("archive.close()")): "close before SHA and resource check",
        },
    }
    permitted_deletion = dump(statement("del contact_archive, integrands"))

    def normalize(node, route):
        signature = dump(node)
        if isinstance(node, ast.Assign) and signature == path_assignment:
            removals[route].append("relocated output path")
            return None
        if route == "original" and isinstance(node, ast.Assign) and signature == old_archive:
            removals[route].append("in-memory archive dictionary construction")
            return None
        if isinstance(node, ast.Expr) and signature in permitted_expressions[route]:
            removals[route].append(permitted_expressions[route][signature])
            return None
        if route == "revised" and isinstance(node, ast.Delete) and signature == permitted_deletion:
            removals[route].append("release already written observation arrays")
            return None
        for field, value in ast.iter_fields(node):
            if isinstance(value, list):
                kept = []
                for child in value:
                    revised_child = normalize(child, route) if isinstance(child, ast.stmt) else child
                    if revised_child is not None:
                        kept.append(revised_child)
                setattr(node, field, kept)
        return node

    old_evolve = normalize(old_evolve, "original")
    new_evolve = normalize(new_evolve, "revised")
    require(len(removals["original"]) == len(removals["revised"]) == 4,
            "Permitted archive-operation inventory differs")
    require(dump(old_evolve) == dump(new_evolve),
            "Remaining evolution, observation, archive schema or arithmetic AST differs")
    return {
        "unchanged_function_count": len(unchanged), "unchanged_functions": unchanged,
        "declarations_gates_and_entrypoint_unchanged": True,
        "provenance_change_only_adds_archive_module": True,
        "permitted_archival_removals": removals,
        "remaining_evolve_AST_identical": True,
        "normalized_evolve_AST_sha256": hashlib.sha256(dump(old_evolve).encode()).hexdigest(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", type=Path, required=True)
    parser.add_argument("--revised", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    require(not args.output.exists(), "Fresh reviewer receipt is required")
    first, second = args.original.read_text(), args.revised.read_text()
    result = check(first, second)
    mutations = [
        ("source amplitude", 'EPSILON = LD("0.0001")', 'EPSILON = LD("0.0002")'),
        ("transport sign", "u+drift*w-forcing_u,E*w-forcing_w", "u+drift*w+forcing_u,E*w-forcing_w"),
        ("Wronskian gate comparison", "amplitude_W <= WRONSKIAN_GATE", "amplitude_W < WRONSKIAN_GATE"),
        ("raw history member", '"history_baseline_contact": history_baseline_contact',
         '"history_baseline_contact_removed": history_baseline_contact'),
        ("saved raw mode copy", 'archive["u_"+str(index)],archive["w_"+str(index)] = u.copy(),w.copy()',
         'archive["u_"+str(index)],archive["w_"+str(index)] = u.real.copy(),w.copy()'),
        ("missing archive closure", "        archive.close()\n", ""),
    ]
    rejected = []
    for name, target, replacement in mutations:
        require(second.count(target) == 1, "Ambiguous mutation: " + name)
        try:
            check(first, second.replace(target, replacement, 1))
        except RuntimeError:
            rejected.append(name)
        else:
            raise RuntimeError("Scientific/archive mutation survived: " + name)
    result.update(status="PASS_STATIC_INVARIANCE", physical_evaluations=0,
                  producer_imports_or_execution=0, python_optimization=sys.flags.optimize,
                  mutations_rejected=rejected,
                  source_pins={"original": {"path": str(args.original),
                                            "sha256": hashlib.sha256(args.original.read_bytes()).hexdigest()},
                               "revised": {"path": str(args.revised),
                                           "sha256": hashlib.sha256(args.revised.read_bytes()).hexdigest()}},
                  scope="Static numeric and archive-schema invariance only; no runtime memory or scientific PASS inference.")
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    print(json.dumps({"status": result["status"], "unchanged_functions": len(result["unchanged_functions"]),
                      "mutations": len(rejected), "physical_evaluations": 0}))


if __name__ == "__main__":
    main()
