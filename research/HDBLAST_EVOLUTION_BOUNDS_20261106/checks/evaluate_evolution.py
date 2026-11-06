"""Finite exact replay from immutable public inputs; no saved producer reports.

Analytic coefficient, harmonic completion, tensor-norm and energy arguments
remain the accompanying proof notes' conditional mathematical premises.
"""
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from pathlib import PurePosixPath
import re

REF = "c4b32c03243e2bcc04b748570e7c1bd5861703fa"
TREE = "8da180daec381dbe79a4dacf0a14d67bf5663625"
REPOSITORY = "https://github.com/maldonado-research/HDblast"


def need(condition, message):
    if not condition:
        raise ValueError(message)


def public_sources(root, check):
    pins = json.loads((root / "SOURCE_PINS.json").read_text(encoding="utf-8"))
    need(set(pins) == {"schema", "repository", "ref", "tree", "scope", "files"}
         and pins["schema"] == "hdblast-evolution-bounds-pins-v1", "Public pin schema mismatch")
    need(pins["ref"] == REF and pins["tree"] == TREE
         and pins["repository"] == REPOSITORY, "Public source repository/ref/tree mismatch")
    files = pins["files"]
    need(isinstance(files, list) and len(files) == 15, "Fifteen original public inputs required")
    records = {}
    for row in files:
        need(set(row) == {"bytes", "git_blob_sha1", "local", "path", "raw_url",
                          "ref", "sha256", "url"}, "Unexpected public pin keys")
        source = PurePosixPath(row["path"])
        local = PurePosixPath(row["local"])
        need(not source.is_absolute() and ".." not in source.parts
             and source.as_posix() == row["path"] and "\\" not in row["path"]
             and source.parts[0] == "research", "Unsafe original public path")
        need(local.parts[0] == "inputs" and len(local.parts) == 2
             and local.as_posix() == row["local"] and "\\" not in row["local"]
             and local.name not in records, "Unsafe or repeated input pin")
        need(row["ref"] == REF
             and row["url"] == REPOSITORY + "/blob/" + REF + "/" + source.as_posix()
             and row["raw_url"] == "https://raw.githubusercontent.com/"
                 "maldonado-research/HDblast/" + REF + "/" + source.as_posix(),
             "Immutable public URL mismatch")
        raw = (root / row["local"]).read_bytes()
        check("source_size:"+local.name, type(row["bytes"]) is int and row["bytes"] == len(raw))
        check("source_sha256:"+local.name,
              re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) is not None
              and hashlib.sha256(raw).hexdigest() == row["sha256"])
        blob = b"blob " + str(len(raw)).encode("ascii") + b"\0" + raw
        check("source_git_blob:"+local.name,
              hashlib.sha1(blob).hexdigest() == row["git_blob_sha1"])
        records[local.name] = row
    actual = {p.relative_to(root).as_posix() for p in (root / "inputs").iterdir()}
    need(actual == {row["local"] for row in files}, "Public input roster mismatch")
    return records


def load_module(root, name):
    path = root / "checks" / (name+".py")
    spec = importlib.util.spec_from_file_location("evolution_"+name, path)
    need(spec is not None and spec.loader is not None, "Pinned evaluator loader unavailable")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate(root):
    data = json.loads((root / "DATA.json").read_text(encoding="utf-8"))
    need(set(data) == {"schema", "coefficients", "gauge", "energy", "scope"}
         and data["schema"] == "HDBLAST_EVOLUTION_BOUNDS_DATA_V1", "Data schema mismatch")
    for section in ("coefficients", "gauge"):
        for key, value in data[section].items():
            if key.endswith("_decade") or key.endswith("_order"):
                need(type(value) is int and value >= 0,
                     "Nonnegative integer order/decade required: " + section + ":" + key)
    suites = {name: [] for name in ("sources", "coefficients", "gauge_symbolic", "gauge_bounds", "energy")}

    def checker(suite):
        def check(name, condition):
            need(condition, "Arithmetic failed: " + suite + ":" + name)
            need(isinstance(name, str) and name
                 and name not in {row["name"] for row in suites[suite]}, "Repeated or invalid check name")
            suites[suite].append({"name":name, "status":"PASS"})
        return check

    records = public_sources(root, checker("sources"))
    # Bind the one inherited I2 premise to its literal public statement.
    physical = (root / "inputs/PHYSICAL_CAUCHY_SECOND_JET_06E.md").read_text(encoding="utf-8")
    checker("gauge_bounds")("retained physical second-data norm premise",
        "I2=||h2||_H14+||K2||_H13+||phi2||_H14+||Pi2||_H13 <10^"
        +str(data["gauge"]["physical_I2_input_decade"])+"." in physical)
    exact = {}
    exact["coefficients"] = load_module(root, "coefficient_bounds").evaluate(
        root, data["coefficients"], checker("coefficients"))
    exact["gauge_symbolic"] = load_module(root, "gauge_symbolic").symbolic(
        checker("gauge_symbolic"), data["gauge"])
    exact["gauge_bounds"] = load_module(root, "gauge_bounds").evaluate(
        root, data["gauge"], checker("gauge_bounds"))
    exact["energy"] = load_module(root, "energy_bounds").evaluate(
        root, data["energy"], checker("energy"))
    names = [suite+":"+row["name"] for suite, rows in suites.items() for row in rows]
    return {"schema":"hdblast-evolution-bounds-v1", "status":"PASS",
        "named_check_count":len(names), "named_checks":names,
        "suite_counts":{suite:len(rows) for suite,rows in suites.items()},
        "checks":suites, "source_count":len(records), "exact":exact,
        "scope":data["scope"]}
