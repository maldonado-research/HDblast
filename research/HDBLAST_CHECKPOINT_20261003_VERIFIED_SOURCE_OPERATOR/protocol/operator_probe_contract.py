"""Stdlib validator for narrow registered moment/work output; no source calls.

Bootstrap/source proof/runtime/resource authentication belongs to the outer
registered driver. This validator checks exact output shape, membership,
real-target identities, containment compatibility and complete radii.
"""
from __future__ import annotations
from fractions import Fraction as Q
import json
import re

SOURCES = ("positive_B", "signed_uB")
MOMENTA = (Q(0), Q(1, 2**40), Q(1, 2**12), Q(1, 4),
           Q(1), Q(16), Q(64), Q(128), Q(256))
CHANNELS = ("M0", "Mexp", "Mu")
GATE = Q(1, 10**26)
SCOPE = "UNIFORM_ANALYTIC_MODEL_PLUS_FIXED_RATIONAL_MOMENT_PROBES"
ROW_KEYS = {"source", "momentum", "interval", "moments", "total_absolute_radii"}
TOP_KEYS = {"schema_version", "scope", "configuration", "panel_rows",
            "whole_rows", "source_work_panel_rows", "source_work_whole_rows",
            "archive_arrays_decoded"}
CONFIGURATIONS = ("PRIMARY_ARB256_SOURCE24_PHASE96", "INDEPENDENT_DYADIC512_SOURCE24_ODE96")
WORK_KEYS = {"source", "interval", "moment", "total_absolute_radius"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rat(value):
    value = Q(value)
    return f"{value.numerator}/{value.denominator}"


def fraction(text):
    require(type(text) is str and len(text) <= 4096 and
            re.fullmatch(r"-?(?:0|[1-9][0-9]*)/[1-9][0-9]*", text),
            "Canonical exact rational required")
    value = Q(text)
    require(rat(value) == text, "Reduced canonical rational required")
    return value


def interval(box):
    require(type(box) is dict and set(box) == {"lo", "hi"},
            "Exact lower/upper endpoints required")
    lo, hi = fraction(box["lo"]), fraction(box["hi"])
    require(lo <= hi, "Reversed interval")
    return lo, hi


def encode(lo, hi):
    require(lo <= hi, "Reversed interval")
    return {"lo": rat(lo), "hi": rat(hi)}


def complex_interval(box):
    require(type(box) is dict and set(box) == {"real", "imag"},
            "Two component enclosures required")
    real, imag = interval(box["real"]), interval(box["imag"])
    return real, imag


def absolute_radius(box):
    real, imag = complex_interval(box)
    return (real[1]-real[0]+imag[1]-imag[0])/2


def overlap(first, second):
    return max(first[0], second[0]) <= min(first[1], second[1])


def validate_chronology(value):
    require(type(value) is dict and set(value) == {
        "new_target_evaluations_before_public_verification",
        "prior_known_outcomes", "not_blind_confirmation"},
        "Canonical chronology schema required")
    marker = value["new_target_evaluations_before_public_verification"]
    require(type(marker) is dict and set(marker) == {"primary", "independent"},
            "Exactly two named real-source routes required")
    require(all(type(v) is int and v == 0 for v in marker.values()),
            "Counts must be exact integer zeros")
    require(value["not_blind_confirmation"] is True and value["prior_known_outcomes"] == [
        "ORIGINAL_METRIC_ENDPOINT_29_AND_REFINEMENT_30_FAILURES",
        "COMPLETED_EMPIRICAL_ACTIVE_SOURCE_ATTRIBUTION",
        "ORIGINAL_REGISTRATION_MARKER_FAILURE_AND_DISCLOSED_REPAIR"],
        "Disclosed known-outcome chronology required")


def panel_interval(panel):
    start = Q(-9, 2)+Q(panel, 64)
    return [rat(start), rat(start+Q(1,64))]


def validate_row(row, source, momentum, endpoints, panel=None):
    require(type(row) is dict and set(row) == ROW_KEYS | ({"panel"} if panel is not None else set()),
            "Exact row schema required")
    require(row["source"] == source and row["momentum"] == rat(momentum) and
            row["interval"] == endpoints, "Row target/ordering mismatch")
    if panel is not None:
        require(type(row["panel"]) is int and row["panel"] == panel,
                "Exact panel index required")
    require(type(row["moments"]) is dict and set(row["moments"]) == set(CHANNELS) and
            type(row["total_absolute_radii"]) is dict and
            set(row["total_absolute_radii"]) == set(CHANNELS),
            "Complete moment/radius inventory required")
    consistency = True
    for channel in CHANNELS:
        radius = absolute_radius(row["moments"][channel])
        require(fraction(row["total_absolute_radii"][channel]) == radius,
                "Reported complete radius differs from exact endpoints")
    m0 = complex_interval(row["moments"]["M0"])
    require(m0[1] == (Q(0), Q(0)), "M0 must have exactly zero imaginary component")
    if momentum == 0:
        mexp = complex_interval(row["moments"]["Mexp"])
        mu = complex_interval(row["moments"]["Mu"])
        require(mexp[1] == mu[1] == (Q(0), Q(0)),
                "Zero-momentum moments must be real")
        consistency = overlap(m0[0], mexp[0])
    return consistency


def validate_payload(payload):
    require(type(payload) is dict and set(payload) == TOP_KEYS, "Exact payload schema required")
    require(type(payload["schema_version"]) is int and payload["schema_version"] == 1 and
            payload["scope"] == SCOPE and payload["configuration"] in CONFIGURATIONS,
            "Payload scope/configuration mismatch")
    require(type(payload["archive_arrays_decoded"]) is int and
            payload["archive_arrays_decoded"] == 0,
            "Archived physical arrays are forbidden in this checkpoint")
    panel_rows, whole_rows = payload["panel_rows"], payload["whole_rows"]
    require(type(panel_rows) is list and len(panel_rows) == 1152 and
            type(whole_rows) is list and len(whole_rows) == 18,
            "Exact 1152-panel/18-whole row counts required")
    consistency = True
    index = 0
    for source in SOURCES:
        for panel in range(64):
            m0_boxes = []
            for momentum in MOMENTA:
                row = panel_rows[index]
                consistency = validate_row(row, source, momentum, panel_interval(panel), panel) and consistency
                m0_boxes.append(complex_interval(row["moments"]["M0"])[0])
                index += 1
            consistency = max(b[0] for b in m0_boxes) <= min(b[1] for b in m0_boxes) and consistency
    index = 0
    radius_pass = True
    worst_radius = Q(0)
    for source in SOURCES:
        m0_boxes = []
        for momentum in MOMENTA:
            row = whole_rows[index]
            consistency = validate_row(row, source, momentum, ["-9/2","-7/2"]) and consistency
            m0 = complex_interval(row["moments"]["M0"])[0]
            m0_boxes.append(m0)
            # M0 is momentum independent and additive. Different correct
            # enclosures need only overlap; containment of the wider sum is
            # not required when coefficient dependencies are retained.
            selected = [r for r in panel_rows if r["source"] == source and
                        r["momentum"] == rat(momentum)]
            boxes = [complex_interval(r["moments"]["M0"])[0] for r in selected]
            additive = (sum(b[0] for b in boxes),sum(b[1] for b in boxes))
            consistency = overlap(m0,additive) and consistency
            for channel in CHANNELS:
                radius = absolute_radius(row["moments"][channel])
                worst_radius = max(worst_radius,radius)
                radius_pass = radius <= GATE and radius_pass
            index += 1
        consistency = max(b[0] for b in m0_boxes) <= min(b[1] for b in m0_boxes) and consistency
    work_panels, work_whole = payload["source_work_panel_rows"], payload["source_work_whole_rows"]
    require(type(work_panels) is list and len(work_panels) == 128 and
            type(work_whole) is list and len(work_whole) == 2,
            "Exact 128-panel/2-whole source-work row counts required")
    work_index = 0
    for source in SOURCES:
        selected = []
        for panel in range(64):
            row = work_panels[work_index]
            validate_work_row(row, source, panel_interval(panel), panel)
            selected.append(complex_interval(row["moment"])[0])
            work_index += 1
        row = work_whole[SOURCES.index(source)]
        validate_work_row(row, source, ["-9/2", "-7/2"])
        whole = complex_interval(row["moment"])[0]
        additive = (sum(b[0] for b in selected), sum(b[1] for b in selected))
        consistency = overlap(whole, additive) and consistency
        radius = absolute_radius(row["moment"])
        worst_radius = max(worst_radius, radius)
        radius_pass = radius <= GATE and radius_pass
    return {"status": "CERTIFICATE_CONSISTENCY_FAILURE" if not consistency else
            "PASS_REGISTERED_PROBE_OUTPUT_VALIDATION" if radius_pass else "UNRESOLVED_PROBE_RADIUS",
            "uniform_analytic_proof_validated_by_this_module":False,
            "production_freeze_authorized":False,
            "configuration":payload["configuration"],
            "panel_rows":1152,"whole_rows":18,
            "source_work_panel_rows":128,"source_work_whole_rows":2,
            "worst_whole_interval_absolute_radius":rat(worst_radius),
            "gate":rat(GATE)}


def validate_work_row(row, source, endpoints, panel=None):
    require(type(row) is dict and set(row) == WORK_KEYS | ({"panel"} if panel is not None else set()),
            "Exact source-work row schema required")
    require(row["source"] == source and row["interval"] == endpoints,
            "Source-work target/ordering mismatch")
    if panel is not None:
        require(type(row["panel"]) is int and row["panel"] == panel,
                "Exact source-work panel index required")
    _, imag = complex_interval(row["moment"])
    require(imag == (Q(0), Q(0)), "Lg source work must be exactly real")
    require(fraction(row["total_absolute_radius"]) == absolute_radius(row["moment"]),
            "Source-work radius differs from exact endpoints")


def validate_pair(primary, independent):
    results = [validate_payload(primary), validate_payload(independent)]
    require(primary["configuration"] == CONFIGURATIONS[0] and
            independent["configuration"] == CONFIGURATIONS[1],
            "Exact primary/independent route configurations required")
    consistency = all(r["status"] != "CERTIFICATE_CONSISTENCY_FAILURE" for r in results)
    radii_pass = all(r["status"] == "PASS_REGISTERED_PROBE_OUTPUT_VALIDATION" for r in results)
    max_hull_radius = Q(0)
    for field in ("panel_rows", "whole_rows", "source_work_panel_rows", "source_work_whole_rows"):
        for p_row, i_row in zip(primary[field], independent[field]):
            boxes = [(p_row["moments"][channel], i_row["moments"][channel])
                     for channel in CHANNELS] if "moments" in p_row else [(p_row["moment"], i_row["moment"])]
            for p_box, i_box in boxes:
                p_re, p_im = complex_interval(p_box)
                i_re, i_im = complex_interval(i_box)
                consistency = overlap(p_re, i_re) and overlap(p_im, i_im) and consistency
                hull_radius = (max(p_re[1],i_re[1])-min(p_re[0],i_re[0])+
                               max(p_im[1],i_im[1])-min(p_im[0],i_im[0]))/2
                if "whole" in field:
                    max_hull_radius = max(max_hull_radius, hull_radius)
                    radii_pass = hull_radius <= GATE and radii_pass
    return {"status":"CERTIFICATE_CONSISTENCY_FAILURE" if not consistency else
            "PASS_BOTH_REGISTERED_PROBE_OUTPUTS" if radii_pass else "UNRESOLVED_PROBE_RADIUS",
            "uniform_analytic_proof_validated_by_this_module":False,
            "outer_bootstrap_proof_and_resource_checks_required":True,
            "archive_array_decoding_allowed":False,
            "max_whole_interval_hull_absolute_radius":rat(max_hull_radius),
            "gate":rat(GATE),"route_results":results}


def load_json_no_duplicates(raw):
    def unique(pairs):
        result = {}
        for key,value in pairs:
            require(key not in result,"Duplicate JSON key")
            result[key] = value
        return result
    def invalid(value):
        raise ValueError("Nonfinite JSON literal")
    return json.loads(raw,object_pairs_hook=unique,parse_constant=invalid)
