#!/usr/bin/env python3
"""Exact rational checks of the primary wrapper's fabricated-only route."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import sys

from flint import ctx

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import primary_route as route


def assert_real_contains(encoded, target):
    assert encoded["imag"] == {"lo": "0/1", "hi": "0/1"}
    assert F(encoded["real"]["lo"]) <= target <= F(encoded["real"]["hi"])


def check():
    data, budget = route.run_fabricated()
    assert [len(data[key]) for key in ("panel_rows", "whole_rows", "source_work_panel_rows", "source_work_whole_rows")] == [1152, 18, 128, 2]
    assert budget["status"] == "PASS_CERTIFIED_OPERATOR_WIDTH_GATE"
    assert F(budget["maximum_whole_complete_L1_radius"]) <= F(1, 10**26)
    h, r, rw = F(1, 128), route.v.SOURCE_REMAINDER, route.v.SOURCE_WORK_REMAINDER
    etail = F(2 * 4096 * 8**97, factorial(97))
    qtail = F(4 * 4096 * 8**97, factorial(98))
    coefficient_l1 = {}
    for source in route.SOURCES:
        def authorize(name, detail):
            assert name == "fabricated_source"
        with ctx.workprec(256):
            p = route.fabricated_bundle(F(-9, 2) + h, source, authorize)["g"]
            coefficient_l1[source] = route.polynomial_l1(p)
        sign = 1 if source == "positive_B" else -1
        coefficients = [F(sign**j, (j + 1) * 128**j) for j in range(25)]
        a = sum((coefficients[j] * F(2, j + 1) for j in range(0, 25, 2)), F(0))
        b = sum((coefficients[j] * F(2, j + 2) for j in range(1, 25, 2)), F(0))
        whole_m0 = 64 * h * a
        whole_mu = whole_m0 / 2 - 64 * h * h * b
        for row in data["whole_rows"]:
            if row["source"] != source:
                continue
            assert_real_contains(row["moments"]["M0"], whole_m0)
            if F(row["momentum"]) == 0:
                assert_real_contains(row["moments"]["Mexp"], whole_m0)
                assert_real_contains(row["moments"]["Mu"], whole_mu)
        for row in data["panel_rows"]:
            if row["source"] == source and F(row["momentum"]) == 0:
                assert_real_contains(row["moments"]["M0"], h * a)
                assert_real_contains(row["moments"]["Mexp"], h * a)
                assert_real_contains(row["moments"]["Mu"], h * h * (a - b))
        for row in data["source_work_whole_rows"]:
            if row["source"] == source:
                assert_real_contains(row["moment"], whole_m0 / 4)
        for row in data["source_work_panel_rows"]:
            if row["source"] == source:
                assert_real_contains(row["moment"], h * a / 4)
    for row in budget["panel_rows"]:
        l1 = coefficient_l1[row["source"]]
        k = F(row["momentum"])
        source_expected = {"M0": 2 * h * r, "Mexp": 2 * h * r, "Mu": 2 * h * h * r}
        phase_expected = {"M0": F(0), "Mexp": h * l1 * etail if k else F(0), "Mu": h * h * l1 * qtail if k else F(0)}
        for key, target in source_expected.items():
            assert F(row["source_model_disk_radius"][key]) == target
        for key, target in phase_expected.items():
            assert F(row["phase_model_disk_radius"][key]) == target
    for row in budget["whole_rows"]:
        l1, k = coefficient_l1[row["source"]], F(row["momentum"])
        assert {key: F(value) for key, value in row["source_model_disk_radius"].items()} == {"M0": r, "Mexp": r, "Mu": r / 2}
        expected_exp = 64 * h * l1 * etail if k else F(0)
        expected_mu = 64 * h * h * l1 * qtail if k else F(0)
        if k:
            for j in range(64):
                distance = F(63 - j, 64)
                if 2 * k * distance <= 1:
                    expected_mu += distance * F(3, factorial(98)) * 2 * h * l1
        assert F(row["phase_model_disk_radius"]["Mexp"]) == expected_exp
        assert F(row["phase_model_disk_radius"]["Mu"]) == expected_mu
        assert F(row["phase_model_disk_radius"]["M0"]) == 0
    for row in budget["source_work_panel_rows"]:
        assert F(row["source_model_radius"]) == 2 * h * rw
        assert F(row["phase_model_radius"]) == 0
    for row in budget["source_work_whole_rows"]:
        assert F(row["source_model_radius"]) == rw
        assert F(row["phase_model_radius"]) == 0
    print("FABRICATED_WRAPPER_EXACT_MOMENT_AND_BUDGET_CHECKS_PASS")
    print("rows: 1170 moment rows, 130 work rows; physical callback refused by fabricated entry")


if __name__ == "__main__":
    check()
