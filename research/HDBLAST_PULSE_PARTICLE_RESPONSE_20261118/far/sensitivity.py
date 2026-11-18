"""Check the isolated positive response under changes of state, analyzer and integration accuracy."""
import importlib.util
import json
from pathlib import Path
import tempfile


def calculate():
    far = Path(__file__).resolve().parent
    classical = far.parent / "classical"
    spec = importlib.util.spec_from_file_location("rolling_solver", classical / "rolling_modes.py")
    solver = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(solver)
    tests = [
        ("leading_analyzer", .8, 0, .25, 2e-9),
        ("earlier_start", .5, 2, .25, 2e-9),
        ("refined", .8, 2, .125, 2e-10),
    ]
    results = {}
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        grid = root / "fine"
        grid.mkdir()
        for item in ("checkpoint.npz", "quiet"):
            (grid / item).symlink_to(classical / "runs/fine" / item)
        (grid / "0.001").symlink_to(far / "runs/fine/0.001")
        for name, start, order, spacing, tolerance in tests:
            quiet, pulse = [
                solver.compute(root, "fine", branch, start, spacing, tolerance, order)
                for branch in ("quiet", "0.001")
            ]
            results[name] = {
                t: 100 * (row["comoving_N_8"] / quiet["ends"][t]["comoving_N_8"] - 1)
                for t, row in pulse["ends"].items()
            }
    return results


if __name__ == "__main__":
    print(json.dumps(calculate(), indent=2))
