# Material second-order shell interface sources

This separate continuation derives the complete quadratic matching,
Israel, and scalar junction equations in family-adapted outward Gaussian
normal charts. It leaves the frozen localized checkpoint unchanged.

Read `MATERIAL_SECOND_ORDER_JUNCTIONS.md` for equations, geometric gauge
assumptions, the registered scalar source `−4 ψ₁²`, the necessary
rescaled tensor Neumann condition in the exact reflection-symmetric
even-master class, and the limits of the result.

The equations supply boundary algebra. Evolution through shell contact,
the bulk source calculation, gauge conversion, compatible hyperbolic
boundary conditions, and a well-posedness proof remain to be completed.
No exact nonlinear family or cosmological conclusion is established.

Reproduce the exact checks with Python and the pinned SymPy dependency:

```sh
python check_material_junctions.py --output verification-normal.json
python -O check_material_junctions.py --output verification-optimized.json
python check_material_junctions.py --inject-failure
python -O check_material_junctions.py --inject-failure
```

The last two commands must fail with exit code 1 and `deliberate_failure`.
Normal and optimized runs each completed **20 grouped positive algebra
checks and 11 wrong-formula controls**. These counts are algebra checks,
not numbers of independent physical discoveries. The verifier imports
none of the previous checkpoints' scientific code.

[Technical derivation](MATERIAL_SECOND_ORDER_JUNCTIONS.md) · [Independent internal review](INDEPENDENT_REVIEW.md) · [Completed verification](VERIFICATION_RECEIPT.json). This GitHub-only supplement accompanies the [localized and formal checkpoint](../HDBLAST_SECOND_ORDER_LOCALIZED_20261014/README.md); no Zenodo edition or APS submission was made.
