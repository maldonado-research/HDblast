# Conditional response attribution and a prospective residual method

This separate note uses the retained coupled reduced operator to write
`U=U_hom+W`, where `W=(I-T)^-1(KF+LG)`. The complete bulk/shell source bounds
give `||W||_X<5.228e13` and scalar-mean magnitude `<1.5684e14`. The existing
large mean signs persist under this auxiliary source removal; neither summand
is asserted to be a separate physical constraint solution. These are upper
bounds, with no nonzero increment, profile or local-insignificance claim.

`followon-response/FORCED_INCREMENT_BOUND.md` gives the proof;
`followon-response/METHOD_AND_GAPS.md` describes a future continuous candidate
and complete residual ledger. A coarse total-response error certificate does
not certify resolution of the smaller increment; it does not prove that the
actual approximation error is large or that resolution is impossible. The
full waveform and continuous residual remain unfinished. This package is
separate from the immutable geometric-source capsule.

`followon-response-review/` contains the independent internal review and its
accepted receipt. `PUBLICATION_PINS.json` pins all other 23 files, including
the verifier and all 21 approved scientific/review assets. Its own external
SHA-256 is hardcoded in the dedicated GitHub workflow. After obtaining that
trusted value, choose a new output directory outside this package:

```sh
python -I -S -B verify_followon.py --pins-sha256 "$verified_followon_pins_sha256" --output-dir /tmp/hdblast-response-attribution-verification
```

The standard-library helper authenticates the complete 24-file inventory and
runs the producer and independent checker in normal and optimized modes:
four tiny arithmetic processes with complete typed/raw golden agreement,
no new mutation controls and no repeat of the expensive core scientific
replay. The core capsule's 45/14/84/24 process/check counts remain separate.
