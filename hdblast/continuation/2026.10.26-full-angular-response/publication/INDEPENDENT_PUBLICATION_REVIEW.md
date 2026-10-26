# Independent publication review — 26 October 2026

**Source review passes. GitHub CI was still running at the last snapshot; this is not CI or merge clearance.** No blocking source finding was identified. This review made no repository or remote changes.

| Repository | Reviewed immutable head | Tree | Sole parent and current main at 07:19 UTC |
| --- | --- | --- | --- |
| HDblast | `8ccc766cb36624a7d34c7aa57356d7e7e5357eea` | `e242c877ce0fd63d9ad50694386211a308a3ee6c` | `9a45ee0deaf3496fee93491c2803cb168bb36c92` |
| Shared research site | `e11428415a89f68827079d3d67b450d80c3e7345` | `1da129b8637157872b2239290d918a1b35824a5a` | `eb97aef8ee7f72795dfffdb1668c942f8e09fcc6` |

The complete GitHub recursive trees were fetched independently, without truncation, and matched all local immutable leaf entries, including modes and blob IDs. Both PRs were open and mergeable: HDblast PR86 and shared-site PR77. Their heads and bases matched the table.

HDblast preserves 11,775 old entries exactly. Its three modified old files contain only the reviewed public insertions, and removing each exact inserted block recovers its entire base file byte for byte. All 18 other shared-site entries remain exact; its one modified page passes the same full-file inverse check. There are no deleted files. The 83 added HDblast files are confined to the 77-file scientific capsule, matching ZIP, one workflow, and four declared continuation evidence files. Every one of the 87 changed files across the two repositories matches its expected Git blob, byte count and SHA256. Existing frozen manuscripts, prior research, unrelated project files and the pre-existing symlink are preserved.

All 77 scientific files match the externally pinned manifest. All eight upstream capsule closures were independently authenticated from committed objects, including every declared member's size and SHA256. The ZIP passes CRC, has exactly the same 77 members, and every member equals the committed source bytes.

- Manifest SHA256: `a35f84119e84bf2ea8a712ded86ff79977d3ee91ac25f91ccbd202d379c859ac`.
- ZIP SHA256: `05a39514a6d05b4a5e2881fe5d6abb143f7b9f54c5aa723c22021db4a91fe5f2`.
- Replay entry-point SHA256: `96a2ed48b3cbe572d0d3c0988e0371588be8d218bdb85f59e8fb620e1a8e4c4f`.
- Initial packaged-replay receipt SHA256: `1f4aca6f8d2b0f6c1305a93bba84866f9a6e15163844348524df54e84bfb8f0f`.

The committed initial receipt is byte-identical to the completed packaged replay. It records PASS for 24 unique successful direct runs, 35 scientific golden comparisons, seven rejected source corruptions, eight unchanged upstream closures and 77 unchanged new source files. All 35 golden hashes equal their committed references; 33 retained fresh output hashes were checked directly, while the other two are temporary structure outputs checked by the authenticated control wrapper. Every direct-run log hash and every indexed fresh output's bytes were also verified. This review checks the completed replay evidence and does not rerun the numerical campaign.

The capsule README and derivation match the previously reviewed document hashes. The displayed error, integral-defect, resolvent and TT interval bounds round outward from exact rational data. The text retains the parent coefficient-existence and subsidiary assumptions, distinguishes the zero-TT global approximation from the separate nonzero actual total TT result, and distinguishes intrinsic sphere angular control from physical warp/frame conversion. It claims no pointwise C2 Einstein residual, nonlinear physical family, observed gravitational wave or cosmological origin.

The committed workflow and all four public files also passed a second independent agent's actual-object review. Workflow SHA256 is `97cb093cfbc05cb39c64d20ff38eff4cb96bc50f4fa3e7433774d5c827f9c5d3`. Its three pins, eight upstream sparse paths, capsule and ZIP, parser flags, entry-point authentication, Python 3.12.14, read-only permissions, disabled credential persistence, 40-minute timeout and unconditional artifact upload match the request. The receipt step checks 24/35/7/8/77 and prints every golden hash plus a compact pinned summary to the native job log.

At 07:19 UTC the dedicated exact-head replay, run `40284569389`, remained in progress. Several other exact-head checks had succeeded and none had failed. The shared PR had no check or workflow runs. Completed current-head CI remains a separate requirement. The extra blank line at EOF in the frozen TT certificate is a whitespace-only warning; changing it would disturb the frozen source pin without changing the science.

Machine-readable details are in `INDEPENDENT_PUBLICATION_REVIEW.json`; immutable remote trees, commits and PR/CI snapshots are in `independent-review-evidence/`.
