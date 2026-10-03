# Scientific input provenance

The four files in A1/ are exact copies from Ricardo Maldonado's private HDBLAST research archive, at commit 8f67197b730d4e1c43554b86f224c29cc72629eb, directory research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/runs/main/. Their original filenames, byte counts, Git blob IDs and SHA-256 digests are pinned in ../INPUTS.json.

They are the source-free Y=0 tuned delta=0.1 histories at two radial resolutions. The original archive's later non-finite stopping behavior is preserved. This experiment uses the declared earlier interval; it does not repair or extend those classical evolutions.

independent_geometry.json contains the fine history's 305 shell records, with columns [H0*tau, ln_a, H/H0, phi_b, v/H0], extracted from the same NPZ for the independent JavaScript implementation. It contains no fit to the new quantum result. The primary Python implementation reads the NPZ directly.

The original bulk snapshots remain inside the copied NPZ inputs. No raw chats, personal account metadata or unrelated project files are included.
