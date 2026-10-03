# Bundled saved-trajectory inputs

These are the eleven numeric B3 timeseries NPZ files and matching summaries needed for the nine audited histories, including restart parents. Original states, unrelated runs, chats and working files are excluded. Read NPZ with allow_pickle=False. The data were produced in September, not evolved anew in this checkpoint.

Baseline: `maldonado-research/HDblast-archive@8f67197b730d4e1c43554b86f224c29cc72629eb`; original paths are under `research/HDBLAST_CHECKPOINT_20260930/B3_Y_scan/`. Each summary was already available in the public September checkpoint. SHA-256 values of every copied binary and summary were independently matched to the executed audit's input provenance before bundling. The accompanying B3_RESULTS.json is the byte-identical baseline saved-point report.

The frozen original summaries preserve producer labels. Final audit reliability remains controlling: Y3 conditional, Y5 inconclusive/unreliable; do not upgrade either by reading a producer PASS label.

| Timeseries file | SHA-256 |
|---|---|
| runs/main/main_Y2_dc1e-2_dzf5e-4_timeseries.npz | `8744bc182c0c64651cd1666c64ccb7cdd7aba7ba98fbd39a3bc42808d840b44d` |
| runs/fine/fine_Y2_dc1e-2_dzf2.5e-4_timeseries.npz | `74526055ee928427b690300a5a8d3bd639e68c2b0a459788b0fec4c481622302` |
| runs/main/main_Y2_dc1e-4_dzf5e-4_timeseries.npz | `8b9c9202692436d5d1f09655ef6439f04d9de49614ae8746753817060dfba46a` |
| runs/fine/fine_Y2_dc1e-4_dzf2.5e-4_timeseries.npz | `8478335b2f96cdc3ee1ab5a700d6e37bde77a314b9c176da4f47051b04f7313a` |
| runs/main/main_Y3_dc1e-2_dzf5e-4_timeseries.npz | `c4ae136ef41d1248cab2939b0fbc39dcecee2b3c0f4a6a1e99da341c017c2151` |
| runs/fine/fine_Y3_dc1e-2_dzf1.25e-4_timeseries.npz | `fe8693c76a008fd0b716b2d6841bf000a50d8aac858301fc35a35f5e870b2757` |
| runs/main/main_Y5_dc1e-2_dzf5e-4_cfl0.25_timeseries.npz | `d97fc03c3b17bcba22cc3f620433e971efda1660ba05798e62b2b2923f3ad371` |
| runs/fine/fine_Y5_dc1e-2_dzf1.25e-4_cfl0.25_timeseries.npz | `1f8e4d43b55bf668fa368d56cf7432d2e16c8ab0664db6f3650684b4ec24e7ef` |
| runs/main/main_Y0.5_dc1e-2_dzf5e-4_timeseries.npz | `5b40f1131bf26e7ad92be7071b4633e7f23a9743880878d147610a1a74339416` |
| runs/main/main_Y0.7_dc1e-2_dzf5e-4_timeseries.npz | `269efc673d2d2a5772355eae2ec1cb940a99379ceff0bf5ad6f44b73950412ba` |
| runs/main/main_Y1.5_dc1e-2_dzf5e-4_timeseries.npz | `13b4c3279b131a988408bbe0fdb072dd4aecdddb02dcfd62f9eb7b74f9742b22` |
