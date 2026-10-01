# Next tests (ranked)

1. **Derived energy-exchange channel with an effective Y of ~1.6 or more (B1 -> B3 link).** This is the decisive open question. It needs:
   - a solver that runs beyond the bounded-chart limit (H0tau ~7.4), so that b_cut runs can reach the registered Delta ln a = 0.5 plateau (this would turn B1's descriptive r >= 94 into a registered class);
   - a scan of b in 1e-3..1e-2 (lambda_c 5-12) with backreaction;
   - a self-consistent production model (chi mode functions on the actual shell trajectory, with an adiabaticity-based completion criterion).

   Pre-register the effective-Y comparison against B3's r(Y).
2. **Residual-vacuum recollapse / fine tuning.** Every passing plateau ends within 0.1-0.5 e-folds because Lambda_res < 0 at d*_M8. Two tasks:
   - Run 5D at d = d*_exact = -3.10415 (B4) to check whether the radiation era lengthens as predicted.
   - Build an AdS-sliced static solver so Lambda_res(d) beyond d* needs no extrapolation.

   Any proposal must address the 7e-7 to 7e-80 tolerance (the cosmological-constant problem). In the scanned range (Y = 0.5-5) a larger Y shortens, rather than lengthens, the expansion left after the plateau.
3. **delta = 1e-3 gauge past the light-cone crossing (B2).**
   - Add a second chart stage or a proper-time shell gauge.
   - Use kappa = 10 for all pre-runs, and redo the dc = 1e-4 main pair at a non-artifact x_c (18.7).
   - Test the sensitivity to x_c explicitly, since x_c is not ruled out.
   - Optionally, rescale the friction closure to a fixed Y rho_b in Hubble units (~7.8) so the roll finishes before the crossing.
4. **Consolidate B3 at delta = 0.1.**
   - Run the dc = 1e-4 repeats for Y = 3 and 1.5 with the three-level chain.
   - Use a chart whose shell lapse stays O(1) (this tests the lapse-causation hypothesis), and extend the static table beyond T ~ 10.75.
   - Converge Y = 5.
5. **Calibration rule C2.** Replace the registered C2 with an eigenmode seed, or with the offset-free derivative estimator (pre-registered), so it can pass as written.
6. **Reduced model.** Model the final stop (phi_b 0.95 -> 1.036) with the 5D near-endpoint linear problem or a KK/CFT correction before using any EFT r prediction.
7. **Re-verify the N_eff thresholds and the r -> Delta N_eff mapping** beyond snippet level (these are inherited from A1 and not re-checked here).
8. **Repeat the B3 scan and B1 at delta = 1e-3**, once item 3 is solved.
