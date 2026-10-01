/**
 * C3: prescribed positive-sech-squared oscillator calibration, 2026-10-01.
 * Run: node code/mode_controls.mjs
 * Pure ECMAScript; no dependencies, network access, or filesystem writes.
 *
 * Registration was read before this first execution:
 * maldonado-research/HDblast-archive @
 * ceaeb64adf76b901ded8fa7d1bf1c913a29d22e4,
 * research/HDBLAST_CHECKPOINT_20261001/REGISTRATION.md
 * registration blob: 0e2861a35d94e48e6b837ee420ecb5444d9387e4.
 *
 * The Hamiltonian work check is exploratory, additional to registered C3.
 * This prescribed-history control cannot establish a physical HDBLAST source,
 * backreaction, thermalisation, a coupled 5D solution, or new mathematics.
 */
const LAMBDAS = [0, 0.5, 1, 1.25, 2];
const OMEGAS = [0.5, 1, 2];
const STEPS = [0.02, 0.01, 0.005];
const WINDOWS = [12, 16];
const REGISTRATION = {
  repository: "maldonado-research/HDblast-archive",
  commit: "ceaeb64adf76b901ded8fa7d1bf1c913a29d22e4",
  path: "research/HDBLAST_CHECKPOINT_20261001/REGISTRATION.md",
  blob: "0e2861a35d94e48e6b837ee420ecb5444d9387e4",
  nonzeroMinExact: 1e-8,
  absoluteTolerance: 2e-6,
  relativeTolerance: 0.005,
  nullTolerance: 1e-8,
  fineWronskianTolerance: 1e-6
};
const SOURCE_REFERENCE = {
  repository: "maldonado-research/HDblast-archive",
  commit: "8f67197b730d4e1c43554b86f224c29cc72629eb",
  path: "[private archive reference; path withheld]",
  blob: "17a35192c539412e36981921bec5f1c49425cac7"
};
const ROW_COLUMNS = [
  "lambda", "omega", "T", "requested_h", "steps",
  "exact_n", "n", "abs_n_error", "relative_n_error",
  "alpha_re", "alpha_im", "beta_re", "beta_im",
  "wronskian", "abs_wronskian_error", "direct_wronskian",
  "E_initial", "E_final", "delta_E", "integrated_driver_work",
  "work_minus_delta_E", "omega_times_n",
  "delta_E_minus_omega_times_n", "occupation_pass", "wronskian_pass"
];
function frequencySquared(t, lambda, omega) {
  return omega * omega + lambda * (lambda + 1) / Math.cosh(t) ** 2;
}
function derivative(t, y, lambda, omega) {
  const sech2 = 1 / Math.cosh(t) ** 2;
  const h = lambda * (lambda + 1);
  const q = omega * omega + h * sech2;
  const dq = -2 * h * sech2 * Math.tanh(t);
  return [y[2], y[3], -q * y[0], -q * y[1],
    0.5 * dq * (y[0] ** 2 + y[1] ** 2)];
}
function shifted(y, k, scale) {
  return y.map((value, j) => value + scale * k[j]);
}
function rk4(t, y, h, lambda, omega) {
  const k1 = derivative(t, y, lambda, omega);
  const k2 = derivative(t + h / 2, shifted(y, k1, h / 2), lambda, omega);
  const k3 = derivative(t + h / 2, shifted(y, k2, h / 2), lambda, omega);
  const k4 = derivative(t + h, shifted(y, k3, h), lambda, omega);
  return y.map((value, j) =>
    value + h / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]));
}
function energy(t, y, lambda, omega) {
  return 0.5 * (y[2] ** 2 + y[3] ** 2
    + frequencySquared(t, lambda, omega) * (y[0] ** 2 + y[1] ** 2));
}
function rotate(re, im, phase, scale) {
  const c = Math.cos(phase), s = Math.sin(phase);
  return [scale * (re * c - im * s), scale * (re * s + im * c)];
}
function coefficients(t, y, omega) {
  // u = (alpha exp(-iwt) + beta exp(+iwt))/sqrt(2w).
  // u'/(iw) has real part v_im/w and imaginary part -v_re/w.
  return {
    alpha: rotate(y[0] - y[3] / omega, y[1] + y[2] / omega,
      omega * t, Math.sqrt(omega / 2)),
    beta: rotate(y[0] + y[3] / omega, y[1] - y[2] / omega,
      -omega * t, Math.sqrt(omega / 2))
  };
}
function magnitudeSquared(z) { return z[0] ** 2 + z[1] ** 2; }
function exactOccupation(lambda, omega) {
  // Analytic reflectionless limit, instead of floating-point sin(pi*N).
  if (Number.isInteger(lambda)) return 0;
  return Math.sin(Math.PI * lambda) ** 2 / Math.sinh(Math.PI * omega) ** 2;
}
function tolerance(nExact) {
  return Math.max(REGISTRATION.absoluteTolerance,
    REGISTRATION.relativeTolerance * nExact);
}
function occupationPass(n, nExact) {
  if (nExact === 0) return n <= REGISTRATION.nullTolerance;
  if (nExact >= REGISTRATION.nonzeroMinExact)
    return Math.abs(n - nExact) <= tolerance(nExact);
  // Not encountered in these registered cases; unregistered small-n rule.
  return null;
}
function integrate(lambda, omega, T, requestedH) {
  const count = Math.round(2 * T / requestedH);
  const h = 2 * T / count, normalizer = 1 / Math.sqrt(2 * omega);
  const ur = Math.cos(omega * T) * normalizer;
  const ui = Math.sin(omega * T) * normalizer;
  let y = [ur, ui, omega * ui, -omega * ur, 0];
  const EInitial = energy(-T, y, lambda, omega);
  // Compute t from the integer index, avoiding accumulated time-coordinate drift.
  for (let j = 0; j < count; j++)
    y = rk4(-T + j * h, y, h, lambda, omega);
  const {alpha, beta} = coefficients(T, y, omega);
  const n = magnitudeSquared(beta), nExact = exactOccupation(lambda, omega);
  const w = magnitudeSquared(alpha) - n, EFinal = energy(T, y, lambda, omega);
  const deltaE = EFinal - EInitial;
  return [lambda, omega, T, requestedH, count, nExact, n,
    Math.abs(n - nExact), nExact === 0 ? null : Math.abs(n - nExact) / nExact,
    ...alpha, ...beta, w, Math.abs(w - 1),
    2 * (y[1] * y[2] - y[0] * y[3]),
    EInitial, EFinal, deltaE, y[4], y[4] - deltaE,
    omega * n, deltaE - omega * n, occupationPass(n, nExact),
    Math.abs(w - 1) <= REGISTRATION.fineWronskianTolerance];
}
function coefficientRoundTrip() {
  // An independent algebraic control, no pulse integration:
  // form u,u' from arbitrary complex coefficients and recover them.
  const omega = 0.73, t = 2.31;
  const alpha = [1.2, -0.31], beta = [-0.23, 0.19];
  const a = rotate(...alpha, -omega * t, 1 / Math.sqrt(2 * omega));
  const b = rotate(...beta, omega * t, 1 / Math.sqrt(2 * omega));
  const y = [a[0] + b[0], a[1] + b[1],
    omega * (a[1] - b[1]), omega * (-a[0] + b[0]), 0];
  const out = coefficients(t, y, omega);
  const errors = [...out.alpha.map((v, j) => Math.abs(v - alpha[j])),
    ...out.beta.map((v, j) => Math.abs(v - beta[j]))];
  return {omega, t, alpha, beta, extracted: out, component_errors: errors,
    pass: Math.max(...errors) <= 1e-14};
}
export function runModeCalibration() {
  const rows = [];
  for (const T of WINDOWS)
    for (const lambda of LAMBDAS)
      for (const omega of OMEGAS)
        for (const h of STEPS) rows.push(integrate(lambda, omega, T, h));
  const fineRows = rows.filter(r => r[3] === STEPS[2]);
  const convergence = [];
  for (const T of WINDOWS)
    for (const lambda of LAMBDAS)
      for (const omega of OMEGAS) {
        const trio = rows.filter(r => r[0] === lambda && r[1] === omega && r[2] === T);
        const coarseDiff = Math.abs(trio[0][6] - trio[1][6]);
        const fineDiff = Math.abs(trio[1][6] - trio[2][6]);
        convergence.push([lambda, omega, T, coarseDiff, fineDiff,
          fineDiff === 0 ? null : coarseDiff / fineDiff,
          ...trio.map(r => r[7]), ...trio.map(r => r[14]),
          ...trio.map(r => Math.abs(r[20]))]);
      }
  const domain = [];
  for (const lambda of LAMBDAS)
    for (const omega of OMEGAS) {
      const pair = fineRows.filter(r => r[0] === lambda && r[1] === omega);
      domain.push([lambda, omega, pair[0][6], pair[1][6],
        Math.abs(pair[0][6] - pair[1][6]), Math.abs(pair[0][6] - pair[1][6])
          / (pair[0][5] === 0 ? 1 : pair[0][5]),
        pair[0][14], pair[1][14]]);
    }
  const failures = rows.filter(r => !r[23] || (r[3] === STEPS[2] && !r[24]));
  const algebra = coefficientRoundTrip();
  return {
    control: "C3 prescribed positive sech-squared pulse; known analytic benchmark",
    generated_at_utc: new Date().toISOString(),
    environment: typeof process === "undefined"
      ? "functions.exec isolated V8 ECMAScript; IEEE-754 binary64; no Python or shell"
      : "Node " + process.version + "; ECMAScript IEEE-754 binary64",
    registration: REGISTRATION, archived_source_reference: SOURCE_REFERENCE,
    implementation: {
      integrator: "classical explicit RK4, unchanged for all runs",
      lambda: LAMBDAS, omega: OMEGAS, requested_h: STEPS, T: WINDOWS,
      initial_state: "e^(-i omega t)/sqrt(2 omega), derivative=-i omega u, t=-T",
      extraction: "alpha=sqrt(omega/2)(u-u'/(i omega))*exp(+i omega T); beta=sqrt(omega/2)(u+u'/(i omega))*exp(-i omega T)",
      exact_occupation: "sin(pi lambda)^2/sinh(pi omega)^2; integer limit exactly zero",
      registered_zero_rule: "n<=1e-8",
      registered_nonzero_rule: "abs(n-n_exact)<=max(2e-6,0.005*n_exact), n_exact>=1e-8",
      work_check_status: "additional exploratory control; no registered work tolerance",
      work_definition: "E=.5*(|u'|^2+Omega^2|u|^2), integrated W'=.5*(Omega^2)'*|u|^2",
      finite_domain: "T=12 to T=16 at each step; no separate domain threshold was registered"
    },
    coefficient_round_trip: algebra,
    table_columns: ROW_COLUMNS, table: rows,
    convergence_columns: ["lambda", "omega", "T", "coarse_medium_n_difference",
      "medium_fine_n_difference", "difference_ratio", "coarse_abs_n_error",
      "medium_abs_n_error", "fine_abs_n_error", "coarse_wronskian_error",
      "medium_wronskian_error", "fine_wronskian_error",
      "coarse_abs_work_residual", "medium_abs_work_residual", "fine_abs_work_residual"],
    convergence,
    domain_columns: ["lambda", "omega", "n_T12_fine", "n_T16_fine",
      "abs_domain_difference", "relative_domain_difference_or_abs_for_zero",
      "wronskian_error_T12_fine", "wronskian_error_T16_fine"],
    domain,
    summary: {
      run_count: rows.length,
      fine_count: fineRows.length,
      all_level_occupation_pass: rows.every(r => r[23]),
      fine_occupation_pass: fineRows.every(r => r[23]),
      fine_wronskian_pass: fineRows.every(r => r[24]),
      coefficient_round_trip_pass: algebra.pass,
      registered_controls_pass: failures.length === 0 && algebra.pass,
      failures,
      max_fine_abs_occupation_error: Math.max(...fineRows.map(r => r[7])),
      max_fine_relative_nonzero_occupation_error:
        Math.max(...fineRows.filter(r => r[5] > 0).map(r => r[8])),
      max_fine_null_occupation: Math.max(...fineRows.filter(r => r[5] === 0).map(r => r[6])),
      max_fine_wronskian_error: Math.max(...fineRows.map(r => r[14])),
      max_fine_abs_work_residual: Math.max(...fineRows.map(r => Math.abs(r[20]))),
      max_fine_abs_deltaE_minus_omega_n: Math.max(...fineRows.map(r => Math.abs(r[22]))),
      max_fine_domain_difference: Math.max(...domain.map(r => r[4]))
    },
    limits: [
      "Numerical binary64 diagnostics; no interval certification.",
      "Known positive sech-squared pulse scattering; no new mathematical result.",
      "Prescribed finite-time vacuum and drive; not an inferred HDBLAST source.",
      "No renormalised stress/current, backreaction, thermal bath, or coupled 5D PDE evolution.",
      "Absolute tolerance can dominate for small occupations; passing is not a physical approximation guarantee.",
      "The driver work integral is per normalized oscillator; it is not a completed cosmological energy ledger."
    ]
  };
}
if (typeof process !== "undefined" && process.stdout)
  process.stdout.write(JSON.stringify(runModeCalibration(), null, 2) + "\n");
