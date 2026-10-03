# Independent scientific scope review

Reviewed 2 October 2026, after the registered scalar-response run. This review
read the stored primary and independent results, both validator reports, the
two producer implementations, the frozen analytic protocol, and the analytic
spectral companion. `REVIEW_INPUTS.json` records exact paths and hashes. No
frozen input was edited and no new source, mode, response, stress, or spectral
energy calculation was executed numerically.

## Assessment

The reported PASS is supported as a calibration of the homogeneous linear
variance susceptibility for a minimally coupled scalar on the fixed de Sitter
reference x0=r=2H^2, with an unchanged incoming BD state and the inherited
finite prescription. I found no mathematical or interpretive blocker to that
bounded claim in the reviewed evidence.

The independent implementation evolves `u'=w`, `w'=2ik w-s`, where
`delta v=v u`, from zero response data. Its exact homogeneous step propagator
and local eight-node source quadrature are distinct from the primary global
memory quadratures. It preserves complex response amplitudes, checks rather
than imposes the linear Wronskian, combines mode and subtraction terms before
the momentum sum, and archives the resulting mode data. The primary route
implements the matched finite-cutoff formula and two continuum memory forms
with the correct sign, 2k frequency, normalization, scale factors, and finite
contact. This is a useful independent implementation comparison, while both
routes necessarily inherit the same physical model and subtraction.

The compared arrays contain the intended two pulses, six observations, and
three cutoffs. The 36 comparison records and 17 mutation outcomes agree
between the normal and optimized validators. The exact algebra suite is a
separate body of 46 identities and 14 mutations. The derivative tail module
has an analytic total-variation proof and directed interval evaluation; it
does not enclose every integration and floating-point error in the run.
Reporting the near-machine-level finite-K differences as a certified total
accuracy would overstate the evidence.

Strictly pre-pulse zero, post-pulse memory, excited-state diagnostics, the
initial Bogoliubov boundary diagnostic, and the fixed-reference stationary
contact check test distinct failure modes. They are appropriate controls for
this experiment. They neither establish a unique cosmological initial state
nor supply the missing stress and metric responses.

## Contacts and work

The finite `gamma_E+1` variance contact, the quadratic mass-law current
contact, and the spectral source-work drift are three separate consequences
of the chosen action and source variables. None can be dropped on the grounds
that a canonical plane-wave mode is being used. In particular,

    E_can^(2)=(1/2)integral s' q_M
                -(1/(32pi^2))integral(a'/a)s^2

has the checked positive drift coefficient and the stated Fourier convention.
The completed run contains no numerical test of this identity or evaluated
energy yield. Its canonical energy is also not the full physical minimally
coupled energy density. The reader summary in
`SPECTRAL_CONTACT_READER_NOTE.md` preserves these distinctions.

## Next calculation and boundaries

The inherited `COMMON_ACTION_SMOOTH_FRW.md` supplies an explicit covariant
finite prescription and fourth-order stress subtraction for time-dependent
mass. The sphere-restricted stationary action by itself would not suffice.
Using that fuller inherited prescription, the accompanying post-run analytic
proposal derives stress contacts and specifies finite-K and continuum Ward
comparators. This is new analytic follow-up material, not an additional test
that was preregistered or executed in the completed scalar experiment.

The next numerical task should independently integrate the matched physical
rho and p responses, compare them with the analytic closures, and check both
Ward and trace constraints with their actual finite-K local terms. Derivative
and stress-tail errors need new bounds and a new public registration. A
closure defined through Ward cannot alone serve as its independent test.

No output here can be inserted as the exact quantum Hessian of the shifted
gamma=0.01 stationary endpoint. Its actual propagator, source translation,
metric response, and coupled bulk/boundary conditions remain necessary.
Passing a dimensional scale screen and this scalar kernel calibration does
not establish coupled initial data, stability, relaxation, physical radiation,
heating, or thermalization. These are limits of the calculated observables,
not indications that the completed calibration failed.
