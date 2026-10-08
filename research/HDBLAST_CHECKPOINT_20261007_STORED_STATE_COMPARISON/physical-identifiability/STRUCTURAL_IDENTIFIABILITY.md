# Effective forcing and the limits of dimensional identification

Scientific study date: 7 October 2026, America/Los_Angeles; primary reading and
review occurred on 8 October 2026 UTC. This elementary structural argument
describes the present experiment and a possible physical continuation. It is not
a higher-dimensional origin model, an external novelty assessment, or an
observational detection.

The present HDBLAST calculation fixes a four-dimensional de Sitter background,
H=1 in its declared units, a minimally coupled scalar with mass squared two,
and the prescribed linear response

    U' = W,       W' = 2ikW - g,
    g = T h := 4L²h - 2Lh' - h'',       L = -1/eta.

Its two source profiles are h₁(eta)=B(eta+4) and
h₂(eta)=(eta+4)B(eta+4). The incoming prescription and subsequent direct
operators are fixed. Certifying that calculation establishes statements about
this response problem. It does not derive a five-dimensional Einstein action,
a bulk solution or junction condition, the origin of the source's energy,
a hot radiation era, or the Big Bang.

## Proposition 1: containment of effective inputs prevents origin identification

Let β denote the experimental specification held common between two model
classes: four-dimensional background and evolution operators, initial quantum
state, complete contact/subtraction data or their fixed dependence on the
common effective input, time and momentum domain, observable and
detector/selection map, and all other nuisance/noise laws. Let G be a class
of admissible effective forcings for which the forward problem is well posed.
Write F_β(g) for the predicted observable, or more generally for the COMPLETE
probability law of the observed data. Taking F_β to be a probability law avoids
confusing equality of a mean with equality of fluctuations or correlations.

A proposed D-dimensional family induces g_D: Θ_D -> G. A four-dimensional
comparison family induces g₄: Θ₄ -> G. Assume:

1. Both descriptions use the same β, including the initial state and complete
   contact terms. All their influence on this experiment is represented by g.
2. For every θ in Θ_D there is at least one ν in Θ₄ such that
   g₄(ν)=g_D(θ) on the whole domain relevant to F_β.
3. There is no additional measured channel depending on the bulk description
   except through this same effective input and the common β. If contacts or
   metric perturbations depend on additional data not determined by g, those
   data must be matched explicitly or included in an enlarged effective input.

Then every D-dimensional observable law is also a four-dimensional observable
law:

    { F_β(g_D(θ)) : θ in Θ_D }
       is a subset of { F_β(g₄(ν)) : ν in Θ₄ }.

In particular, the origin label D versus 4 is not identifiable at any paired
predictions covered by these assumptions.

Proof. For a fixed θ choose ν from assumption 2. Substitution of identical
inputs into the same well-posed forward and measurement map gives
F_β(g_D(θ))=F_β(g₄(ν)). This holds for every θ, proving the inclusion. If these
outputs are probability laws, every measurable statistic or decision rule has
the same distribution under the paired hypotheses. A test that has false
positive probability at most α for all four-dimensional alternatives therefore
has power at most α against each exactly replicated D-dimensional alternative.
No precision or sample-size improvement separates identical laws. ∎

For the displayed scalar ODE, this substitution can also be proved directly.
Equal g and equal incoming data give a difference satisfying δU'=δW,
δW'=2ikδW and zero initial values, hence δU=δW=0 for every time and momentum.
Identical direct operators and contacts then give identical responses. The
argument does not depend on having only nine probes or a finite saved grid:
equality can hold at every node, every time and every momentum in the domain.

This is an implication with stated hypotheses, not a claim that every bulk
model reduces to this scalar forcing or that every physically restricted 4D
theory can realize every g. A new channel, distinct preparation state or
background, or a restricted rival that excludes g_D changes the hypotheses
and must be analyzed explicitly. An unrestricted effective forcing class need
not possess a local, unitary ultraviolet completion for every element; it is
an observational nuisance class, not a claimed microscopic theory.

Any Ward identity satisfied by F_β(g) in the shared experiment also holds for
the replicated four-dimensional prediction. A small certified stored-state
error improves confidence in evaluating F_β; it does not determine which
physical mechanism generated g. Likewise, normalization is a state constraint,
not an origin label.

The likelihood conclusion concerns paired laws. If densities exist with a
common reference measure, their likelihoods agree pointwise for paired
parameters, and the maximum likelihood over the encompassing four-dimensional
family is at least that over the D-dimensional subfamily. Bayesian evidences
need not agree: priors, parameter volumes and restrictions can favor one
specified model. Equal pushforward priors on effective observable laws would
give equal prior-predictive distributions; different priors do not manufacture
new observational directions. Support relative to declared rivals is meaningful
when reported with those restrictions and sensitivity analyses.

## Corollary 2: the fixed two-profile linear experiment has rank at most two

If the present profiles are used as a two-amplitude model with real amplitudes
a = (a₁,a₂) in R², write

    h = a₁h₁ + a₂h₂,      g = a₁g₁ + a₂g₂,      g_j = T h_j.

With background, preparation prescription and operators fixed, the unique
first-order response is linear in (a₁,a₂). Thus any collection of linear
observables of that response, with a fixed baseline removed, has the form

    y - y₀ = X a,        X = [R(g₁), R(g₂)],        rank X <= 2.

Complex mode outputs may be represented by their real and imaginary parts;
adding times, nodes, stress components or linear integrals adds rows to X,
not new source directions. For an infinite collection the response still lies
in the span of two fixed response functions. The rank can be less than two if
the chosen measurement cannot distinguish the two responses. This rank bound
requires linear observables at the stated perturbative order; nonlinear
observables can have a larger linear span while still depending on only the
same two amplitudes.

Suppose a bulk proposal supplies only amplitudes a=A(θ), where the parameter
domain Θ is open in Rᵖ and A: Θ -> R² is continuously differentiable (C¹).
Its Jacobian is

    D_θ y = X D_θ A,       rank(D_θ y) <= 2.

Any θ and θ' with A(θ)=A(θ') are exactly indistinguishable in this experiment.
Near a regular point where the observable map has constant rank r, the
constant-rank theorem gives local level sets of dimension p-r. Therefore a
regular p>2 family cannot determine all of its continuous parameters from
these two response directions. This regular-point statement does not claim
the same fiber dimension at singular points, parameter boundaries, or within
discrete or externally restricted parameter sets. Even an injective restricted
two-parameter amplitude map identifies parameters only within that model; it
does not defeat Proposition 1's encompassing 4D rival.

The currently computed two unit-profile responses are controlled experiments,
not two cosmological amplitudes inferred from observations. No fit of a₁,a₂ or
bulk parameters has been performed here. Additional independently correlated
response kernels or channels can enlarge the information available. An action
restriction or external parameter constraint can reduce the model's dimension.
A prior can select among equivalent parameters but its role must be stated;
it does not raise the data's structural rank.

## Memory, boundaries and quantum noise require a larger effective description

A moving boundary, a bulk continuum or an environment can induce memory,
dissipation, nonlocal boundary operators and noise. Matching a single mean
forcing g then does NOT establish equality of observable probability laws.
One must instead specify an effective input such as

    η_eff = (g, retarded kernels, noise and higher cumulants,
             preparation correlations, boundary conditions, contacts,
             other coupled-channel operators),

and use a forward map F_β(η_eff). Proposition 1 still applies if the allowed
4D effective comparison family contains every such induced input, but that
containment is now a separate and stronger premise. It must not be inferred
from the scalar ODE or from a mean response alone.

For a specified open-system action and bath state, response and noise are
jointly constrained; causality, positivity and, when appropriate, equilibrium
KMS/fluctuation-dissipation relations restrict their independent variation.
Gaussian kernels alone do not fix non-Gaussian higher correlations. Conversely,
a spectral shape or a fluctuation-dissipation relation is not by itself a label
of spatial dimension: an ordinary four-dimensional environment can reproduce
an allowed bath spectrum and preparation state over the measured domain.
This is conditional on realizing the required joint kernels and couplings;
no claim is made that an arbitrary kernel has a local, stable 4D completion.
Li's explicit open-EFT construction uses short inflaton modes as the bath and
matches nonlocal, history-dependent kernels within that field-theory setup
[2, Secs.3.2.1 and 5]. It illustrates why memory and noise need a joint model,
without establishing any HDBLAST bulk bath.

The same distinction holds for brane equations. A conserved projected stress
or Weyl term need not have an exhibited regular inducing bulk. The recent
brane-collapse analysis [3, Secs.2.1,4.4,5.3] explicitly leaves bulk realization
and flux orientation assumptions separate. The earlier origin proposal [4,
Sec.4] introduces a thermal bulk atmosphere and equilibrium conditions to
obtain a curvature spectrum. These are model premises that a prescribed pulse
and a nonzero excitation coefficient do not supply.

## A concrete restricted joint prediction, used only as an example

Rao et al. [1] specify a particular parity-odd teleparallel five-dimensional
parent action, a compact spacelike circle and a low-energy KK zero-mode sector.
In their declared flat-FRW electromagnetic-vacuum setting, field variables and
physical branch phi>0, Eqs.(44),(48),(49) give

    Delta omega_EM,A² = 6 Delta omega_GW,A².

This is a coefficient relation in their leading-WKB propagation description,
fixed by that action. It is not a factor-six relation between arbitrary CMB
and gravitational-wave measurements. Their Sec.VI further obtains, at leading
order and ONLY for common emission and detection endpoints,

    Delta alpha = 3 Delta Phi_GW,WKB.

The photon rotation is one half of the photon helicity phase difference. The
fixed-k GW eikonal phase differs from a CBC Fourier-domain waveform phase.
CMB photons and ordinary low-redshift CBC events have different propagation
endpoints; the source also requires inclusion of radion-induced friction in a
quantitative waveform reanalysis. These restrictions are explicitly present
in the inspected paper, Eqs.(52),(55),(56) and the following discussion.

A relation of this kind provides a falsifiable restriction of a specified
model against specified rivals. A 4D action with the same effective coupling
relation can reproduce it; passing it would not uniquely detect a fifth
dimension against every 4D theory. Violating it, after its approximations and
measurement uncertainties are controlled, can reject that particular action
or its premises. The example supplies NO new HDBLAST prediction, amplitude,
constraint or likelihood. Adding those parity operators to HDBLAST would be
a different model requiring a new derivation and test contract.

[5] gives a complementary observational example: a phenomenological GW
amplitude-propagation law is tested with population, selection and catalog
modeling, and its inferred dimension depends on crossover-scale priors. Such
a law is not derived by the present scalar pulse experiment. Its numerical
posterior or dataset is not imported here.

The accompanying `MODEL_TEST_CONTRACT.json` defines the information needed
before a model-specific explanatory comparison. A successful test can support
a restricted model relative to its tested rivals. It is distinct from a
model-independent detection of dimensional origin or an established account
of Big Bang causation. The elementary structural argument has no external
novelty claim, and makes no prediction about prizes or scientific recognition.

## Primary references and reading scope

The statements below are based on verified cached primary-text bytes and
selected passages, not an independent replay of their full calculations.
`PRIMARY_CITATION_PINS.json` records versions, source URLs, PDF/text hashes and
passage locations. The primary PDF/text cache is not part of this public text.

1. H. Rao, Y. Zheng, Q.-Z. Hou, J.-W. Ou and C. Zhu, *Correlated parity violation
   in gravity and electromagnetism from five-dimensional spacetime*,
   [arXiv:2608.09299v3](https://arxiv.org/abs/2608.09299v3), Sec.II.A,
   Sec.V Eqs.(39),(41),(44)-(49), Sec.VI Eqs.(52),(55),(56).
2. Y.-Z. Li, *Stochastic inflation as an open quantum system II: open effective
   field theory and stochastic matching*,
   [arXiv:2605.21929v3](https://arxiv.org/abs/2605.21929v3),
   Secs.3.2.1 and 5; selected passages concerning bath preparation, scale
   separation and nonlocal/non-Markovian matching.
3. R. Sengupta and C. Singha, *Bouncing Dust Collapse and Black-to-White Hole
   Transition on the Brane*, [arXiv:2610.04353v1](https://arxiv.org/abs/2610.04353v1),
   Secs.2.1-2.2,4.4,5.3. Its timelike extra direction and negative tension are
   specific assumptions, not premises of HDBLAST.
4. R. Pourhasan, N. Afshordi and R. B. Mann, *Out of the White Hole: A Holographic
   Origin for the Big Bang*, [arXiv:1309.1487v2](https://arxiv.org/abs/1309.1487v2),
   Eq.(1.1), Sec.4 Eqs.(4.4)-(4.15);
   [doi:10.1088/1475-7516/2014/04/005](https://doi.org/10.1088/1475-7516/2014/04/005).
5. A. Chen and J. Zhang, *Searching for Extra Dimensions with Gravitational
   Waves: Dark-Siren Constraints from GWTC-4*,
   [arXiv:2606.14549v1](https://arxiv.org/abs/2606.14549v1),
   Sec.II Eq.(2.16), Secs.III-VI. No independent data analysis is claimed.

Higher-dimensional cause: NOT_ESTABLISHED. External mathematical novelty:
NOT_ASSESSED. Existing numerical failure and unresolved-gate statuses retain
their separate meaning and are not altered by this structural proposition.
