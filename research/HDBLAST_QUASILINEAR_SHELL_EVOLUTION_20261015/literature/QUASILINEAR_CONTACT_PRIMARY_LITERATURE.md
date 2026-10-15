# Primary literature for the quasilinear shell-contact problem

Retrieved 15 October 2026 UTC. This is a bounded continuation of the
previous angular-uniform literature review. It reports primary passages
actually inspected and separates them from a proposed application to the
registered five-dimensional Einstein–scalar shell model.

## Finding

An accessible full primary paper supplies substantially stronger support
than the earlier abstract-only Shibata–Nakamura citation:
Shibata–Kikuchi, *Journal of Differential Equations* 80 (1989), 154–197,
DOI [10.1016/0022-0396(89)90099-5](https://doi.org/10.1016/0022-0396(89)90099-5).
Its Theorem 1.1 proves a finite-Sobolev quasilinear vector Neumann theorem,
and Section 6 gives the relevant variable-coefficient linear theorem.
The full paper is available through the publisher's open archive.

This does not directly prove the HDBLAST mixed-field theorem. In this
literature, “mixed problem” commonly means an initial–boundary value
problem; it does not mean that some field components obey Dirichlet rows
and others obey Neumann rows. The model's constant complementary
Dirichlet/Neumann projectors, lower-order coupling, exact corner jets,
parameter regularity, and continuation to the chosen finite contact
interval still require a matched proof.

## Shibata–Kikuchi 1989: exact theorem and restrictions

[Publisher landing page](https://www.sciencedirect.com/science/article/pii/0022039689900995)
identifies the article as an open-archive contribution by Yoshihiro
Shibata and Masahiro Kikuchi. The complete 44-page primary PDF was
retrieved; Sections 1, 3–4, 6–7 were read at the theorem, compatibility,
iteration and convergence passages. PDF extraction contains occasional
OCR errors, including confusion between a spatial index and the letter
t; the following is a structural summary rather than a verbatim
transcription of an affected formula.

The displayed problem (N) is a vector equation with normalized leading
time derivative and spatial divergence fluxes depending on
`(t,x,u,du)`. Every boundary row is the associated normal flux plus a
nonlinear lower-order or admissible derivative term.

Theorem 1.1 assumes spatial dimension n at least 2 and hypotheses
(A.1)–(A.5). Their substance is:

1. Smooth coefficient functions on a fixed field/derivative patch, with
   bounded derivatives; a zero-at-zero condition is used on unbounded
   domains for integrability.
2. Spatial flux-Jacobian transpose symmetry and symmetric time-drift
   matrices; the boundary derivative matrices have their specified
   symmetric/skew-symmetric structure.
3. A coercive spatial bulk-plus-boundary quadratic form, up to an L2
   shift. The authors explicitly distinguish this from strong ellipticity
   alone.
4. The extra boundary derivative operator has no normal derivative.
5. The normal-time boundary quadratic form is nonnegative. The
   non-strict inequality permits a reflecting equality case.

For integer K at least `floor(n/2)+3`, Theorem 1.1 gives a unique short-time
solution in

`X^{K,0} = intersection_{j=0}^K C^j([0,T]; H^{K-j})`.

The initial fields are H^K and H^{K-1}. Bulk and boundary sources have
the displayed mixed space/time regularity in (1.2)–(1.6), including
Lipschitz control of the highest indicated time derivatives into L2 and
H^{1/2} respectively. The boundary compatibility conditions hold through
order K−2, as defined in Section 3, (3.4). For n=4 the theorem's numerical
threshold is K at least 5. This threshold belongs to that theorem's full
hypotheses, not to an independently verified HDBLAST existence theorem.

Section 4 fixes compatible initial jets and an extension u^0 in (4.2),
reduces the nonlinear problem to coupled hyperbolic/elliptic problems
for u and its time derivative, and defines successive approximations
in (H)_p and (E)_p. Section 7 proves the contraction (4.11), convergence,
the identity `v = u_t`, and continuity in the top spatial H^K norm.
The introduction expressly warns that naively iterating a fully
nonlinear Neumann boundary operator can lose derivatives. The theorem
is consequently not permission to omit a boundary/iteration estimate.

For the frozen problem, Theorem 6.1 assumes (A.6.1)–(A.6.5): a smooth
bounded coefficient part plus a finite-Sobolev part
`Y^{K-2,1}`, the specified transpose symmetries, coercivity, no extra
normal boundary derivative, and nonnegative normal-time boundary form.
It gives `X^{K-1,0}` existence with compatible initial data and the
stated source norms; (6.10) is the basic energy inequality and (6.11)
starts the higher-order estimate. The inhomogeneous boundary source in
(6.10) is charged in H^{1/2}; it is not an arbitrary bare-L2 Neumann
control theorem.

For the straightened HDBLAST system a common scalar metric principal
part is symmetric in the fiber indices, and near-background positivity
can supply the spatial coercivity. A frozen iterate satisfying the
homogeneous Dirichlet rows has zero normal–tangential metric block at
the shell. This is precisely the structural check needed to obtain
zero normal time drift and make the Neumann rows conormal. A scalar
all-Neumann frozen theorem is therefore a useful building block.

However, Theorem 1.1 is not directly applicable to the whole nonlinear
system merely by calling it “mixed”: its displayed problem has no
Dirichlet rows. In addition, its nonlinear structural assumptions are
imposed on an open state patch, whereas the physical normal–tangential
block condition here is guaranteed on the submanifold on which the
Dirichlet variables vanish. Frozen iterates must actually preserve
those rows. A proof combining the scalar Dirichlet and Neumann frozen
solves and controlling their coupled lower-order terms remains part
of the present argument.

## Nakamura–Shibata 1986: accessible predecessor

Gen Nakamura and Yoshihiro Shibata,
“On a Local Existence Theorem for Quasilinear Hyperbolic Mixed Problems
with Neumann Type Boundary Conditions,” *Proceedings of the Japan
Academy, Series A* 62 (1986), 117–120,
DOI [10.3792/pjaa.62.117](https://doi.org/10.3792/pjaa.62.117).
The complete four-page primary PDF was retrieved and read.

Its Section 2 theorem treats a vector spatial-divergence system with
natural Neumann/nonlinear Robin rows, symmetric coercive flux Jacobian,
and compatible data. It uses the more conservative threshold
`L >= floor(n/2)+8`, gives `E_T^L = intersection C^j H^{L-j}`, and requires
corners through L−2. The introduction already describes the
hyperbolic/elliptic device for avoiding boundary derivative loss.
This is a primary theorem, rather than evidence inferred from a later
citation, but it is again an all-Neumann problem. The stronger 1989
Shibata–Kikuchi theorem explicitly explains the improvement in its
Sobolev threshold.

The 64-page 1989 Shibata–Nakamura article,
DOI [10.1007/BF01180683](https://doi.org/10.1007/BF01180683), was searched
again. No accessible complete body was located in these searches.
No claim about its exact full-body hypotheses is added here.

## Dependence on data and a prescribed reference lifespan

H. Beirão da Veiga, “Structural Stability and Data Dependence for Fully
Nonlinear Hyperbolic Mixed Problems,”
[author-hosted primary PDF](https://people.dm.unipi.it/beiraodaveiga/pdf/hbv-59.pdf).
The complete ten-page primary PDF was retrieved. The theorem, the
literal-applicability qualifications in Remark 1.3, and the finite-time
patching argument in Section 2 were inspected.

For the scalar natural-Neumann class displayed in (1.2), Theorem 1.1
starts with a reference solution on a fixed interval [0,T0]. Under the
specified converging coefficients, compatible H^{k+1}/H^k initial data,
and mixed source norms, sufficiently close problems have solutions on
the whole reference interval and converge strongly in the top solution
norm. Here `k > n/2+1`. The proof explicitly patches a common short-time
result along the reference solution. Remark 1.3 retains additional
source/coefficient assumptions for literal use of Shibata–Kikuchi's
existence theorem; those qualifications must not be dropped.

The paper comments that its dependence argument can be adapted to
vector fields and to a material derivative with boundary-tangential
drift. It also explicitly says that some generality statements concern
the dependence argument and do not themselves extend the existence
theorem. Consequently this supplies a precise precedent for the
background-interval continuation strategy, not an automatic mixed
Einstein–scalar theorem or a twice-differentiable solution map.

## Other screened full primary material

Lu–Xin–Hu, “Local existence of solutions for the Neumann problem of the
nonlinear elastodynamic system outside a domain” (2011),
[publisher primary page](https://www.sciencedirect.com/science/article/pii/S0022247X1100103X).
The complete 17-page PDF was retrieved. Theorem 3.1 and its iteration
(3.9)–(3.10) were inspected. The displayed nonlinear theorem assumes
smooth initial data and all-order compatibility, with lifespan/control
based on the specified finite Sobolev norms (r at least 6 in spatial
dimension 3). It is not cited here as a finite-Sobolev nonsmooth-data
existence theorem. It provides an explicit example of preserving the
same initial compatibility jets at every frozen step.

Two additional retrieved full primary texts were screened at their
introductions and relevant scope statements: Shibata's 1993 EVEQ92
lecture on nonlinear dissipative evolution equations and Ren–Xin's
EJDE article on three-dimensional exterior Neumann quasilinear
waves. Their dissipative/global or dimension-specific hypotheses do
not replace the local mixed-field proof sought here. No theorem from
either is imported into the contact argument.

## Coverage and claim boundary

This subtask added six ordinary web searches and two paper-index
searches to the inherited review. Seven URL retrievals produced six
complete primary PDF extractions and one publisher landing page.
Only the sections described above were closely read. The paper-index
results, including a 2026 two-dimensional exterior-lifespan candidate,
were discovery material; no theorem is inferred from their abstracts.
The search is targeted and is not an exhaustive review of the web.

Raw tool outputs and complete extracted papers remain internal reference
material. This synthesis and its provenance ledger can accompany a
research checkpoint. Literature support does not certify the new PDE
proof, establish external novelty, or constitute external peer review.
