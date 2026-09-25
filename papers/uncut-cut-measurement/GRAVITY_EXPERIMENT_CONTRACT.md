# Physical test contract: source, distance and probe changes

Monty Dabas — v1.0, 25 September 2026.

Status: **PROPOSED_PHYSICAL_CANDIDATE / NO_EMPIRICAL_DATA_INGESTED**.
The version number labels this research packet, not physical validation.

**v1.7 update:** [the atom-interferometer audit](ATOM_INTERFEROMETER_AUDIT.md) now ingests published aggregate phase summaries and implements a separate external numerical reference. No raw-shot dataset or empirical native prediction is available. The historical status above refers to this v1.0 contract; it is not a claim that the package lacks all published observations after v1.7.

## Target and independent calibrations

The first target is the static source-dependent attractive response in U37.
In a regime where the radial approximation is justified, a candidate must
predict the same calibrated instrument reading across several source
strengths and separations, with held-out settings. A different calibration
for every data point is forbidden.

Before evaluating force data, record:

1. The native graph/radial model, weights or capacity-growth rule, cost degree,
   source coupling, outer boundary and finite-source correction.
2. The physical source configuration and independently measured source
   quantity. Mapping source charge Q to mass is a hypothesis to test; a mere
   count of interchangeable tokens is not a mass measurement.
3. The distance/displacement calibration and its relationship to native shell
   labels. Fitting r from the force curve is inadmissible.
4. Probe identity, independently calibrated instrument gain, work/force units,
   sign convention, background subtraction and uncertainty budget. A native
   signed coordinate is not automatically an accelerometer output.
5. Separate training/calibration and held-out configurations. Fit at most the
   declared common amplitude on training data. Retain raw readings and the
   signed source-on minus source-off signal; do not discard its sign when
   testing attraction.

This contract does not name an available instrument or claim one has been
built. No apparatus or raw dataset satisfying these requirements has yet been
attached. Native parameter labels do not provide a physical clock.

## Tests and possible rejection

For a single probe and fixed source geometry in the radial regime, the
quadratic/quadratic-growth candidate predicts the magnitude ratios

    D = F(2r,Q)/F(r,Q) = 1/4,
    S = F(r,4Q)/F(r,Q) = 4.

These statements are conditional on r-independent coupling and the same
background-corrected gain, and on source modification leaving the declared
geometry intact. Finite bodies do not generally obey the point-source ratio
at arbitrary separations. Their response must be computed through the actual
apparatus adapter before applying a rejection gate.

The competing cubic-cost/quartic-capacity model gives D=1/4 but S=2 for a weak
probe. A nonzero onsite cost, a different capacity exponent or a varying
probe coupling gives other falsifiable shapes. The programme has not selected
any one of these alternatives from the primitive; their differences are
candidate-dependent predictions, not a universal prediction of the framework.

For positive response intervals [a-,a+] and [b-,b+] separated from zero, the
observed ratio lies in [b-/a+,b+/a-]. A predicted ratio outside this interval
is rejected under the stated deterministic error bounds. A prediction inside
returns UNRESOLVED. This is a conservative enclosure, not a statistical
confidence level or a proof that the apparatus assumptions are right.
Distance, source-ratio and gain calibration errors must also be propagated
into a predicted-ratio interval. The helper accepts a separately justified
`prediction_error`; it rejects only when observed and predicted intervals are
disjoint. For example, positive ruler intervals R1=[u,v], R2=[w,z] enclose the
ideal inverse-square ratio in [(u/z)^2,(v/w)^2]. Include the source/gain ratio
intervals by positive interval multiplication. Setting prediction_error=0
assumes those calibrations are exact; it is used only for the ideal fixtures.
Correlated errors do not invalidate enclosure if the component bounds are
valid, but can make it unnecessarily wide. Near-zero or sign-ambiguous
readings require a different signed analysis; the executable ratio helper
rejects them as invalid inputs.

The current fixtures include correctly retained ratios under bounded centre
perturbations and rejected 1/r and nonlinear-source controls. **They are
synthetic.** They certify the stated calculations, not gravity in nature.

## Equivalence principle is a separate gate

If an independently established physical inertia adapter converts force to
acceleration, the native prediction contains a factor q_p/m_p for probe p.
Universality requires that ratio to be the same across probe compositions.
Native agreement NP0 does not imply this physical identification. Fitting a
new gain per material to enforce equal accelerations would erase the test.

The MICROSCOPE final result compared titanium and platinum alloy test masses
and reported no detected violation, with
`eta = (-1.5 +/- 2.3 statistical +/- 1.5 systematic) x 10^-15`.
This published result is a future external benchmark, not evidence that the
native candidate has passed. The error components are not repackaged as a
hard deterministic interval. Source: Touboul et al., *Physical Review Letters*
129, 121102 (2022), [arXiv:2209.15487v1](https://arxiv.org/abs/2209.15487v1).

## Existing inverse-square experiment as a benchmark

Lee et al. used a torsion balance and patterned rotating attractor at
separations from 52 micrometres to 3.0 millimetres, and reported agreement
with the Newtonian model. Their quoted gravitational-strength Yukawa range
limit is specific to that deformation and apparatus. It cannot be copied
onto an arbitrary native shell-capacity correction. Source: *Physical Review
Letters* 124, 101101 (2020),
[arXiv:2002.11761v1](https://arxiv.org/abs/2002.11761v1).

Only the authors' published abstracts and bibliographic records were inspected
for these benchmark statements on 25 September 2026. Their raw measurements,
likelihoods and apparatus transfer functions have not been ingested or refit.
No experimental PASS, new numerical exclusion or priority claim is made.
Experimental comparisons are external observations, not imported foundational
field equations.

## Next evidence gate

Choose a physically specified source and probe, justify their mapping into
the native variables and predeclare the ruler/capacity law. Attach a raw
dataset and apparatus model, or acquire such data. Fit shared parameters on
calibration settings, then evaluate the held-out source and distance changes.
Independently deriving the capacity-growth law from a selected native family
would remove a central constitutive assumption. Neither task is completed by
adding further synthetic points to the present certificate.
