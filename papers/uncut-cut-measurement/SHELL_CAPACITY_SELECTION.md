# Can the native law select quadratic shell capacity?

Monty Dabas — focused continuation v1.1, 25 September 2026.

**Verdict under the current primitive:** No. The admissible native carrier and
generator rules of U27 allow different shell-capacity profiles with the same
carrier count, local valence and uniform link stiffness. Consequently those
rules alone cannot prove A(r) proportional to r squared. This manuscript
derives quadratic *leading* capacity from one explicitly supplied
three-direction graph and shows why neither that count nor an information
bound proves the physical pointwise response. All claims concern A; no other
physical law is introduced here.

The [four immutable repository sources](SHELL_GROWTH_SOURCE_PINS.json) bind
the admissible graph family, its existing selection limits and the v1.0
capacity definition. The combinatorial graph counts below are proved directly.

## 1. The exact object whose growth needs explaining

In U36–U37, a declared connected graph Gamma has positive link stiffness
c_uv. Choose a root o. Its native graph-distance labels d(v) are integers;
they are not physical metre readings. For the ball B_n={v:d(v)<=n}, define

    C_n = sum_{u in B_n, v outside B_n} c_uv,
    A_n = ell_n C_n,                                             (39.1)

where ell_n is an independently specified displacement calibration. The
source-cost cut ledger gives outward *total current* Q across every B_n
containing the source and no other sources. Thus current per unit aggregate
capacity is Q/C_n. Only if the observed field is shell-constant, or a lawful
averaged observable obeys the requisite identical gradient, does the local
finite-difference response become q Q/A_n. Equation (39.1) is a capacity
definition from graph/cost/ruler data, not an inverse-square definition.

In the unweighted graphs below all c_uv=1. The same U27 two-component
carrier, individually addressed R and edge-mixing operators, one grading,
and positive source cost are lawful on all of them. The graph is a declared
coupling graph; no embedding in an existing spacetime is a premise.

## U39. Quadratic capacity cannot be inferred from current native axioms

There are two connected graphs on **the same 49 carrier labels**, both
four-regular with exactly 98 links of stiffness one:

* Gamma_T is the 7-by-7 periodic square graph, links in its two cyclic
  directions. For root (0,0), its cut capacities C_0,C_1,C_2 are 4,12,20.
* Gamma_C has vertices modulo 49 and links x to x+/-1 and x+/-2. For root
  zero, its capacities are 4,6,6.

Both have every native step that U27 licenses on a supplied simple graph.
Both are connected, so the U28 whole-component invariant-form line is the
same type, and an exact nonzero readout needs the entire 98-dimensional
real carrier for full addressed-word closure. Neither invariant-form nor
whole-carrier memory rank chooses Gamma_T over Gamma_C. A native topology
selector or extra growth law is needed before (39.1) can have a unique
asymptotic exponent, even when cardinality, degree and weights are held fixed.

**Proof.** The graph constructors enumerate the distinct unordered links.
Each vertex in Gamma_T has two horizontal and two vertical neighbours;
each in Gamma_C has the four neighbours at offsets +/-1,+/-2. Both have
49*4/2=98 distinct unit links. Counting links leaving the root and then its
radius-one and radius-two balls gives the stated profiles. U27's mixing
matrix and its real anti-self-dagger/grading identities are defined for every
one of these links, independently of which graph was supplied. Since the
profiles differ while all stated native rules hold, those rules do not entail
a unique capacity growth law. This is a countermodel, not a physical claim
that either graph represents the world.

The obstruction persists after imposing equal interior degree alone:
on a four-regular rooted tree truncated and grounded at depth three,
`C_n=4*3^n` for n=0,1,2, whereas the square graph gives `4,12,20`.
The infinite tree has exponentially growing spheres. A path gives constant
cut capacity. Graph homogeneity or unit link weight by itself is therefore
not a selection principle for exponent two. Finite truncation is used only
to check the indicated shells and does not establish an infinite carrier
in the U27 finite model.

## U40. What a declared three-direction graph actually derives

Supply vertices x in the integer grid with d native independent commuting
direction labels and nearest-neighbour links x to x+/-e_i, all of unit weight.
The graph-distance ball B_n is specified by sum_i |x_i|<=n. This is a
**declared graph-sector hypothesis**, not a theorem selecting three physical
directions from K, R or a cut. Let S_n^d be its sphere.

For n>=1, choose the k nonzero coordinates, their signs and their positive
parts summing to n. Direct counting gives

    |S_n^d| = sum_{k=1}^{min(d,n)} 2^k binom(d,k) binom(n-1,k-1).

A vertex with z zero coordinates has d+z outward links: one for each nonzero
coordinate and two for each zero coordinate. Counting the zeros by first
choosing their coordinate proves

    C_d(n) = d |S_n^d| + d |S_n^(d-1)|,   n>=1;
    C_d(0) = 2d.                                               (40.1)

The leading term is `[d*2^d/(d-1)!] n^(d-1)`. For d=3,

    |S_n^3| = 4n^2+2,
    C_3(n) = 12n^2+12n+6,   n>=0.                          (40.2)

For d=2, `C_2(n)=8n+4` (and C_2(0)=4); for d=1, `C_1(n)=2`.
The finite exact certificate checks (40.1) against independent edge counts
on balls for d=1,...,4 and radii through four. With uniform link stiffness
c and constant independently calibrated spacing ell, `A_n=ell*c*C_3(n)`
has quadratic *asymptotic* growth when r_n=n ell:

    A_n / r_n^2 -> 12 c/ell,     as n -> infinity.                 (40.3)

The finite formula has nonzero linear and constant corrections; it is not
`A_n=K r_n^2` at every n for that ruler. Choosing r_n after seeing C_n to
eliminate those corrections would repeat W50's fitted-ruler error.

**Selection boundary.** This proves a count in the supplied rank-three
nearest-neighbour graph. It does not derive why three independent direction
labels, this graph, uniform weights, a physical metre, or the required
radial response follow from the uncut primitive. K, R and KR being three
displayed internal matrices does not turn them into three translations:
their composition and grading obey their own algebra, and no theorem in
U1–U38 identifies them with graph directions.

## U41. Counting channels does not prove a local inverse-square field

The cut ledger on the finite three-direction ball still has total current Q
through every shell. Its *average oriented edge current* is Q/C_n. This
statement follows by division and does not require a uniform potential.

The nearest-neighbour cubic graph is not shell-equitable. At distance two,
an axial vertex such as (2,0,0) has one inward neighbour, while a vertex
such as (1,1,0) has two. With a grounded radius-three shell and a unit
source at the origin, solving the exact positive source-cost equations gives
two different values of the field on shell two. Its layer-mean potential
drops are `1/6, 13/378, 5/378`; aggregate inverse capacities are
`1/6, 1/30, 1/78`. The second and third entries differ.

**Proof of the local distinction.** The grounded finite graph has a unique
source solution by U36. Direct rational elimination gives the displayed
means; the certificate checks the source equations and shell currents from
the edge list independently. Shell two contains axial and nonaxial vertices
with different field values. Were the aggregate inverse capacity equal to
the mean potential drop, the second terms would be equal; `13/378 != 1/30`.
Thus deriving C_3(n) is insufficient to deduce a local pointwise inverse-square
force under the current finite assumptions. The finite difference may have
an asymptotic limit, but it is not established here.

An edge-equitable shell condition *is* sufficient for radial stationarity:
each vertex in one layer must see the same total inward and outward link
weights, with shell-constant sources and boundary. Restricting the strictly
convex cost to radial fields then gives a stationary point of the full cost;
uniqueness makes the field radial. The four-regular tree satisfies this
condition and gives exact drops `1/4,1/12,1/36`; its capacity growth is
exponential. Radial response and quadratic growth are separate requirements.

## 5. Three proposed routes, audited only for A(r)

| Route | What it gives under explicit premises | Remaining selection |
|---|---|---|
| Dimensionality | Rank-three integer-direction graph yields (40.2), hence quadratic leading capacity | Select three commuting directions, graph, weights, ruler; additionally justify local radial response |
| Holographic bound | An independently stated upper bound `A_n <= K r_n^2` only limits growth | Constant-capacity path satisfies the same bound; equality or a matching lower bound has to be proved |
| Information capacity | If a native cut carries exactly its independent shell link channels at common cost, its counted channel capacity is (39.1) | U28's full observability rank is identical for the two 49-carrier countermodels; it does not fix their shell edge counts |

Calling a shell law "holographic" or "information capacity" does not change
its proof obligation. None of these routes alone supplies the unique
primitive-to-graph selector requested. A future selector must reject at
least one of the two 49-carrier graphs by a rule fixed **before** reading
the desired exponent. It must also justify link cost, ruler and radial
response separately if the eventual target is measured gravity.

## 6. Certification and scope

The theorem U39 is a finite impossibility for derivation from **the stated
current admissibility rules**; it is not a theorem that no richer uncut law
could select the desired graph. U40 is a conditional combinatorial derivation
of leading quadratic capacity. U41 is an exact finite counterexample to
silently equating aggregate channel count with a pointwise response. All
three results can be checked without importing a spacetime metric or a
classical gravitational field equation. No empirical data or universal
physical inverse-square law is inferred from their finite certificates.
