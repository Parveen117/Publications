# NM1 source audit and lineage

Audit base: Publications `1b4fcd80a0cdb7ac6ad8da28a8a48dbd64957ff5`.
This stage changes one repository only. It adds a new folder and its workflow;
all existing repository blobs are preserved.

## Consumed native sources

| Source | What NM1 consumes | What NM1 does not infer |
| --- | --- | --- |
| `../emk-ugd-algebra/README.md` | Primitive relations; existing winding, grading and ledger results | That every invisible datum is the same kind of memory |
| `../emk-ugd-algebra/certificates/emk1_determinant_seam_ladder.py` | Exact source matrices and matrix arithmetic as a comparison oracle | A diagonal K basis without an explicit change of basis |
| `../emk-ugd-algebra/CONNECTION_CALCULUS.md` | Associative sector, distinct curvature targets and observer-descent contract | That an algebraic section cocycle is automatically connection curvature |
| `../emk-ugd-algebra/UGD_ALGEBRA_CONTRACT.md` | Strict stored-state carry addition is not associative | Permission to import that digit operation as the present algebra addition |
| `../generalized-euler-evolution/GE4_MINIMAL_OBSERVER_MEMORY.md` | Observable recovery has a target and cannot reconstruct every hidden mode | That GE4's dynamical memory dimension counts this central sign |

These five sources are hash-pinned and checked against the audited Git tree.
The current tree was inspected for existing modular/mock-theta folder names;
the native algebra, connection and GE4 sources were read. This is a scoped
prior-work audit, not a claim to have reread every paper in every repository.
The full older source files and certificates remain unchanged.

The new target differs from GE4: it recovers the two matrix lifts of a
projective group element, not a heat generator from response moments. Existing
native winding and ledger ideas motivate it; they are not reclassified as new.

## Classical sources and attribution

J. S. Milne, *Modular Functions and Modular Forms*, version 1.31 (2017),
[course notes](https://www.jmilne.org/math/CourseNotes/MF.pdf), theorem 2.12,
remark 2.14 and the modular-form definitions in section 4. These give the
standard modular generators, projective quotient and transformation laws.
NM1 supplies its own elementary Euclidean generation and sign-lift proofs.
Group extensions, sections, cocycles and automorphy factors are established
mathematical machinery. No priority claim is made for those facts.

Sander Zwegers, *Mock Theta Functions*, doctoral thesis (2002),
[arXiv:0807.4834](https://arxiv.org/abs/0807.4834). Its abstract describes
specific nonholomorphic corrections involving unary theta functions.
This is background for the next analytic gate, not evidence that the NM1 sign
ledger derives a shadow or a mock theta completion. No analytic theorem from
the thesis is used as an unproved step in NM1-T1 through NM1-T5.

## Choices and future gates

The existing native representation has swap K, and the rational S=I+R adapter
produces diagonal K. The lattice S Z^2 and the projective observer are stated
choices. Being available inside the algebra is weaker than being canonically
selected by its native recognition rules.

To advance toward a nontrivial analytic result, choose a native history
function, justify its convergence/domain, derive its transformation defect,
and fix an admissible correction class with existence and normalization.
A finite sign correction must not be substituted for the distinct analytic
obstructions in mock-modular theory. None of these open gates is marked solved.
