# REC-1 — Relational event clock and a derived comparison invariant

## 1. Question and declared sector

Can comparison between observable events supply dynamics, and does obtaining a
45-degree diagonal determine the transport law? We answer in an actual reduction
of LNC-1, rather than selecting a new operator model to enforce a diagonal.

The result is classical and sector-specific. Select periodic canonical fields on
a spatial circle with coordinate period one, homogeneous initial fields,
homogeneous densitized lapse N(s), zero shift, and vanishing matter psi=0.
Use LNC's independent canonical coframe field e>0 and any C1 potential V(phi).
Periodic fields remove the spatial boundary terms in LNC's bulk calculation.
Homogeneity and psi=0 are preserved by its equations with this lapse and shift.
This is a consistent homogeneous solution sector, not a reduction of every
inhomogeneous gauge orbit. It removes the propagating matter sector deliberately.

The spatial proper length is L=int e dx=e in this chart, and its conjugate is
the homogeneous pi. The homogeneous scalar momentum is k=p_phi. The restricted
symplectic potential and constraint are

    pi dL + k dphi,
    C = L pi k + L^2 V(phi) = 0,
    S = int [pi L_dot + k phi_dot - N C] ds.

Thus {L,pi}=1 and {phi,k}=1. The spatial constraint vanishes on this sector.
These conventions follow the selected LNC action, not a physical unit calibration.

**Global cut contract.** This periodic independent-coframe model is an additional
boundary choice. It does not assert a globally periodic single-angle cut lift.
If one also imposes LCF's oriented e=theta_x/2 with a smooth single circle-valued
angle of fixed winding n, then L=pi_math n is fixed. A changing L below cannot
satisfy that extra global contract. Local cut charts still exist; gluing their
global records is a separate problem. Here pi denotes canonical momentum and
pi_math denotes the circular constant.

## 2. Clock reduction theorem

On a connected branch with L>0, pi nonzero and N nonzero, phi is locally a clock:

    phi_dot = N L pi != 0.

Fix its reading to tau, solve C=0 for k, and substitute in the action:

    k = -L V(tau)/pi,
    S_red = int [pi dL/dtau - H_red] dtau,
    H_red(L,pi;tau) = L V(tau)/pi.

The gauge condition chi=phi-tau has {chi,C}=L pi, so the constraint/gauge pair is
regular on this branch. The Dirac bracket of L and pi remains one: both commute
with chi, so both correction terms vanish. No time-independent or positive
Hamiltonian is inferred for arbitrary V and both signs of pi.

Hamilton's equations are

    dL/dtau = -L V(tau)/pi^2,
    dpi/dtau = -V(tau)/pi.                         (REC.1)

For an independent check, the unreduced equations are

    L_dot = N L k,                 phi_dot = N L pi,
    pi_dot = -N(pi k + 2L V),     k_dot = -N L^2 V'.

Dividing by phi_dot and using C=0 gives REC.1. Also differentiating
k=-L V/pi with REC.1 gives dk/dtau=-L V'/pi, the remaining unreduced equation
divided by phi_dot. Thus no momentum equation was discarded. A nonzero change
of N changes the parameter rate but cancels from these comparison derivatives.
This cancellation applies to the selected homogeneous gauge family; it is not a
proof of full local spacetime gauge invariance.

## 3. Two conserved records and complete observables in the sector

Let W be any primitive of V on the clock interval, W'=V. Define

    K = L/pi,                M = pi^2 + 2 W(phi).

Here M is a mathematical conserved record; no physical mass identification is
asserted. REC.1 gives dK/dtau=0 and dM/dtau=0. More strongly, on the unreduced
four-dimensional homogeneous phase space,

    {K,C} = 2 C/pi^2,
    {M,C} = -(2 pi/L) C,
    {K,M} = 2.                                      (REC.2)

Hence K,M are weak Dirac observables for this single homogeneous constraint.
Proof: differentiate the displayed rational functions in the two canonical
pairs; the first two brackets factor by C. The third equals (1/pi)(2pi)=2.
The native word proofs below certify the denominator-cleared identities;
L>0 and pi nonzero remain explicit analytic hypotheses.

At a specified clock value tau, define on a fixed sign branch sigma=sign(pi)

    Pi_tau = sigma sqrt(M - 2 W(tau)),
    L_tau = K Pi_tau,                              (REC.3)

where the radicand is strictly positive and K Pi_tau>0. These functions commute
weakly with C, agree with pi,L at phi=tau, and satisfy {L_tau,Pi_tau}=1:
{K,Pi_tau}=1/Pi_tau follows from {K,M}=2. They are complete observables for the
homogeneous constraint, not for the full local LNC constraint family.

Equivalently, starting at tau0,

    pi(tau)^2 = pi0^2 - 2[W(tau)-W(tau0)],
    L(tau) = L0 pi(tau)/pi0.                       (REC.4)

These expressions solve all reduced equations uniquely on the chosen branch.
For V C1 their vector field is locally Lipschitz away from pi=0; the explicit
formulas exhibit the solution until the radicand reaches zero or the prescribed
clock interval ends. They do not extend the chart through a singular clock.

## 4. What the 45-degree comparison means here

Use positive normalized readings a=L/L0 and b=pi/pi0. The sign of pi cannot
change on the branch. REC.4 implies a=b and the dimensionless comparison

    I = (L/L0)/(pi/pi0) = 1.                       (REC.5)

Thus (a,b) lies on the diagonal, and so does (log a,log b) wherever defined.
With equal Euclidean axis calibration the nonstationary diagonal has angle
45 degrees to either axis. Independent axis rescaling changes the plotted slope.
The invariant is the ratio in REC.5, not that particular plotted angle.

The shared relative rate is

    d log a/dtau = d log b/dtau = -V(tau)/pi(tau)^2.

Consequently a=exp(int -V/pi^2 dtau) on the branch. This is the integrated scalar
rate equation, not a derivation of a constant-rate semigroup or a quantum law.
For constant V=mu the rate is -mu/[pi0^2-2mu(tau-tau0)], generally nonconstant.
For V=0, pi and L are constant even though phi is a valid running clock. Then
the comparison point has no nonzero tangent direction: invariance needs no motion
toward 45 degrees. This reuses IRT-1's calibration theorem, rather than reproving
its general statement as a new result. The new result is its derivation inside
the LNC reduced equations and the corresponding sector Dirac observables.

## 5. What event readings can and cannot identify

If pi(tau) is known differentiably on a clock interval, its rate determines

    V(tau) = -(1/2) d[pi(tau)^2]/dtau.             (REC.6)

For two readings it determines only an integral,

    int_(tau_a)^(tau_b) V(u) du
      = [pi(tau_a)^2-pi(tau_b)^2]/2.               (REC.7)

An invariant or an endpoint pair alone therefore cannot select the dynamics.
Here is an exact example, not an approximation or a fit. Take tau in [0,4],
tau0=0, L0=6, pi0=3 and choose either

    Model A: V_A=1,       W_A=tau,
    Model B: V_B=tau/2,   W_B=tau^2/4.

Both have K=2, M=9, initial (L,pi)=(6,3) and final (L,pi)=(2,1).
Both remain regular with positive pi on this interval and obey a=b throughout.
At tau=2 their pi^2 readings differ: 5 versus 7 (and L^2=20 versus 28).
They also have different initial k fixed by their respective constraints.
Thus this is ambiguity of the specified partial event record, not equality of
all canonical initial data. An additional k reading at the initial event can
distinguish these two models; an interior pi reading does too. No full-state
nonuniqueness for a fixed action is claimed.

Each additional endpoint difference gives another integral of V. Recovering a
pointwise potential from finitely many such integrals requires a chosen function
class or more data. No noise-stability or statistical inverse theorem is supplied.
The straight segment connecting two plotted readings is a visualization of the
comparison, not proof that the physical trajectory in spacetime is straight.

## 6. Failure conditions and the closure distinction

* If pi=0, {phi,C}=0 and division by the clock rate is invalid. At L=1, pi=0,
  k=0, V=0 the constraint is satisfied but this event-clock fails. This does not
  make the entire theory inconsistent; another clock or patch may be needed.
* For finite nonzero K, radicand zero in REC.3 also sends L to zero. This exits
  the positive-coframe sector; no smooth global continuation is asserted here.
* If a fixed-winding global single-angle lift is additionally required, dL/dtau=0
  is demanded, so REC.1 forces V=0 on a regular branch. The varying-potential
  examples do not satisfy that stronger global cut assumption.
* Algebraic closure is not measurement completeness. Exact Pauli operators X,Z
  generate a closed finite matrix algebra although [X,Z]!=0. In the +1
  eigenstate of Y, Var(X)=Var(Z)=1 and |<[X,Z]>|^2/4=1. Algebraic closure and
  nonzero uncertainty coexist. Conversely Var(Z)=0 in a Z eigenstate inside
  the same algebra. Neither uncertainty nor noncommutativity forbids algebraic
  closure. Calling a measurement record 'closed' needs a separate operational
  definition; it cannot replace these algebraic statements.

The finite Pauli calculation is only a control of the terminology. It is not a
quantization of REC.1 or proof that the local LNC quantum constraints close.

## 7. Certification, provenance and next gate

The written theorem proves the interval and branch statements. The existing
native rational/cut engine checks constrained jets, clock cancellation, exact
polynomial potential reconstruction, finite event ambiguity, and the Pauli
control. Native rewriting checks denominator-cleared weak invariants and the
comparison identity with serialized proof replay and a corrupted-proof control.
These certificates do not mechanically formalize continuum existence, square
roots or logarithms; those arguments and their domains are written above.

Relational clocks and complete observables are established methods: see
B. Dittrich, *Partial and Complete Observables for Hamiltonian Constrained
Systems*, [gr-qc/0411013](https://arxiv.org/abs/gr-qc/0411013). That work develops
the general constrained-system construction. The formulas here are derived
directly from the selected LNC homogeneous action; no novelty claim is made for
the general clock method or dilaton-type reduction.

Next: represent or quantize the full local constraints, beginning with a
declared massless/zero-potential diagnostic if useful; calculate central terms,
operator domains and physical inner products. REC-1 supplies an exact reduced
classical comparison target. It does not bypass inhomogeneous reduction,
matter quantization, global cut gluing, or physical 3+1 and empirical obligations.
