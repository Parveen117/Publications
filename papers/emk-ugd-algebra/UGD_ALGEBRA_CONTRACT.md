# UGD-A0 — The associative-carrier gate

Author: Monty Dabas. 2 October 2026.

This is a scoped mathematical boundary result for the existing [UGD-1 numeral implementation](certificates/ugd1_numerals.py), not a replacement number system or a rejection of its certified carry-conservation theorem. It identifies the equality contract needed before importing that numeral container into an associative algebra or tensor calculus.

## A0.1. Strict carry addition is not associative

Use the source's strict equality: phase alphabet, retained digit map and ledger must all agree. At one scale, take phase exponent zero, \(K=4\), zero incoming ledgers, and seam digits
\[
x=-1,\qquad y=1,\qquad z=1.
\]
The source carry rule gives:
\[
(x\oplus y)\oplus z=(s=1,\ell=0),\qquad
x\oplus(y\oplus z)=(s=0,\ell=1).
\]
Both have total seam charge \(1\), but they are different stored UGD states.

**Proof.** The first grouping adds \(-1+1=0\), then \(0+1=1\), without overflow. The second grouping first adds \(1+1\), retains \(s=1\) and ledgers the overflow \(\ell=1\); adding \(-1\) then leaves \(s=0,\ell=1\). Strict source equality distinguishes the results. A single valid counterexample disproves universal associativity of this operation under that equality. No phase carry or approximation is involved. ∎

## A0.2. Charge conservation survives every grouping

At one scale, let \(s_1,s_2\in\{-1,0,1\}\) and let \(\ell_1,\ell_2\) be integer ledgers. The source rule is
\[
s'=\operatorname{clamp}(s_1+s_2,-1,1),\qquad
\ell'=\ell_1+\ell_2+s_1+s_2-s'.
\]
Therefore, with \(q=s+\ell\),
\[
\boxed{q(x\oplus y)=q(x)+q(y).}
\]

**Proof.** Add the displayed expressions for \(s'\) and \(\ell'\). For finite numerals, sum the same identity across scales and include the incoming global ledger. Extra phase-carry slots carry no seam charge, so they do not change the equality. Repeated addition then conserves the sum of all input charges for every grouping, by induction. ∎

Thus a quotient that observes only \(q\) has ordinary associative integer addition. That quotient is target-faithful for **total seam charge**. It is not faithful for the full retained digit/ledger target: A0.1 already supplies two different states it identifies.

Deleting the ledger is not a repair. The \((1,1)\) addition retains \(s=1\) and ledgers another \(1\); deleting that overflow changes charge \(2\) into \(1\).

## A0.3. Consequence for the mathematical programme

An associative algebra over a field has associative addition. Consequently the strict numeral carrier, equipped with this carry addition, cannot simply be declared that algebra without changing its operation or equality contract.

This says nothing against other existing native associative sectors: [EMK-C1](CONNECTION_CALCULUS.md) uses the separately established native operator algebra. It also does not assert that every possible UGD realization is nonassociative.

The next UGD algebra theorem must choose and prove one of the appropriate contracts:

- a declared target quotient, with proof that its product and all required operators descend;
- an associative history/operator lift, with a readout that retains the required seam/ledger information;
- a revised normalization/recognition law, derived and certified before replacing the historical one.

These are mathematical alternatives, not permission to erase ledger data merely because total charge matches. A full algebra still needs product associativity, distributivity, units, typing and compatibility with the chosen recognition relation. Linked-seam composition and infinite numerals retain UGD-1's existing open/declaration boundaries.

## Evidence and reproduction

[UGDA0_RESULT.json](certificates/UGDA0_RESULT.json) records the exact parenthesization witness, all 27 balanced three-digit cases and 225 seam/ledger pair-conservation checks. Its [source pins](certificates/UGDA0_SOURCE_PINS.json) bind the unchanged historical implementation, this proof and the audit adapter.

~~~bash
python3.12 -B papers/emk-ugd-algebra/certificates/ugd_algebra_contract.py --check
~~~

The general no-go follows from the exact counterexample; the general charge identity is proved above. Finite execution checks the source behavior and the published evidence. It does not certify the full UGD algebra, a proof assistant or physical predictions.
