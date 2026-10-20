# EMK-F2 — Discrete observation, multiplicative transport and prime index

Author: Monty. Written development: 21 October 2026.

Continuation of [EMK-F1](FIBONACCI_TRANSPORT_CURVATURE.md).
Source audit: [EMK algebra](README.md) and
[UGD-A0](UGD_ALGEBRA_CONTRACT.md), inspected at
12534fc34ebe8aa15ddc3e8147121640ef79fb87.

## Scope

This note gives a faithful multiplicative seam transport, and proves a
conditional geometric characterization of primes after choosing an
observation lattice. These are elementary algebra and lattice/group
results in the specified representation, not new theorems of prime
distribution. No claim is made that flatness selects primes.

The source audit matters: UGD-A0 starts with integer ledgers and derives
integer charge addition. It does not derive integer quantization from
a continuous EMK carrier. Moreover its strict stored-state carry
addition is not associative; its additive charge quotient forgets
stored history. The lattice model below is an explicitly separate
associative adapter, not a silent replacement for full UGD memory.

## F2.1. Multiplication-compatible seam transport

In the real matrix representation used in F1, take
K=diag(1,-1) and R=[[0,-1],[1,0]]. Set
\[
P_\pm=(I\pm K)/2,\qquad T_x=xP_++P_-,\qquad x>0.
\]
Because the projectors are orthogonal and sum to I,
\[
T_xT_y=T_{xy},\quad T_1=I,\quad T_x^{-1}=T_{1/x}.
\]
In this representation T_x=diag(x,1), so the assignment is injective.
It also has the exact continuous interpolation
\[
T_x=\exp((\log x)P_+).
\]
All these transports commute with K and with one another. This is a
multiplicative seam family, not F1's Fibonacci family Q_a.

In the full positive-real family, every x>1 factors as
\(T_x=T_{\sqrt{x}}T_{\sqrt{x}}\), with both factors nonidentity.
Prime irreducibility therefore requires an admissibility restriction.

For completeness, F1's transfer matrices instead satisfy
\[
Q_aQ_b=\begin{pmatrix}ab+1&a\\b&1\end{pmatrix}\ne Q_{ab}.
\]
They cannot be substituted for T_x in a prime-factorization argument.

## F2.2. Lattice preservation selects integers, conditional on a lattice

Choose the full-rank observation lattice
\(\Lambda=\mathbb Z e_++\mathbb Z e_-\), where
\(e_+=(1,0)\) and \(e_-=(0,1)\). Declare a positive-real transport
admissible when \(T_x\Lambda\subseteq\Lambda\). Then
\[
\boxed{T_x\Lambda\subseteq\Lambda\iff x\in\mathbb Z_{>0}.}
\]

Proof. Applying T_x to e_+ gives (x,0), which belongs to Lambda exactly
when x is an integer. Applying it to e_- leaves e_- fixed.
These conditions are also sufficient on every integer combination.
Positivity excludes zero and negative integers. QED.

This is endomorphism admissibility, not lattice automorphism
admissibility. Requiring equality \(T_x\Lambda=\Lambda\) would leave
only x=1 in this positive family. The nontrivial transports are
invertible over the real carrier but not invertible as lattice
endomorphisms.

The lattice is supplied structure. The theorem does not obtain a
discrete set from the real algebra without that additional input.

## F2.3. Prime observations as simple finite quotients

For an integer n>=1,
\[
T_n\Lambda=n\mathbb Z e_++\mathbb Z e_-,\qquad
\Lambda/T_n\Lambda\cong\mathbb Z/n\mathbb Z.
\]
The isomorphism sends the class of (u,v) to u modulo n: it is
surjective and its kernel is exactly T_n Lambda. In particular,
\[
[\Lambda:T_n\Lambda]=n=\det T_n.
\]

For n>1 the following are equivalent:

1. n is prime.
2. T_n is not a product T_a T_b of two admissible nonidentity
   transports from this family.
3. The additive group Lambda/T_n Lambda is nontrivial and simple.
4. No additive subgroup M lies strictly between T_n Lambda and Lambda.

Proof. Faithfulness and F2.1 identify (2) with the absence of integer
factorizations n=ab with a,b>1. If n is prime, every nonzero element
of Z/nZ generates that group: its integer representative is coprime
to n, and Bezout's identity supplies a multiplicative inverse modulo n.
If n=ab with a,b>1, the multiples of a form a nonzero proper subgroup
of Z/nZ. This proves equivalence with (3). Taking preimages under the
quotient map gives the correspondence in (4). QED.

For example, the index-five quotient has no proper nontrivial subgroup.
At n=6 there is the strict inclusion chain
\[
T_6\Lambda\subsetneq T_2\Lambda\subsetneq\Lambda,
\]
whose indices are 3 and 2. More generally, for n=ab,
\[
[\Lambda:T_a\Lambda]=a,\quad
[T_a\Lambda:T_{ab}\Lambda]=b.
\]

This quotient is a finite algebraic record of residue classes. Calling
it physical memory, a probability distribution, or entropy requires
additional definitions. Its cardinality does not by itself specify
a state probability.

The lattice index is also not a Fredholm index: T_n is an invertible
operator on the real two-dimensional carrier, with zero kernel and
cokernel there, hence Fredholm index zero. Changing from a real vector
space to an integer lattice changes the target.

## F2.4. The primitive symmetries do not select the observation lattice

Both the standard lattice Lambda and
\[
\Lambda'=\mathbb Z(1,1)+\mathbb Z(1,-1)
=\{(u,v)\in\mathbb Z^2:u\equiv v\pmod2\}
\]
are preserved bijectively by K and R.

Proof. On the displayed basis v_1=(1,1), v_2=(1,-1),
K swaps v_1 and v_2, while Rv_1=-v_2 and Rv_2=v_1.
The same invariance on Lambda follows directly from signed coordinate
permutations. QED.

However their T_x admissibility differs. On Lambda' the matrix is
\[
[T_x]_{\Lambda'}=
\frac12\begin{pmatrix}x+1&x-1\\x-1&x+1\end{pmatrix}.
\]
Its entries are all integers exactly when x is an odd integer.
Consequently
\[
T_x\Lambda'\subseteq\Lambda'
\iff x\in\{1,3,5,\ldots\}.
\]
The standard lattice admits T_2, while Lambda' does not, although
both satisfy the primitive K/R symmetry requirement.

This is an explicit obstruction to selecting the desired observation
arithmetic from K/R lattice invariance alone. It does not show that
no stronger native axiom could select a lattice. Such an axiom must be
stated and justified separately.

A general full-rank lattice preserved by T_x necessarily has integer
x: in a lattice basis its endomorphism matrix has integer entries and
determinant x. Sufficiency for every positive integer x is specific to
the chosen lattice, as Lambda' demonstrates.

## F2.5. Joint observer changes preserve the conditional arithmetic

For any invertible real S, transform all data together:
\[
\Lambda_S=S\Lambda,\quad T_{n,S}=ST_nS^{-1},\quad
K_S=SKS^{-1},\quad R_S=SRS^{-1}.
\]
Then S induces an isomorphism
\[
\Lambda/T_n\Lambda
\longrightarrow
\Lambda_S/T_{n,S}\Lambda_S,\qquad [v]\mapsto[Sv].
\]
It is well defined and bijective because S maps T_n Lambda exactly
onto T_{n,S} Lambda_S. Thus the index and the prime simplicity
criterion are invariant under these joint changes of representation.

Changing the lattice while keeping the transport family fixed is a
change in the observation model, not merely this joint observer change.
That is precisely what the comparison in F2.4 tests.

## F2.6. Relation to the native ledger and curvature

UGD-A0 proves \(q(x\oplus y)=q(x)+q(y)\) for the supplied integer
charge q=s+ledger. It also proves that this quotient loses distinctions
between retained digit/ledger histories. It does not supply:

- a faithful identification of full UGD states with Lambda;
- descent of the complete UGD product to this multiplicative adapter;
- a rule selecting Lambda rather than Lambda';
- a physical clock, length, state distribution or coupling.

Accordingly the lattice results above are conditional and must not be
advertised as deriving the UGD discreteness or all its memory.

Every T_n commutes with K, including composite n, so seam compatibility
does not distinguish primes. Complete connection flatness is a separate
target governed by EMK-C1, not merely this commutation condition.
For instance using any invertible T_n as u in
\(A_i=u^{-1}\delta_i u\) gives a pure-gauge flat connection for every n.
In that real algebraic connection construction, prime and composite
labels are equally flat; it is not a lattice automorphism assertion.

## Next obligation and gate accounting

Closed calculations: faithful multiplicative transport, conditional
integer admissibility, finite-index prime characterization and
invariance under joint observer changes.

Resolved limitation: K/R symmetry of a lattice does not uniquely select
the required integer observation family; two explicit compatible
lattices give different admissibility sets. Flatness alone does not
select primes in these constructions.

The next native task is to identify an already specified, history-aware
observation rule that selects an admissible carrier and prove its
compatibility with multiplication. No such rule is derived in this
note. Choosing the standard lattice remains an explicit assumption;
the finite quotient does not recover discarded UGD history.

No new physical prediction, RH positivity result, prime-distribution
theorem or quantum-gravity closure follows here. This is a written
representation and selection audit using elementary algebra.
No scientific computation, numerical simulation, certificate replay,
independent agents or workflow dispatch was performed. Existing
certificate verdicts do not certify this new note.
