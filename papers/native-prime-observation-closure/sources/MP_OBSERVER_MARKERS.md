# Observer-marked spectral separation theorem

## 1. Scope

Let

\[
\mathcal A=M_d(\mathbb C)
\]

with trace pairing

\[
\langle M,D\rangle_{\mathrm{tr}}=\operatorname{Tr}(MD).
\]

All claims in this note are finite-dimensional. No trace-class, Fredholm,
unbounded-domain or gauge-quotient statement is imported silently.

## 2. Spectral blindness inherited from the preceding gate

For invertible \(X,Y\),

\[
YX=X^{-1}(XY)X.
\]

Hence

\[
\det(I-qYX)=\det(I-qXY).
\]

The characteristic determinant therefore cannot distinguish the two-factor
orders, although the endpoints may differ.

Set

\[
U_{YX}=YX,\qquad U_{XY}=XY,\qquad D=U_{YX}-U_{XY}.
\]

The ordinary spectral shadow may vanish while \(D\neq0\).

## 3. Marked trace response

### Definition 3.1

For a marker \(M\in\mathcal A\), define

\[
\tau_M(U)=\operatorname{Tr}(MU).
\]

For two endpoints \(U,V\), define the marked residue

\[
G_M(U,V)=\tau_M(U)-\tau_M(V)=\operatorname{Tr}(M(U-V)).
\]

For a finite marker bank \(\mathcal M=\{M_a\}_{a=1}^r\), define

\[
G_{\mathcal M}(U,V)
=
\bigl(G_{M_a}(U,V)\bigr)_{a=1}^r.
\]

### Theorem 3.2 — two-factor marked detection

For \(U_{YX}=YX\) and \(U_{XY}=XY\),

\[
G_M(U_{YX},U_{XY})
=
\operatorname{Tr}(M(YX-XY)).
\]

Equivalently, by cyclicity of trace,

\[
G_M(U_{YX},U_{XY})
=
\operatorname{Tr}([M,Y]X)
=
\operatorname{Tr}(M[Y,X]).
\]

Therefore

\[
\operatorname{Tr}(M[Y,X])\neq0
\quad\Longrightarrow\quad
M\text{ detects the order reversal}.
\]

Noncommutation alone is not sufficient for a prescribed marker: the trace
pairing may still vanish.

## 4. Complete finite marker banks

### Theorem 4.1 — matrix-unit completeness

Let \(E_{ij}\) be the standard matrix units. Use the marker bank

\[
\mathcal M_{\mathrm{unit}}
=
\{E_{ji}:1\le i,j\le d\}.
\]

Then for every \(D\in M_d(\mathbb C)\),

\[
\operatorname{Tr}(E_{ji}D)=D_{ij}.
\]

Consequently,

\[
G_{\mathcal M_{\mathrm{unit}}}(U,V)=0
\quad\Longleftrightarrow\quad
U=V.
\]

Moreover the endpoint difference is reconstructed exactly entry by entry.

### Corollary 4.2 — finite separation

For every nonzero \(D\in M_d(\mathbb C)\), there exists a marker \(M\) such that

\[
\operatorname{Tr}(MD)\neq0.
\]

This is simply nondegeneracy of the trace pairing.

### Theorem 4.3 — Hermitian marker completeness

The following \(d^2\) Hermitian matrices form a complete marker bank:

\[
E_{ii},
\]

\[
S_{ij}=E_{ij}+E_{ji},
\qquad
A_{ij}=i(E_{ij}-E_{ji}),
\qquad i<j.
\]

Their responses are

\[
\operatorname{Tr}(E_{ii}D)=D_{ii},
\]

\[
\operatorname{Tr}(S_{ij}D)=D_{ji}+D_{ij},
\]

\[
\operatorname{Tr}(A_{ij}D)=i(D_{ji}-D_{ij}).
\]

Thus

\[
D_{ij}
=
\frac{
\operatorname{Tr}(S_{ij}D)
+i\operatorname{Tr}(A_{ij}D)
}{2},
\]

\[
D_{ji}
=
\frac{
\operatorname{Tr}(S_{ij}D)
-i\operatorname{Tr}(A_{ij}D)
}{2}.
\]

Hence the Hermitian bank reconstructs every complex endpoint difference.

## 5. Marked determinants

### Definition 5.1

For \(M,U\in M_d(\mathbb C)\), define

\[
F_M(q;U)=\det(I-qMU).
\]

### Lemma 5.2 — linear coefficient

The first coefficient is

\[
[q]F_M(q;U)=-\operatorname{Tr}(MU).
\]

This follows from the standard first variation of the determinant at the
identity.

### Theorem 5.3 — marked determinant separation

If

\[
\operatorname{Tr}(M(U-V))\neq0,
\]

then

\[
F_M(q;U)\neq F_M(q;V)
\]

as polynomials in \(q\), because their linear coefficients differ.

Thus a marker can restore order information erased by the ordinary
characteristic determinant.

## 6. Similarity covariance and observer relativity

Let \(S\) be invertible and transform

\[
U\mapsto SUS^{-1},
\qquad
M\mapsto SMS^{-1}.
\]

Then

\[
\operatorname{Tr}\!\left(
SMS^{-1}SUS^{-1}
\right)
=
\operatorname{Tr}(MU),
\]

and

\[
\det\!\left(
I-qSMS^{-1}SUS^{-1}
\right)
=
\det(I-qMU).
\]

Therefore marked traces and marked determinants are covariant under
simultaneous transport of the operator and observer frame.

If \(M\) is held fixed while \(U\) is conjugated, the response may change.
A fixed marker is therefore observer-relative and is not automatically a
gauge-invariant observable.

## 7. RSC recognition classification

For a declared marker bank \(\mathcal M\), define the observer recognition map

\[
\mathcal R_{\mathcal M}(U)
=
\bigl(\operatorname{Tr}(M_aU)\bigr)_a.
\]

For two endpoints:

- \(U=V\): exact endpoint closure;
- \(U\neq V\) and \(G_{\mathcal M}(U,V)=0\): marker-bank blindness;
- \(G_{\mathcal M}(U,V)\neq0\): observer-visible residue;
- ordinary determinants equal but \(G_{\mathcal M}\neq0\): repaired spectral
  blindness.

The marker bank must be declared before the comparison. Choosing a marker
post hoc merely to force a desired distinction is not a lawful RSC ledger.

## 8. Negative controls

### 8.1 Scalar blindness

For \(M=cI\),

\[
\operatorname{Tr}(M(YX-XY))
=
c\operatorname{Tr}([Y,X])
=
0.
\]

Scalar markers do not repair the two-factor trace shadow.

### 8.2 Incomplete-bank blindness

A diagonal marker bank misses every purely off-diagonal endpoint residue.
Completeness is not inherited from the word “observer.”

### 8.3 Endpoint equality does not recover history

Let \(N^2=0\). Then

\[
(I+N)(I-N)=I=I\cdot I.
\]

Every endpoint marker gives the same response, although the factor words are
different. A complete endpoint bank reconstructs the endpoint operator, not
the original history.

### 8.4 No robustness without a margin

Exact separation may be arbitrarily small after rescaling the residue.
A quantitative noisy-observer theorem requires a lower separation margin and
conditioning estimates. That gate remains open.

## 9. Exact calibration

The executable model uses

\[
X=I+E_{12},
\qquad
Y=I+E_{23}.
\]

Then

\[
YX-XY=-E_{13}\neq0,
\]

while ordinary trace and characteristic determinant differences vanish.

The marker

\[
M=E_{31}
\]

gives

\[
\operatorname{Tr}(M(YX-XY))=-1
\]

and

\[
\det(I-qMYX)-\det(I-qMXY)=q.
\]

A Hermitian marker \(E_{13}+E_{31}\) gives the same nonzero trace response.

The package also verifies exact reconstruction by matrix-unit and Hermitian
banks, incomplete diagonal-bank blindness, simultaneous similarity covariance,
fixed-marker observer relativity and the endpoint-history no-go result.

## 10. Proof status

```text
MARKED_TRACE_ORDER_DETECTION = PROVED
FINITE_MATRIX_UNIT_BANK_COMPLETENESS = PROVED
FINITE_HERMITIAN_BANK_COMPLETENESS = PROVED
MARKED_DETERMINANT_LINEAR_SEPARATION = PROVED
SIMULTANEOUS_SIMILARITY_COVARIANCE = PROVED
FIXED_MARKER_GAUGE_INVARIANCE = FALSE IN GENERAL
INCOMPLETE_BANK_SUFFICIENCY = FALSE
ENDPOINT_MARKERS_RECONSTRUCT_FACTOR_WORD = FALSE
NOISE_STABILITY = OPEN
INFINITE_DIMENSIONAL_MARKER_COMPLETENESS = OPEN
FREDHOLM_MARKED_DETERMINANT_EXTENSION = OPEN
```

No native RH, Yang–Mills or other Millennium-problem conclusion is asserted.
