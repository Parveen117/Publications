# NSB1 source audit and scope correction

The repository copies of both source PDFs were read directly. The source
manifest pins their bytes. The second paper's printed pages 28–29 were
also rendered to check the finite-section statement visually.

## What the existing papers already establish

| Source | Exact input reused | Limit retained |
| --- | --- | --- |
| SPECTRAL 1.pdf, *When Spectra Forget Order*, T5.1/T5.4 | Characteristic determinants lose two-factor reversal and cyclic order | Similarity does not identify endpoints or histories |
| Same, T15–17 | Complete marker banks and frame-controlled finite reconstruction | A bank must match its declared target and carrier |
| Same, T18.1 | Arbitrary-map target recovery iff observation fibres lie in target fibres | This statement is already nonlinear at the set level |
| Same, T18.2/T19.1 | Linear kernel criterion and endpoint-history obstruction | Nonlinear recovery cannot be justified by an undeclared linear kernel |
| SPECTRAL_2.pdf, T4.1/T7.2 | Minimum scalar lift count on a declared linear target model | NSB1 uses a finite jet-coefficient chart and a differentiable decoder |
| Same, T7.5 and completion sections | Quantitative target-decoder and uniform-margin requirements | Smooth injectivity alone does not give global or uniform stability |
| Same, Section 16 | Complement information is necessary for global operator sign | The generic restricted-eigenvalue tail bound needs the correction below |
| NT3/NT4/NT7 | Selected native response connection, Hessian intertwiner, typed process flatness | Different connections carry different curvature targets |
| NSR1-T2/T3/T4/T5 | Coupling/impedance curvature, wave reconstruction, symmetry constraints, one-point blindness | NSB1 strengthens the examples to entire spectral fields |
| NSR2 | Conjugate diagram transfer and the six lambda/z properties | Physical units, held variables and constitutive assumptions remain explicit |
| YM82 | Native SU(2) calculus, exact parent construction and Haar cap | NSB1's gauge example is a separate model, not the original continuum theory |

The new mathematical content within this programme is the common-filtration
spectral obstruction for nonlinear jet maps, its exact coefficient-marker
repair, the explicit all-bias constitutive examples, and the finite-jet
compact-gauge counterfamily. The usual Carleman triangular mechanism,
Poincare integrability and dual-flat Hessian geometry are acknowledged
inputs. No claim of external priority follows from an exact certificate.

## Required correction to SPECTRAL_2, Section 16

The printed source defines

\[
 \mu_n=\sup\sigma(P_nKP_n|_{P_nX})
\]

for a compact self-adjoint K, then states

\[
 \sup\sigma(K)\le\mu_n+\|K-P_nKP_n\|.
\]

That bound is not valid without an additional condition such as
mu_n>=0. Its zero-extended finite-rank operator has an extra zero
spectral value on the omitted complement.

For a proper finite-rank cut, the correct generic bound is

\[
 \boxed{\sup\sigma(K)
 \le\max\{0,\mu_n\}+\|K-P_nKP_n\|.}
 \tag{A1}
\]

Thus a sufficient certificate for I-K>=0 is

\[
 \boxed{\max\{0,\mu_n\}+\|K-P_nKP_n\|\le1.}
 \tag{A2}
\]

**Proof.** The zero-extended compression has spectral supremum
max(0,mu_n). Apply the self-adjoint norm perturbation inequality
sup sigma(A+B)<=sup sigma(A)+||B|| with
A=P_nKP_n and B=K-P_nKP_n. If P_n=I, use mu_n itself;
the displayed max bound remains a valid, weaker bound.

**Exact counterexample.** On l2 let

\[
 K=\operatorname{diag}(-100,2,0,\ldots)
\]

and let P project onto the first coordinate. Then mu=-100 and the
matched tail norm is 2. The printed expression gives -98<=1 and
would certify I-K>=0. In fact I-K has the negative eigenvalue -1.
The corrected expression is 2>1 and properly refuses that certificate.

The issue is confined here to the generic statement and Corollary 16.2
as printed. It does not invalidate Proposition 16.1's exact Schur
complement identity or the finite-dimensional blindness/lift theorems
used by NSB1. This audit does not assert that every numerical application
of the source used a negative finite supremum, and it does not replay
the entire RH programme.

The source PDFs remain unchanged. This explicit correction is published
beside the continuation and is tested adversarially, so the new work
does not silently inherit the faulty implication.

## Necessary distinctions in the nonlinear continuation

- A finite jet pullback is an exact operator on a finite local algebra.
  Its characteristic polynomial is not an infinite-dimensional Koopman
  spectrum and is not a quantum Hamiltonian spectrum.
- Characteristic-polynomial blindness does not imply equality of Jordan
  forms, singular values, full matrices or coordinate response markers.
- State jets at an equilibrium and time-response jets along a nonzero
  trajectory are different observations.
- The all-bias flat/curved pair has gauge-inequivalent selected connections.
  The exponential same-wave-speed pair has constant-fibre-conjugate
  connections and distinguishes calibrated tangent response geometry.
- Maxwell reciprocity expresses potential integrability. Arbitrary
  rotational, gauge or Onsager–Casimir symmetry does not force every
  selected curvature to vanish.
- Exactness of work, a dissipation potential, response-metric curvature
  and spacetime gauge curvature retain their separate physical meanings.
- A correct observation obstruction identifies missing information; it
  does not turn that information into a proved Yang–Mills mass gap.

## Primary literature checked

- Amini, Zheng, Sun and Motee, *Carleman Linearization of Nonlinear
  Systems and Its Finite-Section Approximations*,
  https://arxiv.org/abs/2207.07755.
- Shima, *Hessian manifolds of constant Hessian sectional curvature*,
  J. Math. Soc. Japan 47 (1995),
  https://www.jstage.jst.go.jp/article/jmath1948/47/4/47_4_735/_pdf.
- Mezincescu, Russo and Townsend, *Hamiltonian birefringence and
  Born–Infeld limits*, https://arxiv.org/html/2311.04278v3.
- Saito and Utsumi, *Symmetry in full counting statistics, fluctuation
  theorem, and relations among nonlinear transport coefficients in the
  presence of a magnetic field*, https://arxiv.org/pdf/0709.4128.

These references establish the lineage and physical conventions. The
written proofs and exact examples, rather than the citations alone,
carry the claims made in the continuation.
