# UG-A1: universal-generator draft audit and framework map

9 October 2026. This updates the framework explanation, not the canonical operator engine. The uploaded draft is not republished. T41 already supplies the corrected cut-graded generator theorem:
https://github.com/Parveen117/Recognition-Kernel-Framework/blob/3cc5a33b05c16d59c90994ddda69dedc0d392424/theorum/41_cut_graded_universal_generator_theorem.md

## Structural role

Given a lawful generator G and involution J, G_e=(G+JGJ)/2 and G_o=(G-JGJ)/2 uniquely separate sector-preserving and cross-sector action. The exponential U_t generates evolution. If G is cut odd, JU_tJ=U_-t; its even/odd readings reconstruct U and satisfy E_t^2-O_t^2=I. For general G the cut loop has log(JU_tJU_t)=2tG_e+t^2[G_e,G_o]+O(t^3). This is T41's seam-order coefficient, not automatically NT-3 response curvature.

The recognition split S=R+D needs its event lift and observer cut. T41's generator split is a different decomposition. They connect through a specified flow observer, not by equating their similarly named parts. Likewise shared-record heat needs its averaging/refinement protocol; a single exponential does not select a heat law.

## Draft disposition

| Uploaded claim | Correct treatment |
| --- | --- |
| Displayed infinite coordinate drift is universal | No topology/domain is supplied; angle theta on S1 needs a lift to be a globally real coefficient. On a finite lifted chart the stated drift moves kappa, phi and Sigma, keeping theta and Gamma fixed. It does not generate the later Gamma growth or memory decay. |
| Setting coordinates to zero gives arbitrary D_v | False: the remaining partial_kappa is not an arbitrary vector field; kappa=0 is not invariant under its flow. Restriction requires an invariant carrier and intertwining map. |
| Derivative equals upper-bidiagonal difference matrix | Requires an explicit discretization and limit; not an exact derivative identity. |
| G is diagonal with winding/growth/decay eigenvalues | Not derived from the displayed drift. Integer Fourier labels belong to periodic coordinates, not automatically a real winding coordinate. |
| Both cut/join products have different values | False. Both equal exp(G)-exp(-G); join^2-cut^2=4I. These derived channels are not the primitive involution J. |
| One-step Eye equals ker G | False in resonant sectors. All-time fixed vectors equal ker G for a differentiable group. A one-step fixed vector may have nonzero frequency 2pi. |
| Negative spectrum gives a nonzero Eye attractor | Strict decay gives zero. A nonzero limit needs a zero sector and decay on its complement, under specified hypotheses. |
| 8 lambda=2pi i m gives both closure signs | Only + closure. Minus closure requires odd multiples of pi i. Whole-operator closure is stronger than closure on one eigenvector. |
| Hyperoperator labels automatically form spectral decomposition | Completeness/invariance must be proved. Orthogonal complete component rows preserve G*G; arbitrary named parts need not. |
| G_(n+1)=2sinh(G_n) converges to zero | False for nonzero real scalars: magnitude grows by more than factor two. Zero is repelling. +/-i pi are not fixed points either. |
| Trace(Eye)=R/lambda | No typed projector, finite trace, or cost theorem supplied; T41 uses target-relative decoder burden. |
| Time-dependent attenuation is exp(tG) on the original state | Generally fails composition. Can be made autonomous with an additional decaying rate state; see below. |

## A constructive repaired cost/memory flow

On a finite lifted chart introduce coordinates (kappa,phi_1,...,phi_m,r, Sigma,memory,z,amplitude), r=log Gamma, and fixed real rates omega_j,alpha,b,gamma,lambda. Define the smooth vector field

    G = partial_kappa + sum omega_j partial_phi_j + alpha partial_r
        + b partial_Sigma - gamma memory partial_memory
        - lambda z partial_z - z amplitude partial_amplitude.

Its complete real flow is

    kappa(t)=kappa0+t; phi_j(t)=phi_j0+omega_j t;
    r(t)=r0+alpha t; Sigma(t)=Sigma0+b t;
    memory(t)=memory0 exp(-gamma t); z(t)=z0 exp(-lambda t);
    amplitude(t)=amplitude0 exp[-z0(1-exp(-lambda t))/lambda], lambda!=0.

At lambda=0 use amplitude0 exp(-z0 t). Gamma(t)=Gamma0 exp(alpha t). The integrated-rate identity
c_lambda(t+s)=c_lambda(t)+exp(-lambda t)c_lambda(s) proves composition because the second step starts with the updated z, not z0. All these equations follow by solving the stated scalar ODEs; local pullback differentiation gives the generator. This is a declared finite-dimensional model of the draft's desired flow, not a derivation of physical rates, primitive geometry or universality. Positive Gamma is implemented by its log coordinate; angles may be quotiented modulo 2pi after the lifted construction. Real kappa is not an integer winding label until a separate return/ledger contract is supplied.

## Certificate scope

`python3.12 verify.py` runs seven exact symbolic/rational groups and five explicit false-identity rejection controls. It recomputes RESULT.json with source hashes. This freshly verifies the stated finite controls; it does not rerun T41's historical CI or certify all bounded/unbounded conclusions there. T41 has written proofs and a source-reported exact/CI certificate. Its original Hilbert-space formulation is preserved as such; a native-first completion theorem is a separate bridge, not earned by renaming symbols.
