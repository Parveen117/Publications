# PCR-1 observer meaning and replay

## What is observed

| Quantity | Meaning | How it enters the result |
|---|---|---|
| X_i=M(w_i) | Existing positive smeared density | m_i=<X_i> determines theta in the fixed family |
| P_i=D(k_i) | Existing Cameron--Martin current | p_i=<P_i> supplies beta_i=p_i/m_i |
| q_i=log m_i | Logarithm of the density mean | Coordinate in F=(1/h)sum dq_i wedge dp_i |
| H=Cov(X) | Density-response covariance | Determines the response tensor; independent of beta |
| T=M(S1) | Total density | Controls current domains; it is not time |

The pulse exp(-i beta.X/h) preserves every density probability while
changing current means. The two currents commute with one another, but
each fails to commute with its own density. The stated means are ensemble
expectations from separate preparations, not a single joint sharp reading.
Identification assumes known windows and the fixed prepared family.
Finite-sample reconstruction and a uniform inverse-conditioning bound
over the whole parameter space are not supplied.

## Native adapter

From the repository root:

```javascript
const PCR = require('./papers/native-phase-current-response/phase_current_adapter.cjs');
const QRT = require('./papers/native-quantum-response-transport/response_adapter.cjs');

const moments = QRT.moments(
  ['1/3', '1/3', '1/3'], [[1, 1], [2, 1], [1, 2]]);
const geometry = PCR.quantumGeometry(moments.H, '1/2');
const observed = PCR.readoutGeometry(
  moments.means, moments.H, ['3/5', '-2/7'], '1/2');

const source = QRT.finiteSource(
  ['1/3', '2/3', '2/3'], [[1, 1], [2, 1], [1, 2]], 3,
  'phase-current-example');
const path = PCR.phaseCurrentPath(source, {
  h: '1/2', c: '3/5', d: '-2/7',
  thetaPath: {a: 1, b: 1, p: 1, q: 2}
});
```

All coefficient arithmetic uses the frozen native F/Cut engine. Supply
exact rational strings and a real positive-definite 2-by-2 H. The
normalization satisfies 0<h<1, and logarithmic mean coordinates require
positive means. The phase parameters are real.

`quantumGeometry` returns the quantum tensor `QGT`, real metric `g`,
complex structure `J`, `BerryF`, `omega` and the Gram factorization.
The control order is (theta_1,theta_2,beta_1,beta_2). The complex tensor
has rank two; the real metric and Berry form have rank four.
The convention is A=i<psi,dpsi>, F=dA=-2 omega.

`readoutGeometry` returns `p`, `betaRecovered`, the (m,p) and (q,p)
Jacobians, their determinants, the canonical Berry form and its pullback.
It retains q=log(m) as a coordinate contract; it does not manufacture an
exact rational logarithm. It also does not supply matrices for P_i.

`phaseCurrentPath` uses beta(t)=(c t,d t²) and the supplied theta path.
It returns the phase `V`, state `psi`, lift `lift=V K`, its Hamiltonian
`liftGenerator`, Berry connection, horizontal derivative and norm, mean
and phase readout paths, and a Berry-horizontal scalar phase/state.
The comparison starts at beta(0)=0. Use a finite source with retained
order at least two. Path values are finite jets, and differentiation
reduces the certified order by one; they do not by themselves evaluate
a finite-time continuum evolution.

## Transport and phase

The actual source lift is V(t)K(t,r)V(r)*, with instantaneous Hamiltonian

    H_control=dot(beta).X+i h V B_theta V*.

The Schrödinger equation is proved on Dom X_1 intersect Dom X_2.
For each original current A, the uncut lift preserves Dom A intersect
Dom T. This does not assert preservation of every full current domain
without total-mass regularity.

For a rectangle +eta theta_i,+zeta beta_i,-eta theta_i,-zeta beta_i,

    gamma_B=(zeta/h)[m_i(theta+eta e_i)-m_i(theta)]<0

when eta,zeta>0. It is a mixed amplitude/phase Berry loop. The uncorrected
lift returns the base state with zero total phase because its dynamical
phase cancels this geometric phase. A coherent reference and dynamical-
phase calibration are needed to expose the geometric phase. This is
distinct from the theta-only TVSP response-plane holonomy.

## Material comparison and validation

The separate coherent family in PCR-T5 recovers smeared classical density
and conjugate material data as h tends to zero. In its one-channel positive
homogeneous specialization, beta=<D(1)>/<T> and r=Y^(-1)(beta). The length
exp(r) is an inferred classical family label, not a sharp quantum ratio or
pointwise quantum material-metric operator. The DCD curvature enclosure
has not been assigned to this separate coherent family.

```bash
python3.12 papers/native-phase-current-response/run.py --check
python3.12 scripts/replay_quantum_geometry.py
```

The packet replays 247 exact checks, 14 symbolic certificates and 35
negative controls. Python 3.12 handles I/O and source hashes; native JS
handles scientific arithmetic and proof replay. The integrated chain
contains 47 stages. Written analytic proofs supply the continuum claims.
