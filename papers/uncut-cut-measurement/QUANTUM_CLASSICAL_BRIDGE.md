# Quantum and classical readings of one native graph

Status: **conditional finite bridge**, with written proofs and exact rational source-bound controls. The paired carrier, positive form, Cayley continuation and source cost come from U27/U36 **under their declared graph and cost inputs**. A represented quarter-turn is not automatically the central primitive `iota`; an operational detector, frequency interpretation, physical clock and `hbar` calibration are separate.

## 1. Common native system and classical equilibrium

Take one connected grounded graph with positive stiffness `H` on `m` retained labels. On the paired carrier `W=R^m tensor R^2`, set `R=[[0,-1],[1,0]]`, `Z=I_m tensor R`, `A=H tensor R`, `K=H tensor I_2`, and declare a two-component source `b`. U36 constructs the affine Cayley continuation

`x_next=x_*+T_h(x-x_*)`, `x_*=K^{-1}b`, `T_h=(I+hA)(I-hA)^{-1}`.

It preserves both Euclidean norm of the centered fluctuation `xi=x-x_*` and its source excess cost `xi^T K xi/2`. The **stationary source reading** is the classical field `x_*`; U36's graph flux law holds at that stationary field. If two real aligned field readings are independently calibrated as entropy and volume, U36/TF-2 also produce an equilibrium thermodynamic `U(S,V)` with Maxwell identities. Stationarity or equilibration in nature is *not* forced by conservative Cayley iteration.

### QB-1 (paired phase evolution)

Identify each real pair `(u_j,v_j)` of `xi` with the complex coordinate `psi_j=u_j+i v_j`. Then `Z` acts as multiplication by `i`, `A` as `iH`, and exactly

`psi_next=(I+i h H)(I-i h H)^{-1} psi =: U_h psi`.

`U_h` is unitary, so `||psi_next||^2=||psi||^2`; `H` is the **same** supplied grounded graph stiffness appearing in the stationary source problem. The native parameter `h` labels steps, not laboratory time. Writing an energy-time Schrödinger equation requires an independently calibrated action unit and a time adapter. This finite statement does not use the withdrawn scalar energy dynamics in Information Invariance v2.

**Proof.** `A^T=-A`, `[A,Z]=0`, and `I-hA` is invertible since `||(I-hA)v||²=||v||²+h²||Av||²`. Real Cayley orthogonality gives norm preservation. Complex identification turns the two real factors into the displayed rational complex matrix. □

### QB-2 (native norm share and the detector assumption)

For a **declared** orthogonal detector split into sites `j`, the natural nonnegative reading is `w_j=u_j²+v_j²`; on `xi!=0` set

`p_j=w_j/sum_k w_k`, `rho=psi psi†/(psi†psi)`.

Then `p_j>=0`, `sum_j p_j=1`, and `p_j=rho_jj`. Within the **declared class of homogeneous quadratic local readings**, invariance under each site's represented quarter-turn and equality under carrier permutations force every `w_j` to be the same positive constant times `u_j²+v_j²`, up to common normalization. An orthogonal change of detector basis gives squared amplitudes in the chosen basis. This is a mathematical Born-*form* rule derived from the positive native metric plus **detector covariance and additive quadratic-readout assumptions**. Identification of `p_j` with experimental event frequencies, and selection of an actual detector basis, remain physical adapters.

**Proof.** Any real quadratic local weight has a symmetric `2x2` coefficient matrix `B_j`. Invariance under `R` implies `R^T B_j R=B_j`, hence `B_j=c_j I_2`; permutation equivalence fixes `c_j=c>0`. Orthogonality gives `||xi||²=sum_j (u_j²+v_j²)`, proving positivity and normalization. `rho` is positive rank one and trace one by construction. □

### QB-3 (classical shadow and exact seam-memory residue)

Let `D(rho)=diag(rho)` be the **declared phase-erasing cut** in the site detector basis. For a dephased state with `p_j=rho_jj`, one native step followed by this cut gives the classical probability update

`p'_k=sum_j P_kj p_j`, `P_kj=|(U_h)_kj|²`.

`P` is doubly stochastic. After *each* step's phase erasure, repeated readings follow `p^(n+1)=P p^n`. Without that intervention, generally `diag(U_h^n rho U_h†^n) != P^n diag(rho)`: the cut loses phase information needed by later dynamics. Precisely, with `M=rho-D(rho)` and `L(rho)=U_h rho U_h†`,

`D(L(rho))-D(L(D(rho))) = D(L(M))`.

This identity locates the quantum/classical seam in one common model. Classical dynamics is the **repeatedly phase-erased protocol**; it is not the full coherent trajectory under a relabelled name. Reconstructing arbitrary coherent futures from a diagonal alone is impossible when the right-hand side is nonzero. A cut-history or off-diagonal memory variable is required.

**Proof.** Expand the diagonal of `U_h diag(p) U_h†`; unitarity makes every row and column sum of `P` equal one. The residue follows by writing `rho=D(rho)+M` and using linearity. □

## 2. One graph, numerical separation in exact fractions

Use the connected two-site graph with one unit link and unit grounding at each site:

`H=[[2,-1],[-1,2]]`, `det H=3`, `b=(1,2)` on the real quadratures.

Its stationary field is `x_*=(4/3,5/3)` in the real quadrature and zero in the conjugate one. Its real equilibrium cost, with visible coordinates `(S,V)`, is

`U(S,V)=S²-SV+V²-S-2V`.

On the chart point `(S,V)=(2,1)`, the calibrated equilibrium quantities are `T=2`, `P=2`, `C_V=1`, `C_P=4/3`, `K_S=2`, `K_T=3/2`, giving `C_P/C_V=K_S/K_T=4/3`. The entropy/volume/energy interpretation and measured units are still external identifications.

For the **same graph** and a centered fluctuation prepared as `psi_0=(1,0)`, take the native parameter `h=1`. Exact Cayley arithmetic gives

`U_1=[[(-2+4i)/5,(2+i)/5],[(2+i)/5,(-2+4i)/5]]`,

`P=[[4/5,1/5],[1/5,4/5]]`.

After two coherent steps, site 2 has weight **`16/25`**. Erase phase after step 1 and continue under the *same* `U_1`: site 2 has weight **`8/25`**. The difference **`8/25`** is the off-diagonal seam-memory contribution, not numerical noise or an unrelated classical model. Starting from a site basis state, both protocols agree at step 1. A different relative phase at fixed first-step populations can change the next coherent readout; the diagonal is therefore not a closed state of the uncut phase dynamics.

The [exact rational certificate](certificates/quantum_classical_bridge.py) checks `H` as a grounded graph, the affine fixed point, the same Hessian's thermo identities, Cayley defining equation and unitarity, Born-form normalization, stochasticity, two-step interference, residue equality, and a negative control for false diagonal closure.

## 3. Claim boundary and next physical test

QB-1–QB-3 construct **quantum-form and classical-form descriptions of the same finite native continuation**, with the cut's information loss quantified. The names “quantum” and “classical” here identify mathematical/operational forms; there is as yet no evidence that physical detectors follow the declared quadratic weights, that native steps equal laboratory time, or that this graph is selected by uncut ground. `hbar` is neither assumed in the finite step nor numerically derived. The stationary thermodynamic reading and dynamic quantum-form reading share `H`, but stationary preparation does not follow from the unitary trajectory.

The next decisive physical task is a **predeclared detector-and-clock adapter**: identify a reproducible preparation, detector projectors and readout frequencies, calibrate the step parameter and energy unit, and compare the coherent `16/25`-type and deliberately phase-erased `8/25`-type predictions within uncertainty on a chosen physical platform. A match would validate that specific adapter; it would not by itself derive a universal spacetime metric or quantum gravity.

Sources: [U27's interacting paired carrier](NATIVE_INTERACTION_AND_LAWFUL_CUTS.md), [U36's common graph, source cost and affine Cayley law](SOURCE_RESPONSE_AND_GRAVITY_TESTS.md), [thermodynamic foundation TF-1/TF-2](THERMODYNAMIC_FOUNDATION_BRIDGE.md), and [Information Invariance's conditional complex-flow criterion and D2b verdict](https://github.com/Parveen117/Publications/blob/a3ff8e9f8e6afdea44c55557346e8b5fea0c5d22/papers/information-invariance/LEDGER.md). Squared-amplitude rules and phase erasure are familiar mathematical forms; no global priority is asserted for them.
