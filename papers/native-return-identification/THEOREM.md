# Identifying native interaction strength from a symmetric response probe

Monty Dabas · NI-1–NI-4 · 30 September 2026

This continues the [alpha selection audit](../native-alpha-selection/THEOREM.md).
It separates two questions: does an existing symmetry select the remaining
cell parameter b, and can a complete response observation identify b without
an absolute probe-gain calibration? For the canonical R2 family the first
candidate symmetries below do not select it; the second question has an exact
constructive answer. This is identification from response data, not a
parameter-free prediction of b or of alpha.

The source carrier and solver are the canonical RKF R2 paired-depth EMK
family at commit `3cc5a33b05c16d59c90994ddda69dedc0d392424`. Its presentation
has R²=-1, K²=1, KR=-RK and L=KR. At the exact cut, cells alternate between
a=b+1 and b>0, and the completed response is I+L. A supplied common cell
multiplier z>0 gives the unique positive return x(z) satisfying

\[
bz x^2+x-(b+1)z=0,\qquad x(1)=1. \tag{1}
\]

These are the original source hypotheses, not a photon identification.
The complete native aperture memory is retained through the source solver.

## NI-1. Two plausible symmetry constraints do not select a finite b

**One-cell translation equality.** Requiring the two cell products to be
identical gives a=b, incompatible with exact cut closure a=b+1. Equivalently,
requiring both phase responses x=y=1 contradicts the source recurrence
y=b/(1+bx): at x=1 it gives y=b/(1+b)<1 for every finite positive b.
Thus exact closure at every one-cell-shifted aperture is impossible in this
finite-coupling periodic class. This does not rule out other carriers, finite
approximation, or a controlled infinite-coupling limit.

**Balanced directional scalar amplitudes.** Each cell contains two matched
edge products. Requiring each individual opening scalar to equal its closing
scalar still realizes every positive product a_j by choosing both equal to
sqrt(a_j) in the completed positive radial field. It leaves b arbitrary.
This scalar balance is not a theorem of physical reciprocity or a statement
that the R and K edge arrows are daggers of one another. In fact their native
dagger rules must not be erased by equating their scalar magnitudes.

The source cut-square theorem factors supplied positive weights; it does not
select these edge products by minimizing an unprovided physical action.
The canonical thermo response-cost theorem likewise starts with a declared
source reconstruction and stable response metric. Neither statement supplies
a value of b for this excursion family. This is a scoped dependency audit,
not an assertion about every possible future native source law.

## NI-2. A differential shape observable removes affine probe gain

AS-1 gives s=x'(1)=1/(1+2b) and
x''(1)=-2b(3+4b)/(1+2b)^3. Therefore define

\[
\mathcal I=-\frac{x(1)x''(1)}{x'(1)^2}
=\frac{2b(3+4b)}{1+2b}=\frac2s-1-s. \tag{2}
\]

A multiplicative readout y(t)=A x(1+g t), A>0 and g≠0, has the same
I=-y(0)y''(0)/y'(0)^2. Thus I is unchanged by either readout amplitude or
an unknown affine probe gain. Additive backgrounds must be subtracted; a
nonlinear probe z(t) contributes a z'' term and is not covered.

On b>0 the function I is strictly increasing from zero to infinity. This
follows from s decreasing from one to zero and d(2/s-1-s)/ds=-2/s²-1<0.
Solving gives the unique inverse

\[
\boxed{b=\frac{\mathcal I-3+\sqrt{\mathcal I^2+2\mathcal I+9}}8.} \tag{3}
\]

It identifies b once that shape observable is supplied. It does not make
I or b into alpha. The derivatives and input normalization must refer to the
specified paired-cell probe, not an arbitrary laboratory control variable.

## NI-3. Exact two-probe reconstruction without numerical derivatives

A stronger finite experiment uses normalized readings
u=x(1+δ), v=x(1-δ), where the same unknown δ obeys 0<δ<1. The reference
reading is x(1)=1. AS-2 proves strict increase and strict concavity, so

\[
u>1,\quad 0<v<1,\quad u+v<2. \tag{4}
\]

Conversely every pair satisfying (4) is realized by exactly one b>0 and one
δ∈(0,1) in this family. To prove and construct this, set A=1-u²<0 and
B=1-v²>0. Rearranging (1) gives

\[
z(x)=\frac{x}{1+b(1-x^2)}.
\]

The symmetry condition z(u)+z(v)=2 is equivalent, within positive
denominators, to

\[
\boxed{2AB b^2+[2(A+B)-uB-vA]b+2-u-v=0.} \tag{5}
\]

Its leading coefficient is negative and its constant term is positive, so
it has exactly one positive root (and one negative root). At b=0 its value
is positive. At b*=1/(u²-1) it is -u(1+b*B)<0. Hence the positive root
lies strictly below b*, making both original denominators positive; no
spurious root from denominator clearing is admitted. Define

\[
\boxed{\delta=\frac{u}{1+b(1-u^2)}-1
=1-\frac{v}{1+b(1-v^2)}.} \tag{6}
\]

The first expression is positive because its denominator is less than one;
the second is less than one because v and its denominator are positive.
Both original equations hold. Their unique positive solutions are u and v,
by the canonical return theorem. This proves existence and uniqueness of
both reconstructed quantities. No absolute probe-gain value was inserted.

**Exact witness.** Input only u=6/5 and v=2/5. Equation (5) reconstructs
b=5/7, a=12/7 and δ=3/4. Thus the probes were at z=7/4 and z=1/4.
The native solver independently encloses the two exact returns. A second
example u=8/7,v=11/14 gives b=7/15 and δ=1/3. These are chosen arithmetic
fixtures, not measured experimental data or predicted constants.

The protocol assumes the actual edge-product multipliers are symmetric
about z=1. Symmetry of arbitrary instrument dial positions alone does not
establish that assumption. It also assumes the measured return coefficient
has its additive offset removed and is normalized by the reference cut.
A multiplicative readout gain cancels in these ratios.

## NI-4. Retained-tail and observation-error intervals

No floating-point root or fitted optimiser is required. Bisection of the
quadratic (5) on [0,1/(u²-1)] gives exact rational brackets for b. The signs
at the endpoints and the unique positive root prove every bracket. The
algorithm stops only when its declared width tolerance is met; a budget
exhaustion is reported, not promoted to an exact root.

For uncertain normalized readings u∈[uL,uU], v∈[vL,vU], require the entire
box to satisfy (4): uL>1, 0<vL≤vU<1 and uU+vU<2. The reconstructed b is
strictly decreasing in each argument. Here is a sign proof: put
G(b,u,v)=z(u)+z(v)-2. At the root, Q=-D G, where Q is (5) and
D=(1+bA)(1+bB)>0. The derivative Q_b at its positive root is negative,
so G_b>0. Also z_x=(1+b+bx²)/(1+b(1-x²))²>0. Implicit differentiation
then gives b_u=-G_u/G_b<0 and b_v=-G_v/G_b<0. Consequently

\[
\boxed{b(u_U,v_U)\le b\le b(u_L,v_L).} \tag{7}
\]

Compute outward rational root brackets at the two corners to get a certified
parameter interval. R2's native finite-aperture intervals are legitimate
inputs when the positive periodic and nonnegative-tail source contracts
hold. Increasing aperture depth reduces that error. An independent instrument
error would additionally need its own warranted intervals and correlations;
a rectangle is a conservative enclosure when such bounds are available.

The code obtains u,v intervals from the unchanged native solver, then uses
(7) to enclose the independently known fixture b. It rejects data boxes
crossing the admitted domain. As u+v approaches two the inference approaches
a boundary and relative conditioning can deteriorate; exact uniqueness alone
does not guarantee an experimentally precise identification.

## What this closes for the constants programme

One source-calibration obstacle is removed: within this return family, two
symmetric response probes determine the cell parameter and unknown probe
amplitude together. The independent source law that predicts these responses
from primitive data is still missing. Inputting measured u,v estimates b;
it is not a parameter-free derivation.

The electromagnetic bridge must additionally identify the charged source,
photon pole and action normalization from AS-3–AS-5. No electromagnetic
observable or measured alpha entered this packet. **Alpha remains NOT DERIVED.**
This is a concrete response-identification experiment and a theorem about
its mathematical admissibility; it is not a report that the experiment has
been carried out. No private philosophical text is reproduced.
