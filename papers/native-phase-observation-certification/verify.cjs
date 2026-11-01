'use strict';
// Certificate construction uses only the unchanged, source-pinned RKF engines.
// JS supplies traversal/control; native F, Cut, Expr and matrix operations do math.
const fs = require('node:fs'), path = require('node:path'), crypto = require('node:crypto');
const assert = require('node:assert/strict');
const root = __dirname, args = process.argv.slice(2);
const value = flag => { const i = args.indexOf(flag); assert.ok(i >= 0 && args[i + 1], flag); return args[i + 1]; };
const nativeRoot = path.resolve(value('--native-root'));
const upstream = JSON.parse(fs.readFileSync(value('--upstream-report'), 'utf8'));
const write = args.includes('--write');
assert.ok(write !== args.includes('--check'), 'Choose exactly one of --write and --check');
const pins = JSON.parse(fs.readFileSync(path.join(root, 'SOURCE_PINS.json'), 'utf8'));
const sha256 = data => crypto.createHash('sha256').update(data).digest('hex');
for (const [file, expected] of Object.entries(pins.native_engine.files))
  assert.equal(sha256(fs.readFileSync(path.join(nativeRoot, file))), expected, file + ': native source pin');
const N = require(path.join(nativeRoot, 'native_operator.cjs'));
const M = require(path.join(nativeRoot, 'paninian_operator.cjs'));
const {F, Cut, IOTA} = N, f = x => F.of(x);
const checks = [], controls = [], certificates = [], presentations = [], derivations = [];
const check = (name, ok, evidence = {}) => { assert.ok(ok, name); checks.push({name, passed: true, ...evidence}); };
const control = (name, ok, evidence = {}) => { assert.ok(ok, name); controls.push({name, false_inference_rejected: true, ...evidence}); };
const rejects = fn => { try { fn(); return false; } catch { return true; } };
const systems = new Map();
function presentation(name, tokens, extra = []) {
  const rules = [];
  for (let j = 1; j < tokens.length; j++) for (let i = 0; i < j; i++)
    rules.push({id: tokens[j] + '_past_' + tokens[i], lhs: [tokens[j], tokens[i]],
      rhs: [[[tokens[i], tokens[j]], 1]], source: 'Declared real central comparison coefficients'});
  const S = new M.Presentation({tokens, rules: rules.concat(extra)});
  const w = Object.fromEntries(tokens.map(t => [t, S.word([t])]));
  const audit = S.audit(), dagger = M.daggerGate(S, w);
  check(name + ':confluence', audit.status === 'CONFLUENT_BY_CHECKED_DIAMONDS');
  check(name + ':dagger', dagger.status === 'DAGGER_DESCENDS');
  presentations.push({name, spec: S.spec(), audit, dagger});
  systems.set(S.hash(), S);
  return {S, w};
}
function prove(name, lhs, rhs) {
  const S = lhs.system, result = M.proveEquality(lhs, rhs);
  check(name, result.status === 'EQUAL_IN_DECLARED_QUOTIENT');
  const certificate = JSON.parse(JSON.stringify(result.certificate));
  // Bind the witness to the advertised proposition, including the zero target.
  assert.ok(S.fromJSON(certificate.input).equals(lhs.minus(rhs)));
  assert.ok(S.fromJSON(certificate.output).equals(S.zero(lhs.from, lhs.to)));
  check(name + ':serialized_replay', S.replay(certificate) === true);
  certificates.push({name, lhs: lhs.toJSON(), rhs: rhs.toJSON(), certificate});
}
const sq = x => x.times(x), times = (...xs) => xs.reduce((a, b) => a.times(b));
const pow = (x, k) => { let y = x.system.one(); for (let j = 0; j < k; j++) y = y.times(x); return y; };
const derivative = (e, values) => M.polynomialDerivative(e, values);
const circle = (a, b) => ({id: b + '_unit_square', lhs: [b, b],
  rhs: [[[], 1], [[a, a], -1]], source: 'Declared native unit pair: ' + a + '^2+' + b + '^2=1'});
const I = N.identity(2), R = N.matrix([[0, -1], [1, 0]]), K = N.matrix([[0, 1], [1, 0]]);
const J = N.scale(K, IOTA), Z = N.scale(N.mul(R, K), -1), iR = N.scale(R, IOTA);
function symbolicMatrices(S) {
  const zero = S.zero(), one = S.one();
  const constant = A => A.map(row => row.map(x => one.scale(x)));
  const add = (A, B) => A.map((row, i) => row.map((x, j) => x.plus(B[i][j])));
  const sub = (A, B) => A.map((row, i) => row.map((x, j) => x.minus(B[i][j])));
  const scale = (A, x) => A.map(row => row.map(a => a.times(x)));
  const mul = (A, B) => A.map(row => B[0].map((_, j) => row.reduce((z, x, k) => z.plus(x.times(B[k][j])), zero)));
  const mm = (...xs) => xs.reduce(mul), trace = A => A[0][0].plus(A[1][1]);
  const eq = (name, A, B) => { for (let i = 0; i < 2; i++) for (let j = 0; j < 2; j++) prove(name + ':' + i + j, A[i][j], B[i][j]); };
  return {constant, add, sub, scale, mul, mm, trace, eq, id: constant(I), r: constant(R), j: constant(J)};
}

// PM observation contract: both oriented coefficients, not two signs alone.
{
  const {S, w: {a, b, c, d}} = presentation('PM_unit_pairs', ['a', 'b', 'c', 'd'], [circle('a', 'b'), circle('c', 'd')]);
  const s = symbolicMatrices(S), U = s.add(s.scale(s.id, a), s.scale(s.r, b));
  const Ud = s.sub(s.scale(s.id, a), s.scale(s.r, b));
  const V = s.add(s.scale(s.id, c), s.scale(s.j, d)), Vd = s.sub(s.scale(s.id, c), s.scale(s.j, d));
  const X = sq(a).minus(sq(b)), Y = a.times(b).scale(2), P = sq(c).minus(sq(d)), Q = c.times(d).scale(2);
  prove('PM:channel_X_Q', s.trace(s.mm(s.constant(K), Ud, s.constant(K), U)).scale('1/2'), X);
  prove('PM:channel_Y_Q', s.trace(s.mm(s.constant(Z), Ud, s.constant(K), U)).scale('1/2'), Y);
  prove('PM:channel_X_G', s.trace(s.mm(s.constant(iR), Vd, s.constant(iR), V)).scale('1/2'), P);
  prove('PM:channel_Y_G', s.trace(s.mm(s.constant(Z), Vd, s.constant(iR), V)).scale('1/2'), Q);
  prove('PM:doubled_cut_phase', sq(a.plus(b.scale(IOTA))), X.plus(Y.scale(IOTA)));
  prove('PM:unit_channel', sq(X).plus(sq(Y)), S.one());
  const C = X.times(P).plus(Y.times(Q)), T = Y.times(P).minus(X.times(Q));
  prove('PM:relative_unit_channel', sq(C).plus(sq(T)), S.one());
}
{
  const {S, w: {C, T}} = presentation('PM_relative_Cayley', ['C', 'T'], [circle('C', 'T')]);
  // r=T/(1+C). This division-free identity requires 1+C != 0 when divided.
  prove('PM1:relative_Cayley_cross_multiplication', C.plus(T.scale(IOTA)).times(S.one().plus(C).minus(T.scale(IOTA))),
    S.one().plus(C).plus(T.scale(IOTA)));
}
{
  const {w: {p0, p1, m0, k, pi}} = presentation('PM_record_update', ['p0', 'p1', 'm0', 'k', 'pi']);
  prove('PM2:record_update_identity', p0.plus(pi.times(m0)).plus(p1.minus(p0).plus(pi.times(k))), p1.plus(pi.times(m0.plus(k))));
  // pi here is a formal scalar; integer admission and inequalities are separate.
}

// AC1: general polynomial identities first, then an admitted tangent derivation.
{
  const {S, w: {a, b, da, db}} = presentation('AC_unrestricted_pair', ['a', 'b', 'da', 'db']);
  const X = sq(a).minus(sq(b)), Y = a.times(b).scale(2), radius = sq(a).plus(sq(b));
  const deriv = {a: da, b: db, da: S.zero(), db: S.zero()};
  const angular = a.times(db).minus(b.times(da)), radial = a.times(da).plus(b.times(db));
  prove('AC1:channel_norm_factorization', sq(X).plus(sq(Y)), sq(radius));
  prove('AC1:channel_current_numerator', X.times(derivative(Y, deriv)).minus(Y.times(derivative(X, deriv))), radius.times(angular).scale(2));
  prove('AC1:tangent_reconstruction_da', a.times(radial).minus(b.times(angular)), radius.times(da));
  prove('AC1:tangent_reconstruction_db', b.times(radial).plus(a.times(angular)), radius.times(db));
}
{
  const {S, w: {a, b, t}} = presentation('AC_unit_tangent', ['a', 'b', 't'], [circle('a', 'b')]);
  const deriv = {a: b.times(t).scale(-1), b: a.times(t), t: S.zero()};
  const gate = M.derivationGate(S, deriv);
  check('AC1:tangent_derivation_descends', gate.status === 'DERIVATION_DESCENDS');
  derivations.push({name: 'native_unit_tangent', presentation: S.hash(), images: Object.fromEntries(Object.entries(deriv).map(([k, v]) => [k, v.toJSON()])), gate});
  const X = sq(a).minus(sq(b)), Y = a.times(b).scale(2), s = symbolicMatrices(S);
  prove('AC1:unit_current_is_t', X.times(derivative(Y, deriv)).minus(Y.times(derivative(X, deriv))), t.scale(2));
  const Ud = s.sub(s.scale(s.id, a), s.scale(s.r, b)), dU = s.add(s.scale(s.id, deriv.a), s.scale(s.r, deriv.b));
  s.eq('AC1:operator_current', s.mul(Ud, dU), s.scale(s.r, t));
  prove('AC1:trace_current', s.trace(s.mm(s.r, Ud, dU)).scale('-1/2'), t);
  control('an_arbitrary_derivative_need_not_descend_to_the_unit_circle', M.derivationGate(S, {a: S.zero(), b: S.one(), t: S.zero()}).status === 'RELATION_NOT_PRESERVED');
}
{
  const {S, w: {q, dq}} = presentation('AC_Cayley_differential', ['q', 'dq']);
  const A = S.one().minus(sq(q)), B = q.scale(2), D = S.one().plus(sq(q)), deriv = {q: dq, dq: S.zero()};
  const dA = derivative(A, deriv), dB = derivative(B, deriv), dD = derivative(D, deriv);
  prove('AC1:Cayley_unit_identity', sq(A).plus(sq(B)), sq(D));
  prove('AC1:Cayley_current_cross_multiplication', A.times(dB.times(D).minus(B.times(dD))).minus(B.times(dA.times(D).minus(A.times(dD)))), dq.times(sq(D)).scale(2));
}
{
  const {S, w: {x, y, dx, dy, rho, drho}} = presentation('AC_common_gain', ['x', 'y', 'dx', 'dy', 'rho', 'drho']);
  const X = rho.times(x), Y = rho.times(y), dX = drho.times(x).plus(rho.times(dx)), dY = drho.times(y).plus(rho.times(dy));
  prove('AC1:varying_gain_numerator', X.times(dY).minus(Y.times(dX)), sq(rho).times(x.times(dy).minus(y.times(dx))));
  prove('AC1:varying_gain_denominator', sq(X).plus(sq(Y)), sq(rho).times(sq(x).plus(sq(y))));
  prove('AC3:real_gain_has_zero_current_numerator', rho.times(S.zero()).minus(S.zero().times(drho)), S.zero());
}
{
  const {S, w: {x, y, u, v, dx, dy, du, dv}} = presentation('AC_cut_product', ['x', 'y', 'u', 'v', 'dx', 'dy', 'du', 'dv']);
  const deriv = {x: dx, y: dy, u: du, v: dv, dx: S.zero(), dy: S.zero(), du: S.zero(), dv: S.zero()};
  const Nu = x.times(dy).minus(y.times(dx)), Nv = u.times(dv).minus(v.times(du));
  const Du = sq(x).plus(sq(y)), Dv = sq(u).plus(sq(v));
  const P = x.times(u).minus(y.times(v)), Q = x.times(v).plus(y.times(u));
  prove('AC3:product_current_numerator', P.times(derivative(Q, deriv)).minus(Q.times(derivative(P, deriv))), Nu.times(Dv).plus(Nv.times(Du)));
  prove('AC3:product_current_denominator', sq(P).plus(sq(Q)), Du.times(Dv));
  prove('AC3:dagger_current_numerator', x.times(dy.scale(-1)).minus(y.scale(-1).times(dx)), Nu.scale(-1));
  const Pr = x.times(u).plus(y.times(v)), Qr = y.times(u).minus(x.times(v));
  prove('AC3:relative_current_numerator', Pr.times(derivative(Qr, deriv)).minus(Qr.times(derivative(Pr, deriv))), Nu.times(Dv).minus(Nv.times(Du)));
  prove('AC3:relative_current_denominator', sq(Pr).plus(sq(Qr)), Du.times(Dv));
}

// AC4 is the exact mixed coefficient of the displayed frozen word.
{
  const {S, w: {ax, ay, bx, by, eps, eta}} = presentation('AC_frozen_loop', ['ax', 'ay', 'bx', 'by', 'eps', 'eta'], [
    {id: 'eps_squared_zero', lhs: ['eps', 'eps'], rhs: [], source: 'Declared first-order increment eps^2=0'},
    {id: 'eta_squared_zero', lhs: ['eta', 'eta'], rhs: [], source: 'Declared first-order increment eta^2=0'}]);
  const s = symbolicMatrices(S), Ax = s.add(s.scale(s.r, ax), s.scale(s.j, bx)), Ay = s.add(s.scale(s.r, ay), s.scale(s.j, by));
  const area = ax.times(by).minus(ay.times(bx)), bracket = s.sub(s.mul(Ax, Ay), s.mul(Ay, Ax));
  const coefficient = s.scale(s.mul(s.r, s.j), area.scale(2));
  s.eq('AC4:current_area_bracket', bracket, coefficient);
  const H = s.mm(s.sub(s.id, s.scale(Ay, eta)), s.sub(s.id, s.scale(Ax, eps)), s.add(s.id, s.scale(Ay, eta)), s.add(s.id, s.scale(Ax, eps)));
  s.eq('AC4:ordered_frozen_loop', H, s.sub(s.id, s.scale(coefficient, eps.times(eta))));
  const reverse = s.mm(s.sub(s.id, s.scale(Ax, eps)), s.sub(s.id, s.scale(Ay, eta)), s.add(s.id, s.scale(Ax, eps)), s.add(s.id, s.scale(Ay, eta)));
  s.eq('AC4:reversed_frozen_loop', reverse, s.add(s.id, s.scale(coefficient, eps.times(eta))));
  control('the_opposite_loop_sign_is_not_the_displayed_word', M.proveEquality(H[0][0], s.add(s.id, s.scale(coefficient, eps.times(eta)))[0][0]).status === 'DISTINCT_IN_DECLARED_QUOTIENT');
}

// Exact finite polynomial instances of the endpoint-jet obstruction.
// The general r argument and native circle periodicity remain written premises.
const jetCases = [];
{
  const {S, w: {x}} = presentation('AC_endpoint_polynomials', ['x']);
  const evaluate = (e, at) => M.matrixValue(e, {x: N.matrix([[at]])}, 1)[0][0];
  const integratePolynomial = e => new M.Expr(S, '*', '*', S.reduce(e).normal.entries().map(([word, c]) => [word.concat('x'), c.div(f(word.length).add(1))]));
  for (let r = 1; r <= 6; r++) {
    const density = S.reduce(pow(x, r).times(pow(S.one().minus(x), r))).normal;
    const primitive = integratePolynomial(density), norm = evaluate(primitive, 1).sub(evaluate(primitive, 0));
    check('jets:' + r + ':positive_normalizer', norm.turn.zero() && !norm.rad.le(0));
    const h = primitive.scale(norm.inv()), dh = derivative(h, {x: S.one()});
    prove('jets:' + r + ':derivative_factorization', dh, density.scale(norm.inv()));
    if (r === 1) prove('jets:1:displayed_cubic', h, sq(x).scale(3).minus(pow(x, 3).scale(2)));
    const endpoints = [];
    let current = h;
    for (let k = 0; k <= r; k++) {
      const at0 = evaluate(current, 0), at1 = evaluate(current, 1);
      check('jets:' + r + ':endpoint_derivative:' + k, at0.zero() && at1.eq(k === 0 ? 1 : 0));
      endpoints.push({order: k, at0, at1});
      current = derivative(current, {x: S.one()});
    }
    const currentPrimitive = integratePolynomial(dh), integral = evaluate(currentPrimitive, 1).sub(evaluate(currentPrimitive, 0));
    check('jets:' + r + ':integrated_current_in_pi_units', integral.eq(1));
    const {S: T, w: {e}} = presentation('AC_endpoint_jet_' + r, ['e'], [
      {id: 'jet_truncation', lhs: Array(r + 1).fill('e'), rhs: [], source: 'Finite jet through order ' + r + ': e^' + (r + 1) + '=0'}]);
    const substitute = at => h.entries().reduce((z, [word, c]) => z.plus(pow(at, word.length).scale(c)), T.zero());
    prove('jets:' + r + ':initial_phase_jet_zero', substitute(e), T.zero());
    prove('jets:' + r + ':final_phase_jet_one', substitute(T.one().plus(e)), T.one());
    control('jets:' + r + ':matching_endpoint_jets_do_not_determine_accumulation', integral.eq(1) && !integral.zero());
    jetCases.push({order: r, phase_polynomial_in_pi_units: h.toJSON(), normalizer: norm, endpoints, integral_in_pi_units: integral});
  }
}

// PM finite execution. Rational phase labels below mean multiples of native pi.
// No numerical trigonometry, floating pi or root approximation is introduced.
const cutPower = (z, n) => { assert.ok(Number.isSafeInteger(n)); if (n < 0) { z = z.dagger(); n = -n; } let y = Cut.of(1); for (let k = 0; k < n; k++) y = y.mul(z); return y; };
const quarterChannel = theta => { const n = f(theta).mul(4); assert.equal(n.d, 1n); return cutPower(IOTA, Number(n.n)); };
const canonicalQuarter = z => { const j = [0, 1, 2, 3].find(j => cutPower(IOTA, j).eq(z)); assert.notEqual(j, undefined); return f(j).div(4); };
function relativeChannel(before, after, boundAttested) {
  assert.equal(boundAttested, true, 'An independent actual-step bound is required');
  assert.ok(before.norm2().eq(1) && after.norm2().eq(1), 'Unit channel required');
  const z = after.mul(before.dagger()), den = f(1).add(z.rad);
  assert.ok(!den.zero(), 'The strict half-turn boundary has no admitted short step');
  const r = z.turn.div(den);
  assert.ok(new Cut(1, r).div(new Cut(1, r.neg())).eq(z));
  return {C: z.rad, S: z.turn, channel_cayley: r};
}
function chartStep(before, after, boundAttested) {
  assert.equal(boundAttested, true);
  [before, after].forEach(x => assert.ok(f(0).le(x) && !f(1).le(x), 'Canonical phase must lie in [0,1) pi units'));
  const candidates = [-1, 0, 1].map(k => ({k: f(k), delta: after.sub(before).add(k)})).filter(x => !f('1/2').le(x.delta.abs()));
  assert.equal(candidates.length, 1, 'The open short-step interval must have exactly one representative');
  return candidates[0];
}
function admitMemory(lift, phi) { const m = f(lift).sub(phi); assert.equal(m.d, 1n, 'Winding residual must be an exact integer'); return m; }
function decodeQuarterPath(channels, anchor, boundAttested) {
  let phi = canonicalQuarter(channels[0]), m = admitMemory(anchor, phi), lift = f(anchor);
  const rows = [{channel: channels[0], phi_pi: phi, m, lift_pi: lift}];
  for (let j = 1; j < channels.length; j++) {
    const observation = relativeChannel(channels[j - 1], channels[j], boundAttested), next = canonicalQuarter(channels[j]);
    const step = chartStep(phi, next, boundAttested); m = m.add(step.k); lift = lift.add(step.delta); phi = next;
    assert.ok(lift.eq(phi.add(m))); admitMemory(lift, phi);
    rows.push({channel: channels[j], phi_pi: phi, m, lift_pi: lift, delta_pi: step.delta, crossing: step.k, ...observation});
  }
  return rows;
}
const paths = [];
for (const [name, phases] of [
  ['draft_forward', Array.from({length: 9}, (_, j) => f(j).div(4))],
  ['reverse', Array.from({length: 9}, (_, j) => f(j).div(-4))],
  ['mixed_with_zero_steps', ['0', '1/4', '1/4', '0', '-1/4', '-1/2', '-1/4', '0', '1/4'].map(f)]
]) {
  const channels = phases.map(quarterChannel), rows = decodeQuarterPath(channels, phases[0], true);
  rows.forEach((row, j) => check('PM2:' + name + ':anchored_lift:' + j, row.lift_pi.eq(phases[j])));
  paths.push({name, phase_unit: 'native pi, retained symbolically', rows});
}
check('PM2:displayed_winding_record', paths[0].rows.every((row, j) => row.m.eq(['0', '0', '0', '0', '1', '1', '1', '1', '2'][j])));
check('PM2:all_three_crossing_corrections_exercised', new Set(paths.flatMap(p => p.rows.slice(1).map(row => row.crossing.toString()))).size === 3);
control('PM:missing_independent_bound_is_rejected', rejects(() => relativeChannel(Cut.of(1), IOTA, false)));
control('PM:half_turn_boundary_is_rejected', rejects(() => relativeChannel(Cut.of(1), Cut.of(-1), true)) && rejects(() => chartStep(f(0), f('1/2'), true)));
control('PM:nonunit_channel_is_rejected', rejects(() => relativeChannel(Cut.of(2), IOTA, true)));
control('PM:incompatible_anchor_is_rejected', rejects(() => decodeQuarterPath([Cut.of(1), IOTA], '1/4', true)));
control('PM:endpoint_does_not_supply_the_missing_bound', quarterChannel('1/4').eq(quarterChannel('5/4')));
control('AC2:noninteger_memory_is_rejected_without_rounding', rejects(() => admitMemory('21/10', 0)));
const fractionalAliases = [];
for (const n of [1, 2, 3, 5]) for (const k of [-2, -1, 1, 2]) {
  const exposure = f(n).div(12), aliasDelta = f(12).mul(k), channelTurns = exposure.mul(aliasDelta);
  check('PM:rational_exposure_alias:' + n + ':' + k, channelTurns.d === 1n && channelTurns.eq(f(n).mul(k)));
  fractionalAliases.push({exposure, alias_delta_in_pi: aliasDelta, complete_channel_turns: channelTurns});
}
const localPairs = [];
const transportPair = q => { q = f(q); const d = f(1).add(q.pow(2)); return new Cut(f(1).sub(q.pow(2)).div(d), q.mul(2).div(d)); };
for (const q0 of ['0', '1/4', '1/3', '1/2', '2/3', '3/4']) for (const q1 of ['0', '1/4', '1/3', '1/2', '2/3', '3/4']) {
  const u0 = transportPair(q0), u1 = transportPair(q1), local = relativeChannel(u0.mul(u0), u1.mul(u1), true);
  check('PM1:rational_channel_Cayley:' + q0 + ':' + q1, f(1).add(local.C).mul(local.channel_cayley).eq(local.S));
  localPairs.push({q0, q1, ...local});
}

// AC finite jets use the native cut product and its exact product rule.
const jetMul = (a, b) => ({v: a.v.mul(b.v), d: a.d.mul(b.v).add(a.v.mul(b.d))});
const jetDagger = a => ({v: a.v.dagger(), d: a.d.dagger()});
const angular = a => { assert.ok(!a.v.norm2().zero(), 'Zero channel has no angular current'); return a.v.dagger().mul(a.d).turn.div(a.v.norm2().mul(2)); };
function pairJet(q, dq) {
  q = f(q); dq = f(dq); const d = f(1).add(q.pow(2)), u = transportPair(q);
  const du = new Cut(q.mul(dq).mul(-4).div(d.pow(2)), f(1).sub(q.pow(2)).mul(dq).mul(2).div(d.pow(2)));
  return {v: u.mul(u), d: u.mul(du).mul(2)};
}
const currentCases = [];
for (const q of ['0', '1/3', '1', '2', '-1/2']) for (const p of ['-2', '0', '1/2']) for (const g of ['-1/3', '2/3']) {
  const source = pairJet(q, '3/2'), reference = pairJet(p, '-2/3'), frame = pairJet(g, '5/7');
  const gainS = {v: Cut.of('8/3'), d: Cut.of(-4)}, gainR = {v: Cut.of('7/5'), d: Cut.of(3)};
  const seenS = jetMul(gainS, jetMul(frame, source)), seenR = jetMul(gainR, jetMul(frame, reference));
  const expected = angular(source).sub(angular(reference)), observed = angular(seenS).sub(angular(seenR));
  const tag = q + ':' + p + ':' + g;
  check('AC1:Cayley_current:' + tag, angular(source).eq(f(3).div(f(1).add(f(q).pow(2)))));
  check('AC3:common_frame_and_varying_gains:' + tag, observed.eq(expected));
  check('AC3:relative_product:' + tag, angular(jetMul(seenS, jetDagger(seenR))).eq(expected));
  currentCases.push({q, p, frame_q: g, source, reference, frame, gainS, gainR, expected_relative_current: expected, observed_relative_current: observed});
}
check('AC1:transport_chart_pole_has_regular_channel_current', angular({v: Cut.of(1), d: new Cut(0, 14)}).eq(7), {unit_pair: ['-1', '0'], tangent: ['0', '-7']});
control('AC1:zero_channel_is_rejected', rejects(() => angular({v: Cut.of(0), d: Cut.of(1)})));
const baseJet = {v: new Cut(1, 1), d: new Cut(0, 2)};
control('AC1:unequal_coordinate_gains_do_not_cancel', !angular(baseJet).eq(angular({v: new Cut(2, 3), d: new Cut(0, 6)})));
control('AC1:unknown_offsets_do_not_cancel', !angular(baseJet).eq(angular({v: new Cut(2, 1), d: baseJet.d})));
const sJet = pairJet('1/3', 1), rJet = pairJet('1/2', 1), rotating = pairJet('2/3', 1);
control('AC3:noncommon_frame_rotation_does_not_cancel', !angular(jetMul(rotating, sJet)).sub(angular(rJet)).eq(angular(sJet).sub(angular(rJet))));
control('AC3:unknown_reference_current_cannot_be_set_to_zero', !angular(sJet).sub(angular(rJet)).eq(angular(sJet)));
const nativeU = N.cayley(R, '1/3'), phaseU = N.scale(nativeU, new Cut('3/5', '4/5'));
control('PM_AC:extra_global_phase_is_channel_invisible', N.equal(N.mul(N.mul(N.dagger(nativeU), K), nativeU), N.mul(N.mul(N.dagger(phaseU), K), phaseU)) && !N.equal(nativeU, phaseU));
check('AC4:dependent_directions_have_zero_commutator', N.isZero(N.commutator(N.add(R, N.scale(J, 2)), N.add(N.scale(R, 3), N.scale(J, 6)))));
control('AC4:frozen_commutator_is_not_the_full_variable_field_curvature', N.isZero(N.commutator(N.zeros(2), N.scale(R, 2))) && !N.isZero(R),
  {field: 'A_x=0, A_y=x R in a commuting coordinate frame', frozen_commutator: '0', derivative_curvature: 'partial_x A_y=R'});
const circulation = [];
for (const nu of [-3, -2, -1, 0, 1, 2, 3]) {
  const sign = nu % 2 === 0 ? 1 : -1, unitPair = Cut.of(sign);
  check('AC2:closed_channel_parity:' + nu, unitPair.mul(unitPair).eq(1) && unitPair.eq(nu % 2 === 0 ? 1 : -1));
  circulation.push({integral_in_pi_units: f(nu), operator_endpoint_sign: f(sign), channel_endpoint: unitPair.mul(unitPair)});
}
check('AC:displayed_three_turn_memory', admitMemory(3, 0).eq(3));

// Adversarial proof controls. These are rejection checks, not new theorems.
{
  const item = certificates.find(x => x.certificate.steps.length > 0), S = systems.get(item.certificate.presentation);
  const badStep = JSON.parse(JSON.stringify(item.certificate)); badStep.steps[0].coefficient = ['991', '0'];
  control('certificate:changed_rewrite_coefficient_is_rejected', rejects(() => S.replay(badStep)));
  const badSource = JSON.parse(JSON.stringify(item.certificate)); badSource.presentation = '0'.repeat(64);
  control('certificate:wrong_presentation_hash_is_rejected', rejects(() => S.replay(badSource)));
  const badOutput = JSON.parse(JSON.stringify(item.certificate)); badOutput.output = S.one().toJSON();
  control('certificate:nonzero_claimed_normal_form_is_rejected', rejects(() => S.replay(badOutput)));
  control('native:floating_point_input_is_rejected', rejects(() => f(0.25)));
}

const output = {
  schema: 'native-phase-observation-certification.v1', status: 'PASS_SCOPED_NATIVE_CERTIFICATION',
  source_commit: pins.source_commit,
  scope: 'Unchanged native engines: algebraic rewrite certificates, serialized replay, exact finite controls, and six inherited packet reruns. Analytic/path-lifting claims retain the written hypotheses in CLAIM_COVERAGE.md; no general-purpose formal proof or physical gravity validation is claimed.',
  source_pins_sha256: sha256(fs.readFileSync(path.join(root, 'SOURCE_PINS.json'))),
  native_engine: pins.native_engine, upstream,
  presentations, derivations, certificates, checks, negative_controls: controls,
  examples: {phase_paths: paths, local_rational_channel_pairs: localPairs, fractional_aliases: fractionalAliases,
    current_cases: currentCases, endpoint_jet_obstructions: jetCases, circulation_parity: circulation},
  summary: {native_presentations: presentations.length, symbolic_certificates: certificates.length,
    serialized_replays: certificates.length, exact_checks: checks.length, negative_controls: controls.length,
    inherited_packets: upstream.packets.length, phase_paths: paths.length, rational_channel_pairs: localPairs.length,
    moving_gain_frame_cases: currentCases.length, endpoint_jet_orders: jetCases.length}
};
const serialized = JSON.stringify(output, null, 2) + '\n', digest = sha256(serialized);
const certificatePath = path.join(root, 'CERTIFICATE.json'), expectedPath = path.join(root, 'EXPECTED.sha256');
if (write) {
  fs.writeFileSync(certificatePath, serialized);
  fs.writeFileSync(expectedPath, digest + '  CERTIFICATE.json\n');
} else {
  assert.equal(fs.readFileSync(certificatePath, 'utf8'), serialized, 'Frozen certificate differs from recomputation');
  assert.equal(fs.readFileSync(expectedPath, 'utf8'), digest + '  CERTIFICATE.json\n', 'Frozen certificate checksum differs');
}
console.log(JSON.stringify({status: output.status, mode: write ? 'explicit certificate generation' : 'read-only replay', ...output.summary, certificate_sha256: digest}, null, 2));
