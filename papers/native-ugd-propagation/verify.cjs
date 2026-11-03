'use strict';
/** Publication witness importing the canonical engine; no second arithmetic engine. */
const fs = require('node:fs'), path = require('node:path'), crypto = require('node:crypto');
const cp = require('node:child_process'), assert = require('node:assert/strict');
const here = __dirname, repo = path.resolve(here, '../..');
const argv = process.argv.slice(2), at = argv.indexOf('--rkf-root');
if (at < 0 || !argv[at+1] || (argv.includes('--write') && argv.includes('--check'))) {
  throw new Error('Usage: node verify.cjs --rkf-root /path/to/Recognition-Kernel-Framework [--check|--write]');
}
const rkf = path.resolve(argv[at+1]), engine = path.join(rkf, 'operator_foundation');
const sha = x => crypto.createHash('sha256').update(x).digest('hex');
const pins = JSON.parse(fs.readFileSync(path.join(here, 'SOURCE_PINS.json'), 'utf8'));
for (const [root, files] of [[rkf, pins.rkf_files], [repo, pins.publications_files]]) {
  for (const f of files) assert.equal(sha(fs.readFileSync(path.join(root, f.path))), f.sha256,
                                    'Consumed source drift: '+f.path);
}
const o = require(path.join(engine, 'core/native_operator.cjs'));
const p = require(path.join(engine, 'core/paninian_operator.cjs'));
const wb = require(path.join(engine, 'core/workbench.cjs'));
const paths = require(path.join(engine, 'core/native_paths.cjs'));
const gate = require(path.join(engine, 'core/dependency_gate.cjs'));
const {WeightedCompletion} = require(path.join(engine, 'core/weighted_completion.cjs'));
const {F, Cut} = o;
const checks = [], proofs = {};
function check(name, f) { f(); checks.push(name); }
function equal(s, a, b) {
  const r = p.proveEquality(a,b);
  assert.equal(r.status,'EQUAL_IN_DECLARED_QUOTIENT');
  assert(s.replay(r.certificate));
  return r.certificate;
}
const base = JSON.parse(fs.readFileSync(path.join(engine,'examples/emk_job.json'),'utf8'));
const model = wb.buildModel(base.presentation);
check('canonical_engine_derives_EMK_basis', () => {
  assert.equal(model.basis.status,'COMPLETE_FINITE_NORMAL_BASIS');
  assert.equal(model.basis.dimension,4);
});
const s = model.presentation, one = s.one(), R = s.word(['R']), K = s.word(['K']);
const P = one.plus(K).scale(new Cut('1/2')), Q = one.minus(P);
check('native_cut_before_representation', () => {
  equal(s,P.times(P),P); equal(s,Q.times(Q),Q); equal(s,P.times(Q),s.zero());
  assert.equal(p.daggerGate(s,{R:R.scale(-1),K}).status,'DAGGER_DESCENDS');
});
check('recognition_compression_retains_returning_memory', () => {
  proofs.invisible_step = equal(s,P.times(R).times(P),s.zero());
  proofs.visible_return = equal(s,P.times(R).times(R).times(P),P.scale(-1));
  proofs.memory_corner = equal(s,P.times(R).times(Q).times(R).times(P),P.scale(-1));
  assert.equal(p.proveEquality(P.times(R).times(R).times(P),s.zero()).status,
               'DISTINCT_IN_DECLARED_QUOTIENT');
});
const coefficients = e => model.coefficients(e);
const stateBasis = [P,R.times(P)], columns = stateBasis.map(coefficients);
const lift = o.matrix(Array.from({length:4},(_,i)=>columns.map(v=>v[i])));
function descendedOnLeftIdeal(a) {
  return o.matrix([0,1].map(i => stateBasis.map(b => {
    const sol = p.solveCoefficients(a.times(b),stateBasis);
    assert.equal(sol.status,'UNIQUE_SOLUTION'); return sol.particular[i];
  })));
}
const actionR = descendedOnLeftIdeal(R), actionK = descendedOnLeftIdeal(K);
const future = o.rowClosure([[1,0]],[actionR,actionK]);
check('state_action_is_derived_and_intertwined', () => {
  assert(o.equal(o.mul(model.leftAction(R),lift),o.mul(lift,actionR)));
  assert(o.equal(o.mul(model.leftAction(K),lift),o.mul(lift,actionK)));
  assert.equal(future.rank,2); assert.equal(future.extra,1);
});
check('native_derivation_must_preserve_relations', () => {
  assert.equal(p.derivationGate(s,{R:p.bracket(K,R),K:s.zero()}).status,'DERIVATION_DESCENDS');
  assert.equal(p.derivationGate(s,{R:one,K:s.zero()}).status,'RELATION_NOT_PRESERVED');
});
// This Laurent sheet extension is not the proper one-sided seam algebra.
const spec = JSON.parse(JSON.stringify(base.presentation));
spec.tokens = ['R','K',{name:'Z',from:'*',to:'*',charge:1},
                        {name:'Zi',from:'*',to:'*',charge:-1}];
for (const [id,lhs,rhs] of [
  ['ZZi',['Z','Zi'],[]], ['ZiZ',['Zi','Z'],[]],
  ['ZR',['Z','R'],['R','Z']], ['ZK',['Z','K'],['K','Z']],
  ['ZiR',['Zi','R'],['R','Zi']], ['ZiK',['Zi','K'],['K','Zi']]]) {
  spec.rules.push({id,lhs,rhs:[[rhs,1]],source:'NC2: declared product of EMK coefficients and UGD integer translations'});
}
const weighted = new WeightedCompletion(spec,{R:1,K:1,Z:1,Zi:1});
const h = weighted.presentation, hr = h.word(['R']), hk = h.word(['K']), z = h.word(['Z']);
const hp = h.one().plus(hk).scale(new Cut('1/2')), hq = h.one().minus(hp);
const U = hr.times(z), pow = n => weighted.power(U,n);
check('full_sheet_carrier_is_infinite_and_charge_preserving', () => {
  assert.equal(wb.finiteBasis(h).status,'INFINITE_IRREDUCIBLE_LANGUAGE');
  assert.throws(() => h.addRule({id:'erase_sheet',lhs:['Z','Z','Z','Z'],rhs:[[[],1]],
                                source:'planted illegal sheet erasure'}),/integer residue/);
});
check('native_dagger_on_full_sheet_algebra', () => {
  const star = {R:hr.scale(-1),K:hk,Z:h.word(['Zi']),Zi:z};
  assert.equal(weighted.certifyDagger(star).status,'DAGGER_EXTENDS_ISOMETRICALLY');
  equal(h,p.star(U,star).times(U),h.one());
  equal(h,U.times(p.star(U,star)),h.one());
});
const recurrence = [];
check('sheet_resolved_return_and_recurrence', () => {
  equal(h,hp.times(U).times(hp),h.zero());
  proofs.sheet_return = equal(h,hp.times(pow(2)).times(hp),hp.times(z).times(z).scale(-1));
  proofs.four_steps = equal(h,pow(4),weighted.power(z,4));
  assert.equal(p.proveEquality(pow(4),h.one()).status,'DISTINCT_IN_DECLARED_QUOTIENT');
  for(let n=0;n<=8;n++) {
    equal(h,hp.times(pow(n+2)),z.times(z).times(hp).times(pow(n)).scale(-1));
    recurrence.push(weighted.normal(hp.times(pow(n)).times(hp)).toJSON());
  }
});
check('clock_free_cut_differential_and_composition', () => {
  const differential = a => hp.times(a).minus(a.times(hp));
  const a=U,b=pow(3);
  equal(h,differential(b.times(a)),differential(b).times(a).plus(b.times(differential(a))));
  equal(h,hp.times(b).times(a).times(hp),hp.times(b).times(hp).times(a).times(hp)
           .plus(hp.times(b).times(hq).times(a).times(hp)));
});
check('Morphic_clock_second_derivative_retains_clock_curvature', () => {
  // a=tau', b=tau'' at a regular point. These are admitted clock data.
  for (const [a0,b0] of [[1,0],[2,0],[2,2]]) {
    const a=new F(a0),b=new F(b0),a2=a.pow(2),a3=a.pow(3);
    const first=hp.times(U).scale(new Cut(new F(1).div(a)));
    const second=hp.times(pow(2)).scale(new Cut(new F(1).div(a2)))
      .minus(hp.times(U).scale(new Cut(b.div(a3))));
    const restoring=z.times(z).times(hp).scale(new Cut(new F(1).div(a2)));
    equal(h,second.plus(first.scale(new Cut(b.div(a2)))).plus(restoring),h.zero());
    if(b0) assert.equal(p.proveEquality(second.plus(restoring),h.zero()).status,
                       'DISTINCT_IN_DECLARED_QUOTIENT');
  }
});
const cutIota = o.IOTA;
const pathStep = paths.arrow('cut','cut',1,cutIota);
let pathFour=paths.unit(['cut']);
for(let n=0;n<4;n++) pathFour=pathFour.mul(pathStep);
check('independent_native_path_backend_keeps_the_fourth_sheet', () => {
  assert(pathFour.eq(paths.arrow('cut','cut',4)));
  assert(!pathFour.eq(paths.unit(['cut'])));
  assert(o.equal(pathFour.endpointShadow(['cut']),o.identity(1)));
});
// The scalar quarter-turn is a checked coefficient-sector image, not the whole EMK algebra.
const properSpec = {tokens:[{name:'S',from:'*',to:'*',charge:1},
                          {name:'T',from:'*',to:'*',charge:-1}],
  rules:[{id:'TS',lhs:['T','S'],rhs:[[[],1]],source:'RKF R1 proper depth seam: TS=1 only'}]};
const proper = new p.Presentation(properSpec), boundary=proper.one().minus(proper.word(['S','T']));
check('proper_depth_seam_is_not_replaced_by_bilateral_sheet', () => {
  equal(proper,boundary.times(boundary),boundary);
  assert.equal(p.proveEquality(boundary,proper.zero()).status,'DISTINCT_IN_DECLARED_QUOTIENT');
  assert.equal(wb.finiteBasis(proper).status,'INFINITE_IRREDUCIBLE_LANGUAGE');
});
const series = p.formalExp(U,12).map(e=>weighted.normal(hp.times(e).times(hp)));
const t = new F(1,2);
function sumRange(lo,hi) {
  let out=h.zero();
  for(let j=lo;j<=hi;j++)out=out.plus(series[j].scale(new Cut(t.pow(j))));
  return weighted.normal(out);
}
const M2=sumRange(0,2),correction=sumRange(3,4),M4=sumRange(0,4);
check('Madhava_correction_transfers_exactly_from_Smriti', () => {
  equal(h,M4,M2.plus(correction));
  const tail2=sumRange(3,12),tail4=sumRange(5,12);
  equal(h,tail2,correction.plus(tail4));
  assert(weighted.mass(correction).eq(new F(1,384)));
  assert(!weighted.mass(tail4).zero());
  // This coefficient check is only modulo t^13; the infinite tail bound is below.
  assert.equal(p.proveEquality(M2.plus(correction).plus(tail2),M2.plus(tail2)).status,
               'DISTINCT_IN_DECLARED_QUOTIENT');
});
const tailBound=t.pow(6).div(720).div(new F(1).sub(t.pow(2).div(56)));
check('native_weighted_tail_bound_and_bilateral_parity', () => {
  assert(weighted.mass(U).eq(1) && weighted.mass(hp).eq(1));
  assert(tailBound.eq(new F(7,321120)));
  assert(weighted.mass(sumRange(5,12)).le(tailBound));
  for(let j=0;j<series.length;j++) {
    if(j%2) assert.equal(series[j].terms.size,0);
    equal(h,hk.times(p.formalExp(U,j)[j]).times(hk),p.formalExp(U,j)[j].scale(j%2?-1:1));
  }
});
const dependencies={
  CUT:{role:'NATIVE_SCALAR',requires:[]},
  RULES:{role:'NATIVE_PATH',requires:['CUT']},
  MEMORY:{role:'NATIVE_ALGEBRA_DERIVATION',requires:['RULES']},
  FORMAL:{role:'DECLARED_FORMAL_EXTENSION',requires:['MEMORY']}
};
check('no_physical_adapter_in_native_identity_proof', () => {
  assert.deepEqual(gate.nativeDependencies('FORMAL',dependencies),['CUT','RULES','MEMORY','FORMAL']);
  assert.throws(()=>gate.nativeDependencies('FORMAL',{...dependencies,
    CUT:{role:'PHYSICAL_ADAPTER',requires:[]}}),/Non-native premise/);
});
const geometryRun=cp.spawnSync('python',[path.join(here,'geometry_checks.py')],{encoding:'utf8'});
assert.equal(geometryRun.status,0,geometryRun.stderr);
const geometry=JSON.parse(geometryRun.stdout);
const sourceFiles=fs.readdirSync(here).filter(f=>/\.(md|cjs|py|json)$/.test(f)&&f!=='CERTIFICATE.json').sort();
const report=p.stable({protocol:'NATIVE_UGD_PROPAGATION_CORRECTION_R1',
  status:'PASS_SOURCE_BOUND_NATIVE_CONTROLS',checks,
  source_sha256:Object.fromEntries(sourceFiles.map(f=>[f,sha(fs.readFileSync(path.join(here,f)))])),
  consumed_rkf_commit:pins.rkf_commit,consumed_sources_checked:pins.rkf_files.length+pins.publications_files.length,
  native:{EMK_basis:model.basis.words,derived_left_ideal_actions:{R:actionR,K:actionK},
          minimum_added_state_channels:future.extra,coefficient_field:'CUT_COMPLEX; RADIAL FIXTURES ARE REAL',
          full_sheet_normal_language:wb.finiteBasis(h).status,weighted_contract:weighted.contract(),
          recognition_recurrence:'x_(n+2)=-Z^2 x_n; U=R Z is the declared example',
          Morphic_clock_equation:'D_tau^2 x+(tau_double_prime/tau_prime^2) D_tau x+(Z^2/tau_prime^2) x=0',
          recognized_seed_returns:recurrence,proofs,
          proper_seam_boundary:proper.reduce(boundary).normal.toJSON(),
          path_backend_four_steps:pathFour.toJSON()},
  madhava_smriti:{formal_check_modulus:'t^13',correction_at_half:'1/384',
                  infinite_tail_mass_bound_at_half_after_degree_four:tailBound.toString(),
                  tail_bound_source:'WC1/WC2 completion plus the even native return coefficients',
                  physical_clock_interpretation:false},
  geometry,
  scope:{native_arithmetic_reimplemented:false,primitive_Hilbert_space:false,
         chosen_finite_matrices_as_premise:false,canonical_engine_full_master_recertified:false,
         arbitrary_UGD_numeral_multiplication_equals_arrow_composition:false,
         complete_proper_seam_identified_with_sheet_translation:false,
         unique_physical_propagator_selected:false,SI_c_predicted:false,
         quantum_gravity_proved:false,independent_formal_proof:false}});
const content=JSON.stringify(report,null,2)+'\n',digest=sha(content);
if(argv.includes('--write')) {
  fs.writeFileSync(path.join(here,'CERTIFICATE.json'),content);
  fs.writeFileSync(path.join(here,'EXPECTED.sha256'),digest+'\n');
} else {
  assert.equal(fs.readFileSync(path.join(here,'CERTIFICATE.json'),'utf8'),content,'Native certificate mismatch; nothing overwritten');
  assert.equal(fs.readFileSync(path.join(here,'EXPECTED.sha256'),'utf8').trim(),digest);
}
console.log(JSON.stringify({status:report.status,mode:argv.includes('--write')?'write':'check',
                           checks:checks.length,sha256:digest,runtime:process.version,
                           native_engine_imported_from:engine,tail_bound:tailBound.toString()}));
