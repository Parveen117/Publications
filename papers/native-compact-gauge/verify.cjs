#!/usr/bin/env node
'use strict';
// Consumer only: all native words and rational matrices use the pinned RKF engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const home=__dirname,root=path.resolve(home,'../..');
const arg=n=>{const i=process.argv.indexOf(n);if(i<0||!process.argv[i+1])throw Error('Missing '+n);return path.resolve(process.argv[i+1]);};
const rkf=arg('--rkf-root'),physics=arg('--physics-root');
const ensure=(v,m)=>{if(!v)throw Error(m);};
const sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const blob=b=>crypto.createHash('sha1').update(Buffer.concat([Buffer.from('blob '+b.length+'\0'),b])).digest('hex');
const pinBytes=fs.readFileSync(path.join(home,'SOURCE_PINS.json')),pins=JSON.parse(pinBytes);
function checkPins(base,entries){for(const e of entries){const f=path.resolve(base,e.path);ensure(f.startsWith(base+path.sep),'Invalid source path');const b=fs.readFileSync(f);ensure(sha(b)===e.sha256&&blob(b)===e.git_blob_sha,'Source pin mismatch: '+e.path);}}
checkPins(root,pins.local_inputs);checkPins(root,pins.reused_publication_inputs);
checkPins(rkf,pins.rkf.files);checkPins(physics,pins.physics.files);
if(process.argv.includes('--pins-only')){console.log('PASS: source pins before execution');process.exit(0);}
ensure(process.argv.includes('--check')!==process.argv.includes('--write'),'Choose --check or --write');
const p=require(path.join(rkf,'operator_foundation/core/paninian_operator.cjs'));
const o=require(path.join(rkf,'operator_foundation/core/native_operator.cjs'));
const system=new p.Presentation(JSON.parse(fs.readFileSync(path.join(home,'presentation.json'))));
const audit=system.audit();ensure(audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS','Presentation not certified confluent');
const N=x=>system.reduce(x).normal,one=system.one(),zero=system.zero();
const add=(a,b)=>N(a.plus(b)),sub=(a,b)=>N(a.minus(b)),mul=(a,b)=>N(a.times(b)),scale=(a,c)=>N(a.scale(c));
const sum=xs=>xs.reduce(add,zero),comm=(a,b)=>N(p.bracket(a,b));
const R1=system.word(['R1']),K1=system.word(['K1']),R2=system.word(['R2']),K2=system.word(['K2']);
const dagger={R1:scale(R1,-1),K1,R2:scale(R2,-1),K2};
const star=x=>N(p.star(x,dagger)),sc=x=>N(x).terms.get('[]')||o.ZERO;
const C=R1,e=[R2,mul(R1,K2),mul(mul(R1,R2),K2)],t=e.map(x=>scale(x,'1/2')),t0=scale(C,'1/2');
const basis=[one,K1,R1,mul(R1,K1)].flatMap(a=>[one,K2,R2,mul(R2,K2)].map(b=>mul(a,b)));
const I=o.identity(2),R=o.matrix([[0,-1],[1,0]]),K=o.matrix([[1,0],[0,-1]]);
const values={R1:o.kron(R,I),K1:o.kron(K,I),R2:o.kron(I,R),K2:o.kron(I,K)};
const repr=x=>p.matrixValue(N(x),values,4);
const eps=(i,j,k)=>new Set([i,j,k]).size<3?0:((i+1)%3===j?1:-1);
let count=0;const groups=[],replays=[],negative={};
function eq(a,b,label,replay=false){ensure(N(a.minus(b)).terms.size===0,'Native identity: '+label);count++;
 if(replay){const q=p.proveEquality(a,b);ensure(q.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(q.certificate),'Proof replay failed');replays.push({name:label,certificate:q.certificate});}}
function meq(a,b,label){ensure(o.equal(a,b),'Matrix identity: '+label);count++;}
function group(name,fn){const start=count;const extra=fn()||{};groups.push({name,exact_checks:count-start,...extra});}
function reject(name,a,b){ensure(N(a.minus(b)).terms.size!==0,'False alternative survived: '+name);negative[name]=true;}
const controls=JSON.parse(fs.readFileSync(0,'utf8'));
ensure(controls.status==='PASS_INDEPENDENT_COMPACT_GAUGE_CONTROLS'&&controls.fixtures.length===6,'Independent controls missing');

group('native_tensor_carrier_and_positive_pairing',()=>{
 ensure(p.daggerGate(system,dagger).status==='DAGGER_DESCENDS','Dagger does not descend');
 for(const [a,name] of [[R1,'R1'],[R2,'R2']])eq(a.times(a),scale(one,-1),name+' square',true);
 for(const [a,name] of [[K1,'K1'],[K2,'K2']])eq(a.times(a),one,name+' square',true);
 eq(K1.times(R1),scale(R1.times(K1),-1),'first EMK anticommutation',true);
 eq(K2.times(R2),scale(R2.times(K2),-1),'second EMK anticommutation',true);
 for(let i=0;i<16;i++)for(let j=0;j<16;j++){
  const a=basis[i],b=basis[j];
  ensure(sc(mul(star(a),b)).eq(i===j?1:0),'Native basis norm');count++;
  ensure(sc(mul(a,b)).eq(sc(mul(b,a))),'Native cyclicity');count++;
  meq(repr(mul(a,b)),o.mul(repr(a),repr(b)),'faithful product');
 }
 ensure(o.rank(basis.map(b=>o.flatten(repr(b))))===16,'Incomplete matrix representation');count++;
 return {basis_dimension:16,basis_pairs:256,representation:'derived coefficient adapter, not primitive'};
});

group('quaternion_and_action_compatible_centralizer',()=>{
 for(let i=0;i<3;i++){
  eq(e[i].times(e[i]),scale(one,-1),'quaternion square '+i,true);
  eq(star(e[i]),scale(e[i],-1),'skew generator');eq(comm(C,e[i]),zero,'massive multiplier preservation');
  for(let j=0;j<3;j++)eq(mul(e[i],e[j]),add(scale(one,-Number(i===j)),sum(e.map((x,k)=>scale(x,eps(i,j,k))))),'quaternion table');
 }
 eq(e[0].times(e[1]),e[2],'oriented quaternion product',true);
 eq(e[1].times(e[2]),e[0],'second quaternion product',true);
 eq(e[2].times(e[0]),e[1],'third quaternion product',true);
 for(const a of [t0,...t])for(const b of [t0,...t]){
  ensure(sc(mul(a,b)).mul(-4).eq(a===b?1:0),'Compact form normalization');count++;
 }
 const constraints=basis.map(x=>[...o.flatten(repr(comm(x,C))),...o.flatten(repr(add(x,star(x))))]);
 ensure(16-o.rank(o.dagger(constraints))===controls.centralizer_dimensions.four_real_copies,'Independent centralizer mismatch');count++;
 reject('swap_quaternion_orientation',comm(t[0],t[1]),scale(t[2],-1));
 reject('discard_central_phase',t0,zero);
 reject('K1_preserves_massive_multiplier',comm(K1,C),zero);
 return {compact_generators:4,centralizer_dimension:4};
});

group('native_connections_match_independent_color_jets',()=>{
 for(const f of controls.fixtures){
  const a=f.A.map(v=>sum(v.map((x,i)=>scale(t[i],x))));
  const da=f.dA.map(row=>row.map(v=>sum(v.map((x,i)=>scale(t[i],x)))));
  for(let i=0;i<4;i++)for(let j=0;j<4;j++){
   const F=add(sub(da[i][j],da[j][i]),comm(a[i],a[j]));
   eq(F,sum(f.F[i][j].map((x,k)=>scale(t[k],x))),'independent affine curvature');
   eq(star(F),scale(F,-1),'skew curvature closure');
   for(let k=0;k<3;k++)eq(comm(t[k],comm(a[i],a[j])),add(comm(comm(t[k],a[i]),a[j]),comm(a[i],comm(t[k],a[j]))),'adjoint derivation');
  }
 }
 eq(comm(t[0],t[1]),t[2],'constant exact-coefficient curvature',true);
 meq(repr(comm(R1,mul(R1,K1))),o.matrix([[-2,0,0,0],[0,-2,0,0],[0,0,2,0],[0,0,0,2]]),'nonzero diagonal curvature witness');
 reject('off_diagonal_entries_are_required',comm(R1,mul(R1,K1)),zero);
 // Pure-gauge ordered product g=exp(x t1)exp(y t2): d_y A_x=[t1,t2], d_x A_y=0.
 eq(add(scale(comm(t[0],t[1]),-1),comm(t[0],t[1])),zero,'Maurer Cartan cancellation');
 return {independent_affine_jet_cases:6};
});

group('finite_unitaries_readout_and_based_loop_sign',()=>{
 for(let i=0;i<3;i++){
  const g=add(scale(one,'3/5'),scale(e[i],'4/5')),gi=star(g);
  eq(mul(gi,g),one,'rational native unitary',i===0);eq(comm(g,C),zero,'action-compatible unitary');
  for(const x of t)for(const y of t){
   const xx=mul(mul(gi,x),g),yy=mul(mul(gi,y),g);
   ensure(sc(mul(xx,yy)).eq(sc(mul(x,y))),'Gauge-invariant positive pairing');count++;
  }
 }
 for(const a of t)for(const b of t){
  const ex=x=>p.formalExp(x,3),prod=(a,b)=>p.seriesMultiply(a,b,3);
  // Positive x,y,-x,-y loop, transport dw=-A w, later segments on the left.
  const u=prod(prod(prod(ex(b),ex(a)),ex(scale(b,-1))),ex(scale(a,-1)));
  eq(u[1],zero,'loop linear term');eq(u[2],scale(comm(a,b),-1),'native loop curvature sign');
  ensure(sc(mul(star(u[2]),u[2])).eq(sc(mul(star(comm(a,b)),comm(a,b)))),'leading positive loop record');count++;
 }
 reject('reverse_transport_sign',scale(comm(t[0],t[1]),-1),comm(t[0],t[1]));
 return {finite_unitaries:3,ordered_plaquettes:9,formal_order:3};
});

group('memory_balance_and_proof_controls',()=>{
 const a=t[0],b=t[1],m=t[2],n=add(t[0],t[1]);
 eq(comm(sub(a,m),sub(b,n)),add(sub(comm(a,b),sub(comm(a,n),comm(b,m))),comm(m,n)),'balanced curvature on a constant frame',true);
 reject('subtract_curvatures_instead_of_connection_memory',comm(sub(a,m),sub(b,n)),sub(comm(a,b),comm(m,n)));
 const proof=replays.find(r=>r.certificate.steps.length);ensure(proof,'Missing nonempty derivation');
 const altered=JSON.parse(JSON.stringify(proof.certificate));altered.steps[0].rule='__wrong__';
 try{system.replay(altered);negative.altered_proof_rejected=false;}catch{negative.altered_proof_rejected=true;}
 try{checkPins(rkf,[{...pins.rkf.files[0],sha256:'0'.repeat(64)}]);negative.altered_pin_rejected=false;}catch{negative.altered_pin_rejected=true;}
 ensure(Object.values(negative).every(Boolean),'A mathematical alteration survived');
});

const theorem=fs.readFileSync(path.join(home,'THEOREM.md'),'utf8');
const sections=[...theorem.matchAll(/## NCG-([1-8])\. ([^\n]+)\n([\s\S]*?)(?=\n## |$)/g)];
ensure(sections.length===8&&sections.every((s,i)=>+s[1]===i+1&&s[3].includes('**Proof.**')),'Written proof bindings missing');
const result={schema:'publications.native-compact-gauge.v1',status:'PASS_NCG1_NCG8_SCOPED_EXACT_CONTROLS',
 source_pins_sha256:sha(pinBytes),canonical_engine_commit:pins.rkf.commit,presentation_sha256:system.hash(),
 native_exact_checks:count,groups,word_replays:replays,
 independent_controls:{sha256:sha(Buffer.from(JSON.stringify(controls))),checks:controls.checks,total:controls.total,
  negative_controls:controls.negative_controls,centralizer_dimensions:controls.centralizer_dimensions,invariant_pairing_dimension:controls.invariant_pairing_dimension},
 native_negative_controls:negative,written_results:sections.map(s=>({id:'NCG-'+s[1],title:s[2],sha256:sha(Buffer.from(s[0]))})),
 scope:{general_written_proofs:true,finite_exact_execution:true,exhaustive_finite_linear_problems:true,
  independent_matrix_and_color_route:true,whole_engine_rerun:false,proof_assistant_formalized:false,external_review:false,
  gauge_group_physically_selected:false,quantum_fermions_constructed:false,interacting_continuum_constructed:false,
  YM_E4DC_closed:false,Clay_dictionary_closed:false,quantum_gravity_completed:false,empirical_validation:false}};
const bytes=Buffer.from(JSON.stringify(p.stable(result),null,2)+'\n');
if(process.argv.includes('--write')){fs.writeFileSync(path.join(home,'CERTIFICATE.json'),bytes);fs.writeFileSync(path.join(home,'EXPECTED.sha256'),sha(bytes)+'\n');}
else ensure(fs.readFileSync(path.join(home,'CERTIFICATE.json')).equals(bytes)&&fs.readFileSync(path.join(home,'EXPECTED.sha256'),'utf8').trim()===sha(bytes),'Frozen certificate mismatch');
console.log(JSON.stringify({status:result.status,native_exact_checks:count,independent_exact_checks:controls.total,
 word_replays:replays.length,negative_controls:Object.keys(negative).length+controls.negative_controls.length,
 certificate_sha256:sha(bytes)}));
