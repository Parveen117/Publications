#!/usr/bin/env node
'use strict';
// Consumer of the canonical RKF engine. No native arithmetic/rewrite fork.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const home=__dirname,root=path.resolve(home,'../..');
const arg=n=>{const i=process.argv.indexOf(n);if(i<0||!process.argv[i+1])throw Error('Missing '+n);return path.resolve(process.argv[i+1]);};
const rkf=arg('--rkf-root'),thermo=arg('--thermo-root');
const ensure=(v,m)=>{if(!v)throw Error(m);};
ensure(process.argv.includes('--check')!==process.argv.includes('--write'),'Choose exactly one of --check or --write');
const sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const blob=b=>crypto.createHash('sha1').update(Buffer.concat([Buffer.from('blob '+b.length+'\0'),b])).digest('hex');
const pinsBytes=fs.readFileSync(path.join(home,'SOURCE_PINS.json')),pins=JSON.parse(pinsBytes);
function checkPins(base,entries){for(const e of entries){const f=path.resolve(base,e.path);ensure(f.startsWith(base+path.sep),'Invalid source path');const b=fs.readFileSync(f);ensure(sha(b)===e.sha256&&blob(b)===e.git_blob_sha,'Source pin mismatch: '+e.path);}}
checkPins(root,pins.local_inputs);checkPins(root,pins.reused_publication_inputs);
checkPins(rkf,pins.rkf.files);checkPins(thermo,pins.thermo.files);
const p=require(path.join(rkf,'operator_foundation/core/paninian_operator.cjs'));
const o=require(path.join(rkf,'operator_foundation/core/native_operator.cjs'));
const system=new p.Presentation(JSON.parse(fs.readFileSync(path.join(rkf,'operator_foundation/examples/emk_job.json'))).presentation);
ensure(system.audit().status==='CONFLUENT_BY_CHECKED_DIAMONDS','Presentation audit failed');
const N=x=>system.reduce(x).normal,one=system.one(),zero=system.zero();
const add=(a,b)=>N(a.plus(b)),sub=(a,b)=>N(a.minus(b));
const mul=(a,b)=>N(a.times(b)),scale=(a,c)=>N(a.scale(c));
const sum=xs=>xs.reduce(add,zero),comm=(a,b)=>N(p.bracket(a,b));
const R=system.word(['R']),K=system.word(['K']),L=mul(R,K);
const basis=[one,K,R,L],star=x=>N(p.star(x,{R:scale(R,-1),K}));
const sc=x=>N(x).terms.get('[]')||o.ZERO;
const f=x=>o.F.of(x),half=f('1/2');
const e=scale(add(one,K),half),q=scale(sub(one,K),half),e12=scale(sub(L,R),half);
const sym=h=>sum([scale(one,f(h[0][0]).add(h[1][1]).div(2)),scale(K,f(h[0][0]).sub(h[1][1]).div(2)),scale(L,h[0][1])]);
const inverseH=h=>{
  const a=f(h[0][0]),b=f(h[0][1]),c=f(h[1][1]),d=a.mul(c).sub(b.mul(b));
  ensure(!d.zero()&&!a.le(0)&&!d.le(0),'Not a stable response element');
  return scale(sum([scale(one,a.add(c).div(2)),scale(K,a.sub(c).div(-2)),scale(L,b.neg())]),f(1).div(d));
};
const repr=x=>{const n=N(x),co=w=>n.terms.get(JSON.stringify(w))||o.ZERO;
 const a=co([]),k=co(['K']),r=co(['R']),l=co(['R','K']);
 return o.matrix([[a.add(k),l.sub(r)],[l.add(r),a.sub(k)]]);};
let identities=0;
const groups=[],replays=[],negative={};
function eq(a,b,label,keep=false){ensure(N(a.minus(b)).terms.size===0,'Identity: '+label);identities++;
 if(keep){const v=p.proveEquality(a,b);ensure(v.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(v.certificate),'Replay: '+label);replays.push({name:label,certificate:v.certificate});}}
function neq(a,b,label){ensure(N(a.minus(b)).terms.size!==0,'Missing counterexample: '+label);negative[label]=true;}
function meq(a,b,label){ensure(o.equal(a,b),'Matrix identity: '+label);identities++;}
function group(name,fn){const before=identities;groups.push({name,...(fn()||{}),exact_identities:identities-before});}
// verify.sh runs the independent Python route afresh and pipes its result here.
const controls=JSON.parse(fs.readFileSync(0,'utf8'));
ensure(controls.status==='PASS_INDEPENDENT_RATIONAL_CONTROLS','Independent controls failed');
ensure(controls.fixtures.length===48,'Independent fixture coverage changed');

group('native_carrier_inverse_and_positive_factor',()=>{
 ensure(p.daggerGate(system,{R:scale(R,-1),K}).status==='DAGGER_DESCENDS','Dagger gate');
 eq(R.times(R),scale(one,-1),'RR',true);eq(K.times(K),one,'KK',true);
 eq(p.bracket(K,L),scale(R,-2),'KL orientation',true);
 eq(p.bracket(R,L),scale(K,-2),'diagonal commutator curvature',true);
 for(const a of basis)for(const b of basis){
  ensure(sc(mul(a,b)).eq(sc(mul(b,a))),'Cyclic scalar part');identities++;
  meq(repr(mul(a,b)),o.mul(repr(a),repr(b)),'faithful multiplication');
 }
 for(let n=1;n<=12;n++){
  const a=f(n+1),b=f(n%3-1),c=f(n%4+1);
  const B=sum([scale(e,a),scale(e12,b),scale(q,c)]);
  const H=sym([[a.mul(a),a.mul(b)],[a.mul(b),b.mul(b).add(c.mul(c))]]);
  eq(mul(star(B),B),H,'native Cholesky factor');
  const x=sum([scale(one,n),K,scale(R,-2),L]);
  ensure(sc(mul(star(x),x)).eq(n*n+6),'Native positive scalar part');identities++;
  const hi=inverseH([[a.mul(a),a.mul(b)],[a.mul(b),b.mul(b).add(c.mul(c))]]);
  eq(H.times(hi),one,'response inverse',n===1);eq(hi.times(H),one,'response inverse reverse');
 }
 return {factor_cases:12,representation_basis_products:16};
});

const nativeFixtures=[];
group('native_response_connection_and_independent_Riemann_adapter',()=>{
 for(const z of controls.fixtures){
  const H=sym(z.H),hi=inverseH(z.H),hs=z.dH.map(sym),x=hs.map(h=>mul(hi,h));
  const A=x.map(v=>scale(v,half));
  const dx=(i,j)=>add(scale(mul(x[i],x[j]),-1),mul(hi,sym(z.ddH[i][j])));
  const direct=add(scale(sub(dx(0,1),dx(1,0)),half),comm(A[0],A[1]));
  const F=scale(comm(x[0],x[1]),'-1/4');
  eq(direct,F,'native differentiated connection');
  meq(repr(F),o.matrix(z.curvature),'independent full Christoffel route');
  eq(F,scale(mul(R,H),f(z.kappa).neg()),'native curvature line');
  ensure(sc(mul(mul(R,H),F)).div(z.delta).eq(z.kappa),'signed native marker');identities++;
  ensure(sc(F).zero(),'Curvature scalar trace');identities++;
  for(let i=0;i<2;i++){
   eq(add(mul(star(A[i]),H),mul(H,A[i])),hs[i],'metric compatibility');
   eq(mul(star(A[i]),H),mul(H,A[i]),'strain selection');
   meq(repr(A[i]),o.matrix(z.Gamma[i]),'connection intertwiner');
  }
  for(const t of ['0','1/3','1/2','1','2']){
   const n=f(t),tangent=add(scale(sub(dx(0,1),dx(1,0)),n),scale(comm(x[0],x[1]),n.mul(n)));
   eq(tangent,scale(comm(x[0],x[1]),n.mul(n.sub(1))),'q(q-1) family');
  }
  nativeFixtures.push({H,x,A,F});
 }
 const t=nativeFixtures[0];
 eq(t.F,sum([scale(K,'5/324'),scale(R,'-35/648'),scale(L,'5/216')]),'curved fixture coefficients',true);
 eq(nativeFixtures[1].F,zero,'quadratic control');
 neq(t.F,nativeFixtures[1].F,'same_point_six_responses_do_not_determine_curvature');
 neq(t.F,scale(comm(t.x[0],t.x[1]),'1/4'),'wrong_curvature_sign_rejected');
 neq(t.F,comm(t.A[0],t.A[1]),'connection_commutator_without_derivatives_rejected');
 negative.scalar_ratio_and_trace_flatness_inference_rejected=!sc(t.F).zero()?false:t.F.terms.size>0;
 return {independent_Christoffel_cases:48,q_family_cases:240,fourth_order_jets_retained:true};
});

group('native_phase_coframe_and_parallel_quarter_turn',()=>{
 for(let n=1;n<=12;n++){
  const m=f(n+1),h=f(n%3-1),t=f(n%4+1),a=m.mul(m),b=m.mul(h),c=h.mul(h).add(t.mul(t));
  const H=sym([[a,b],[b,c]]),hi=inverseH([[a,b],[b,c]]),sqrtD=m.mul(t),d=sqrtD.mul(sqrtD);
  const B=sum([scale(e,m),scale(e12,h),scale(q,t)]);
  const bi=sum([scale(e,f(1).div(m)),scale(e12,h.neg().div(m.mul(t))),scale(q,f(1).div(t))]);
  const J=scale(mul(R,H),f(1).div(sqrtD));
  eq(J.times(J),scale(one,-1),'native represented quarter turn',n===1);
  eq(mul(B,bi),one,'coframe inverse');
  eq(mul(mul(B,J),bi),R,'phase intertwiner');
  for(const jets of [[1,2,-1],[-2,1,3]]){
   const [da,db,dc]=jets.map(f),dh=sym([[da,db],[db,dc]]),dd=da.mul(c).add(a.mul(dc)).sub(b.mul(db).mul(2));
   const A=scale(mul(hi,dh),half);
   const dB=sum([scale(e,da.div(m.mul(2))),
     scale(e12,db.div(m).sub(b.mul(da).div(m.pow(3).mul(2)))),
     scale(q,dd.div(m.mul(m).mul(t).mul(2)).sub(t.mul(da).div(m.mul(m).mul(2))))]);
   const omega=sub(mul(mul(B,A),bi),mul(dB,bi));
   const alpha=db.sub(b.mul(da).div(a)).div(sqrtD.mul(2));
   eq(omega,scale(R,alpha),'native phase one-form');
   const dj=sub(scale(mul(R,dh),f(1).div(sqrtD)),scale(J,dd.div(d.mul(2))));
   eq(add(dj,comm(A,J)),zero,'parallel J');
  }
 }
 const H=nativeFixtures[0].H,F=nativeFixtures[0].F,J=scale(mul(R,H),'1/3');
 ensure(sc(mul(J,F)).neg().eq('-5/108'),'oriented phase marker');identities++;
 neq(F,scale(F,-1),'opposite_oriented_curvatures_are_distinct');
 return {coframe_cases:12,derivative_directions:24};
});

// Inner dagger-compatible frame, using the native central cut scalar.
const G=[R,scale(K,o.IOTA),scale(L,o.IOTA)];
const c=Array.from({length:3},()=>Array.from({length:3},()=>[0,0,0]));
for(const [i,j,k]of [[0,1,2],[1,2,0],[2,0,1]]){c[i][j][k]=2;c[j][i][k]=-2;}
const delta=(i,x)=>comm(G[i],x),frame=(i,j,A)=>sum(c[i][j].map((n,k)=>scale(A[k],n)));
const curv=(A,i,j)=>sub(add(sub(delta(i,A[j]),delta(j,A[i])),comm(A[i],A[j])),frame(i,j,A));
group('noncommuting_frame_and_ledger_balance',()=>{
 const H=nativeFixtures[0].H,hi=inverseH(controls.fixtures[0].H);
 const X=G.map((_,i)=>mul(hi,delta(i,H))),A=X.map(x=>scale(x,half));
 for(let i=0;i<3;i++){
  eq(star(G[i]),scale(G[i],-1),'dagger-compatible frame');
  for(let j=0;j<3;j++){
   eq(comm(G[i],G[j]),frame(i,j,G),'native frame constants');
   eq(curv(A,i,j),scale(comm(X[i],X[j]),'-1/4'),'nonholonomic response curvature');
   eq(curv(X,i,j),zero,'native Maurer-Cartan flatness');
  }
 }
 const u=add(one,R),ui=scale(sub(one,R),'1/2');
 for(let n=0;n<12;n++){
  const A=[sum([R,scale(K,n)]),sum([L,scale(R,n-2)]),sum([K,scale(L,2-n)])];
  const M=[sum([scale(L,n+1),K]),sum([R,scale(K,1-n)]),sum([L,scale(R,n)])];
  const C=A.map((a,i)=>sub(a,M[i]));
  const Ap=A.map((a,i)=>add(mul(mul(ui,a),u),mul(ui,delta(i,u))));
  const Mp=M.map(m=>mul(mul(ui,m),u)),Cp=Ap.map((a,i)=>sub(a,Mp[i]));
  for(let i=0;i<3;i++)for(let j=0;j<3;j++){
   const DM=add(sub(sub(delta(i,M[j]),delta(j,M[i])),frame(i,j,M)),sub(comm(A[i],M[j]),comm(A[j],M[i])));
   eq(curv(C,i,j),add(sub(curv(A,i,j),DM),comm(M[i],M[j])),'full Smriti balance');
   eq(curv(Cp,i,j),mul(mul(ui,curv(C,i,j)),u),'ledger gauge covariance');
  }
 }
 const raw=comm(K,L),ms=zero,mv=L;
 const DM=sub(comm(K,mv),comm(L,ms));
 eq(sub(raw,DM),zero,'balanced open connection',true);
 neq(raw,zero,'balanced_flatness_does_not_erase_raw_curvature');
 neq(sub(raw,comm(ms,mv)),zero,'subtracting_two_raw_curvatures_is_not_ledger_balance');
 const missingFrame=add(sub(delta(0,A[1]),delta(1,A[0])),comm(A[0],A[1]));
 neq(missingFrame,curv(A,0,1),'omitted_frame_term_rejected');
 return {nonholonomic_pairs:9,ledger_samples:12,gauge_pairs:108};
});

group('compression_return_and_moving_cut',()=>{
 for(const a of basis)for(const b of basis){
  const raw=comm(a,b),compressed=comm(mul(mul(e,a),e),mul(mul(e,b),e));
  const hidden=sub(mul(mul(mul(mul(e,a),q),b),e),mul(mul(mul(mul(e,b),q),a),e));
  eq(compressed,sub(mul(mul(e,raw),e),hidden),'native Gauss return formula');
 }
 const raw=comm(R,L),hidden=sub(mul(mul(mul(mul(e,R),q),L),e),mul(mul(mul(mul(e,L),q),R),e));
 eq(mul(mul(e,raw),e),scale(e,-2),'retained raw curvature',true);
 eq(hidden,scale(e,-2),'complement return',true);
 eq(sub(mul(mul(e,raw),e),hidden),zero,'compressed flat verdict');
 neq(mul(mul(e,raw),e),zero,'omitting_return_term_changes_observer_curvature');
 neq(mul(mul(e,R),q),zero,'compression_is_not_quotient_descent');
 // P(s,v)=g e g^-1, g=Exp(s R) Exp(v iota L). Ambient A=0.
 const Z=scale(L,o.IOTA),ps=comm(R,e),pv=comm(Z,e);
 const grass=mul(mul(e,comm(ps,pv)),e);
 const fixedFrame=scale(mul(mul(e,comm(R,Z)),e),-1); // -partial_v a_s
 eq(grass,fixedFrame,'moving-cut two-route curvature');
 eq(grass,scale(e,o.IOTA.mul(2)),'moving cut in flat ambient carrier',true);
 neq(grass,zero,'flat_ambient_connection_need_not_make_compression_flat');
 const rhoR=repr(R),readout=o.matrix([[1,0]]);
 ensure(o.descendedAction(rhoR,readout)===null,'Blind channel must fail descent');identities++;
 return {constant_cut_pairs:16,moving_unitary_cut_cases:1,failed_descent_witnesses:1};
});

group('ordered_plaquette_cut_loop_and_tower_action',()=>{
 for(const a of [R,K,L,add(R,K)])for(const b of [R,K,L]){
  const exp=x=>p.formalExp(x,3);
  const product=(a,b)=>p.seriesMultiply(a,b,3);
  const loop=product(product(product(exp(b),exp(a)),exp(scale(b,-1))),exp(scale(a,-1)));
  eq(loop[1],zero,'first order rectangle cancellation');
  eq(loop[2],scale(comm(a,b),-1),'positive ordered rectangle sign');
  const G=add(a,b),even=scale(add(G,mul(mul(K,G),K)),half),odd=sub(G,even);
  const cut=product(exp(sub(even,odd)),exp(add(even,odd))),log=p.formalLog(cut,3);
  eq(log[1],scale(even,2),'cut-loop linear coefficient');
  eq(log[2],comm(even,odd),'cut-loop second coefficient');
 }
 const ae=repr(K),be=repr(R),at=repr(R),bt=repr(L),I=o.identity(2);
 const da=o.sub(o.kron(ae,I),o.kron(I,o.dagger(at)));
 const db=o.sub(o.kron(be,I),o.kron(I,o.dagger(bt)));
 const expected=o.sub(o.kron(o.commutator(ae,be),I),o.kron(I,o.dagger(o.commutator(at,bt))));
 meq(o.commutator(da,db),expected,'response fibre and dual-slot curvature');
 negative.reversed_plaquette_sign_rejected=comm(R,K).terms.size>0;
 return {ordered_word_cases:12,tensor_dual_slot_cases:1,formal_order:3};
});

group('adversarial_proof_and_artifact_controls',()=>{
 const witness=replays.find(x=>x.certificate.steps.length);ensure(witness,'No nontrivial rewrite witness');
 const broken=JSON.parse(JSON.stringify(witness.certificate));broken.steps[0].rule='__corrupted__';
 try{system.replay(broken);negative.corrupted_rewrite_rejected=false;}catch{negative.corrupted_rewrite_rejected=true;}
 try{checkPins(rkf,[{...pins.rkf.files[0],sha256:'0'.repeat(64)}]);negative.corrupted_pin_rejected=false;}catch{negative.corrupted_pin_rejected=true;}
 try{inverseH([[1,1],[1,1]]);negative.degenerate_response_rejected=false;}catch{negative.degenerate_response_rejected=true;}
 ensure(Object.values(negative).every(Boolean),'A false alternative survived');
 return {rejected_alternatives:Object.keys(negative).length};
});
const text=fs.readFileSync(path.join(home,'THEOREM.md'),'utf8');
const sections=[...text.matchAll(/## NT-([1-8])\. ([^\n]+)\n([\s\S]*?)(?=\n## |$)/g)];
ensure(sections.length===8&&sections.every((s,i)=>+s[1]===i+1&&s[3].includes('**Proof.**')),'Written theorem/proof coverage mismatch');
const result={schema:'publications.native-thermodynamic-curvature.v1',status:'PASS_NT1_NT8_FINITE_CONTROLS',
 source_pins_sha256:sha(pinsBytes),canonical_engine_commit:pins.rkf.commit,presentation_sha256:system.hash(),
 exact_native_and_matrix_identities:identities,groups,word_replay_count:replays.length,word_replays:replays,
 independent_controls:{sha256:sha(Buffer.from(JSON.stringify(controls))),checks:controls.checks,
  rejected_false_alternatives:controls.rejected_false_alternatives},
 rejected_false_alternatives:negative,
 written_theorems:sections.map(s=>({id:'NT-'+s[1],title:s[2],sha256:sha(Buffer.from(s[0]))})),
 fixture:{native_coefficients:{K:'5/324',R:'-35/648',RK:'5/216'},Gaussian_marker:'5/324',phase_curvature:'-5/108',scalar_ratio_closure:'1'},
 scope:{written_proofs:true,finite_exact_native_execution:true,independent_full_Christoffel_route:true,
  original_register_wholesale_recertified:false,proof_assistant_formalized:false,external_review:false,
  new_empirical_validation:false,universal_thermodynamic_constitutive_law_derived:false,
  local_artifact_hashes_checked:true}};
const bytes=Buffer.from(JSON.stringify(p.stable(result),null,2)+'\n'),cert=path.join(home,'CERTIFICATE.json'),expected=path.join(home,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cert,bytes);fs.writeFileSync(expected,sha(bytes)+'\n');}
else ensure(fs.readFileSync(cert).equals(bytes)&&fs.readFileSync(expected,'utf8').trim()===sha(bytes),'Frozen certificate mismatch');
console.log(JSON.stringify({status:result.status,groups:groups.length,exact_identities:identities,
 word_replays:replays.length,independent_Christoffel_cases:controls.fixtures.length,
 negative_controls:Object.keys(negative).length+controls.rejected_false_alternatives.length,
 certificate_sha256:sha(bytes)}));
