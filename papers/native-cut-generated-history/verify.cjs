'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const m=N.matrix,mm=(...xs)=>xs.reduce(N.mul),dg=N.dagger,eq=N.equal,sc=N.scale;
const checks=[],negative_controls=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
const I=N.identity(2),A=m([[1,0],[0,0]]),u=m([['3/5'],['4/5']]),B=mm(u,dg(u));
const vi=m([['3/5'],[new N.Cut(0,'4/5')]]),C=mm(vi,dg(vi));
const cut=(P,b)=>b?P:N.sub(I,P),trace=X=>X.reduce((s,r,j)=>s.add(r[j]),N.ZERO);
const outer=x=>mm(x,dg(x));
const normalized=x=>mm(dg(x),x)[0][0].eq(1);
for(const [name,P]of [['A',A],['B',B],['C',C]]){
 ck('cut:'+name+':self_adjoint',eq(P,dg(P)));
 ck('cut:'+name+':idempotent',eq(mm(P,P),P));
 ck('cut:'+name+':binary_partition',eq(N.add(P,cut(P,0)),I)&&N.isZero(mm(P,cut(P,0))));
}
const dimensions=[N.operatorAlgebra([I]).dimension,N.operatorAlgebra([A]).dimension,N.operatorAlgebra([A,B]).dimension];
ck('algebra:scalar_single_and_two_cut_dimensions',JSON.stringify(dimensions)==='[1,2,4]');
ck('algebra:repeated_cut_does_not_enlarge',N.operatorAlgebra([A,A]).dimension===2);
no('algebra:two_outcomes_force_a_change',eq(mm(A,A),A)&&N.isZero(mm(cut(A,0),A)));
no('algebra:cut_binary_algebra_contains_every_operator',!N.isZero(N.commutator(A,B)));

// Prefix-dependent projection oracle; mathematical operations use exact native matrices.
function history(depth,oracle){
 let branches=[{word:'',K:I}];
 for(let j=0;j<depth;j++)branches=branches.flatMap(({word,K})=>[0,1].map(b=>({word:word+b,K:mm(cut(oracle(j,word),b),K)})));
 return {branches,W:branches.flatMap(b=>b.K)};
}
function extension(branches,j,oracle){
 const V=N.zeros(branches.length*4,branches.length*2);
 branches.forEach(({word},k)=>[0,1].forEach(b=>{
  const P=cut(oracle(j,word),b);
  for(let r=0;r<2;r++)for(let c=0;c<2;c++)V[(2*k+b)*2+r][2*k+c]=P[r][c];
 }));
 return V;
}
const fixed=(j)=>[A,B,A,C][j%4];
const adaptive=(j,word)=>j===0?A:(word.endsWith('1')?C:B);
const sources=[m([['3/5'],['4/5']]),m([['3/5'],['-4/5']]),vi,m([['3/5'],[new N.Cut(0,'-4/5')]])];
sources.forEach((x,k)=>ck('source:unit:'+k,normalized(x)));
const records=[];
for(const [protocol,oracle]of [['fixed',fixed],['adaptive',adaptive]]){
 for(let depth=0;depth<=4;depth++){
  const {branches,W}=history(depth,oracle),tag=protocol+':'+depth;
  ck('history:'+tag+':isometry',eq(mm(dg(W),W),I));
  const R=mm(W,dg(W));
  ck('history:'+tag+':image_projection',eq(mm(R,R),R)&&eq(dg(R),R));
  ck('history:'+tag+':operator_roundtrip',eq(mm(dg(W),W,B,dg(W),W),B));
  const JA=mm(W,A,dg(W)),JB=mm(W,B,dg(W));
  ck('history:'+tag+':product_preserved',eq(mm(JA,JB),mm(W,A,B,dg(W))));
  ck('history:'+tag+':adjoint_preserved',eq(dg(mm(W,A,B,dg(W))),mm(W,B,A,dg(W))));
  for(const [k,x]of sources.entries()){
   ck('history:'+tag+':state_recovery:'+k,eq(mm(dg(W),W,x),x));
   const rho=outer(x),encoded=mm(W,rho,dg(W));
   ck('history:'+tag+':density_recovery:'+k,eq(mm(dg(W),encoded,W),rho));
   ck('history:'+tag+':norm_preserved:'+k,mm(dg(mm(W,x)),W,x)[0][0].eq(1));
  }
  if(depth<4){
   const V=extension(branches,depth,oracle),next=history(depth+1,oracle).W;
   ck('history:'+tag+':extension_isometry',eq(mm(dg(V),V),N.identity(W.length)));
   ck('history:'+tag+':cut_only_intertwiner',eq(mm(V,W),next));
   ck('history:'+tag+':algebra_intertwiner',eq(mm(V,JB,dg(V)),mm(next,B,dg(next))));
  }
  if(depth>0){
   no('history:'+tag+':image_unit_equals_full_unit',!eq(R,N.identity(W.length)));
   const rhoPlus=outer(sources[0]),rhoMinus=outer(sources[1]),rhoImag=outer(sources[2]);
   ck('history:'+tag+':all_classical_branch_states_equal',branches.every(({K})=>eq(mm(K,rhoPlus,dg(K)),mm(K,rhoMinus,dg(K)))&&eq(mm(K,rhoPlus,dg(K)),mm(K,rhoImag,dg(K)))));
   no('history:'+tag+':classical_ledger_reconstructs_phase',!eq(rhoPlus,rhoMinus)&&!eq(mm(W,rhoPlus,dg(W)),mm(W,rhoMinus,dg(W))));
  }
  records.push({protocol,depth,branches:branches.length,carrier_rows:W.length,rank:N.rank(W)});
 }
}
// Source fixture: no interposed Hamiltonian, only actual cut interactions.
ck('fixture:CLD_cut',eq(B,m([['9/25','12/25'],['12/25','16/25']])));
const triple=history(3,fixed),e1=m([[1],[0]]),probabilities={};
for(const {word,K}of triple.branches)probabilities[word]=mm(dg(e1),dg(K),K,e1)[0][0];
for(const [word,p]of Object.entries({'111':'81/625','110':'144/625','101':'256/625','100':'144/625'}))ck('fixture:probability:'+word,probabilities[word].eq(p));
ck('fixture:initial_cut_deterministic',Object.entries(probabilities).filter(([w])=>w[0]==='0').every(([,p])=>p.zero()));
const total=Object.values(probabilities).reduce((s,p)=>s.add(p),N.ZERO);
const survival=Object.entries(probabilities).filter(([w])=>w.endsWith('1')).reduce((s,[,p])=>s.add(p),N.ZERO);
ck('fixture:probability_sum_one',total.eq(1));
ck('fixture:final_return_probability',survival.eq('337/625'));
ck('fixture:changed_reading_probability',N.ONE.sub(survival).eq('288/625'));
no('fixture:all_outcomes_return_to_start',!survival.eq(1));
no('fixture:success_normalization_preserves_branch_weight',!probabilities['111'].eq(1));
const rhoPlus=outer(sources[0]),dephased=N.add(mm(A,rhoPlus,A),mm(cut(A,0),rhoPlus,cut(A,0)));
ck('dephasing:lost_part',eq(N.sub(rhoPlus,dephased),N.add(mm(A,rhoPlus,cut(A,0)),mm(cut(A,0),rhoPlus,A))));
ck('dephasing:commuting_source_fixed',eq(N.add(mm(A,dephased,A),mm(cut(A,0),dephased,cut(A,0))),dephased));
no('dephasing:full_source_preserved_by_classicalization',!eq(rhoPlus,dephased));
// Matched labels in the two orders: (a,b), irrespective of chronological storage.
const orderRecords=[];
for(const [name,D]of [['real_B',B],['complex_C',C],['same_A',A],['complement_A',cut(A,0)]]){
 const forward=[],reverse=[];
 for(const a of [0,1])for(const b of [0,1]){forward.push(...mm(cut(D,b),cut(A,a)));reverse.push(...mm(cut(A,a),cut(D,b)));}
 const delta=N.sub(forward,reverse),comm=N.commutator(D,A),gram=mm(dg(delta),delta);
 ck('order:'+name+':exact_four_commutator_identity',eq(gram,sc(mm(dg(comm),comm),4)));
 ck('order:'+name+':both_histories_isometric',eq(mm(dg(forward),forward),I)&&eq(mm(dg(reverse),reverse),I));
 if(name==='real_B'||name==='complex_C'){
  ck('order:'+name+':actual_budget',eq(gram,sc(I,'576/625')));
  no('order:'+name+':arbitrary_cut_order_invariant',!N.isZero(delta));
 }else ck('order:'+name+':commuting_orders_agree',N.isZero(delta));
 orderRecords.push({name,commutator:comm,discrepancy_gram:gram});
}
// Noncommutative projector proofs: no AB=BA relation is assumed.
const system=new M.Presentation({tokens:['p','b','c'],rules:['p','b','c'].map(t=>({id:t+'_idempotent',lhs:[t,t],rhs:[[[t],1]],source:'Declared orthogonal cut idempotence'}))});
const audit=system.audit();
ck('words:confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
presentations.push({name:'noncommuting_cuts',spec:system.spec(),audit});
const one=system.one(),p=system.word(['p']),b=system.word(['b']),c=system.word(['c']);
function prove(name,lhs,rhs){
 const r=system.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);
 const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',system.replay(certificate));
 proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});
}
const q=one.minus(p);
prove('words:binary_complement',p.times(p).plus(q.times(q)),one);
prove('words:orthogonal_sheets',p.times(q).plus(q.times(p)),system.zero());
function normHistory(depth,oracle){
 let branches=[{word:'',K:one,Kd:one}];
 for(let j=0;j<depth;j++)branches=branches.flatMap(({word,K,Kd})=>[0,1].map(bit=>{const P=oracle(j,word),E=bit?P:one.minus(P);return {word:word+bit,K:E.times(K),Kd:Kd.times(E)};}));
 return branches.reduce((sum,h)=>sum.plus(h.Kd.times(h.K)),system.zero());
}
prove('words:fixed_three_cut_history_isometry',normHistory(3,j=>[p,b,p][j]),one);
prove('words:adaptive_three_cut_history_isometry',normHistory(3,(j,w)=>j===0?p:(w.endsWith('1')?c:b)),one);
let deltaGram=system.zero();
for(const E of [one.minus(p),p])for(const F of [one.minus(b),b]){const d=F.times(E).minus(E.times(F));deltaGram=deltaGram.minus(d.times(d));}
const comm=b.times(p).minus(p.times(b));
prove('words:cut_order_discrepancy',deltaGram,comm.times(comm).scale(-4));
no('words:commutation_follows_from_idempotence',system.reduce(comm).normal.terms.size!==0);
const bad=JSON.parse(JSON.stringify(proofs[0].certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_rewrite_witness',rejects(()=>system.replay(bad)));
const binary=require('./binary_synthesis.cjs');
checks.push(...binary.checks);negative_controls.push(...binary.negative_controls);
proofs.push(...binary.proofs);presentations.push(...binary.presentations);
const circular=require('./circular_observer.cjs');
checks.push(...circular.checks);negative_controls.push(...circular.negative_controls);
proofs.push(...circular.proofs);presentations.push(...circular.presentations);
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'CHD1_CUT_GENERATED_HISTORY_V1',status:'PASS_CUT_GENERATED_HISTORY',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results:{observable_algebra_dimensions:dimensions,history_records:records,source_cut:B,probabilities,final_return_probability:survival,order_records:orderRecords,scope:{uncut_scalar_algebra:'declared observable contract',cut_schedule:'declared',projective_implementation:'selected standard instrument',precut_source_evolution_used:false,physical_time_derived:false,arbitrary_depth:'written induction; finite certificates separately scoped',full_UGD_homomorphism:false,spacetime_dimension_derived:false,quantum_gravity_complete:false}},summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
payload.protocol='CHD3_COMPLETE_CIRCULAR_OBSERVER_VERIFICATION_V1';
payload.results.binary_synthesis=binary.results;
payload.results.circular_observer=circular.results;
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
