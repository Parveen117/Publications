#!/usr/bin/env node
'use strict';
// Publication adapter only. All arithmetic, normal forms and proof replays
// come from the separately pinned, unchanged canonical RKF engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const arg=name=>{const i=process.argv.indexOf(name);if(i<0||!process.argv[i+1])throw Error('Missing '+name);return process.argv[i+1];};
const root=path.resolve(__dirname,'../../..'),rkf=path.resolve(arg('--rkf-root'));
const hash=b=>crypto.createHash('sha256').update(b).digest('hex');
const blob=b=>crypto.createHash('sha1').update(Buffer.concat([Buffer.from('blob '+b.length+'\0'),b])).digest('hex');
const requireTrue=(v,m)=>{if(!v)throw Error(m);};
const pinsPath=path.join(__dirname,'EMKC1_SOURCE_PINS.json'),pinsBytes=fs.readFileSync(pinsPath),pins=JSON.parse(pinsBytes);
function checkFiles(home,rows){
  for(const row of rows){
    const file=path.resolve(home,row.path);
    requireTrue(file.startsWith(home+path.sep),'Source outside declared root');
    const bytes=fs.readFileSync(file);
    requireTrue(hash(bytes)===row.sha256&&blob(bytes)===row.git_blob_sha,'Source pin mismatch: '+row.path);
  }
}
checkFiles(rkf,pins.canonical_engine.files);checkFiles(root,pins.local_inputs);
const p=require(path.join(rkf,'operator_foundation/core/paninian_operator.cjs'));
const o=require(path.join(rkf,'operator_foundation/core/native_operator.cjs'));
const spec=JSON.parse(fs.readFileSync(path.join(rkf,'operator_foundation/examples/emk_job.json'))).presentation;
const system=new p.Presentation(spec);
requireTrue(system.audit().status==='CONFLUENT_BY_CHECKED_DIAMONDS','Uncertified presentation');
const N=x=>system.reduce(x).normal;
const add=(a,b)=>N(a.plus(b)),sub=(a,b)=>N(a.minus(b)),mul=(a,b)=>N(a.times(b)),scale=(a,c)=>N(a.scale(c));
const comm=(a,b)=>N(p.bracket(a,b)),zero=system.zero(),one=system.one();
const R=N(system.word(['R'])),K=N(system.word(['K'])),RK=mul(R,K),basis=[one,R,K,RK];
const G=[R,K,RK],c=Array.from({length:3},()=>Array.from({length:3},()=>[0,0,0]));
function setBracket(t,i,j,k,v){t[i][j][k]=v;t[j][i][k]=-v;}
setBracket(c,0,1,2,2);setBracket(c,1,2,0,-2);setBracket(c,2,0,1,2);
const sum=xs=>xs.reduce(add,zero),linear=(coefficients,values)=>sum(coefficients.map((a,k)=>scale(values[k],a)));
const delta=(g,i,x)=>comm(g[i],x);
const AtoH=(g,A)=>g.map((x,i)=>add(x,A[i]));
const connection=(g,A,i,x)=>add(delta(g,i,x),mul(A[i],x));
const F=(g,t,A,i,j)=>sub(add(sub(delta(g,i,A[j]),delta(g,j,A[i])),comm(A[i],A[j])),linear(t[i][j],A));
const compactF=(g,t,A,i,j)=>sub(comm(add(g[i],A[i]),add(g[j],A[j])),linear(t[i][j],AtoH(g,A)));
const cov=(g,A,i,x)=>add(delta(g,i,x),comm(A[i],x));
const cyclic=[[0,1,2],[1,2,0],[2,0,1]];
function bianchi(g,t,A){
  const derivative=sum(cyclic.map(([i,j,k])=>cov(g,A,i,F(g,t,A,j,k))));
  const frame=sum(cyclic.flatMap(([i,j,k])=>t[j][k].map((v,l)=>scale(F(g,t,A,i,l),v))));
  return {derivative,frame,total:add(derivative,frame)};
}
let identityCount=0;
const witnesses=[],checks=[];
function equal(a,b,label,keep=false){
  requireTrue(N(a.minus(b)).terms.size===0,'Failed identity: '+label);identityCount++;
  if(keep){
    const proof=p.proveEquality(a,b);
    requireTrue(proof.status==='EQUAL_IN_DECLARED_QUOTIENT'&&system.replay(proof.certificate),'Failed replay: '+label);
    witnesses.push({name:label,certificate:proof.certificate,replay:'REPLAY_MATCH'});
  }
}
const different=(a,b,label)=>requireTrue(N(a.minus(b)).terms.size!==0,'Counterexample failed: '+label);
function group(name,fn){const before=identityCount,details=fn()||{};checks.push({name,passed:true,exact_identities:identityCount-before,...details});}
function frameChecks(g,t){
  for(let i=0;i<3;i++)for(let j=0;j<3;j++)equal(comm(g[i],g[j]),linear(t[i][j],g),'frame bracket');
  for(let i=0;i<3;i++)for(let j=0;j<3;j++)for(let k=0;k<3;k++)for(let m=0;m<3;m++){
    let v=0;for(let l=0;l<3;l++)v+=t[j][k][l]*t[i][l][m]+t[k][i][l]*t[j][l][m]+t[i][j][l]*t[k][l][m];
    requireTrue(v===0,'Frame Jacobi constants');
  }
}
group('EMK_frame_and_derived_direction_operators',()=>{
  frameChecks(G,c);
  equal(p.bracket(R,K),scale(RK,2),'R-K bracket',true);
  equal(p.bracket(K,RK),scale(R,-2),'K-RK bracket',true);
  equal(p.bracket(RK,R),scale(K,2),'RK-R bracket',true);
  equal(R.times(R),scale(one,-1),'primitive RR',true);
  equal(K.times(K),one,'primitive KK',true);
  for(let i=0;i<3;i++)for(let j=0;j<3;j++)for(const x of basis)
    equal(sub(delta(G,i,delta(G,j,x)),delta(G,j,delta(G,i,x))),sum(c[i][j].map((v,l)=>scale(delta(G,l,x),v))),'derived direction commutator');
  return {frame_directions:3,derivative_commutator_instances:36};
});
group('native_Leibniz_on_algebra_basis',()=>{
  for(let i=0;i<3;i++)for(const a of basis)for(const b of basis)
    equal(delta(G,i,mul(a,b)),add(mul(delta(G,i,a),b),mul(a,delta(G,i,b))),'inner Leibniz');
  return {basis_product_instances:48};
});
const trials=[];
for(const a of [-1,0,1])for(const b of [-1,0,1])for(const d of [-1,0,1])
  trials.push([sum([scale(R,a),scale(K,b),scale(RK,o.IOTA.mul(d))]),
               sum([scale(RK,b),scale(R,d),scale(K,a)]),
               sum([scale(K,d),scale(RK,a),scale(one,b)])]);
group('regular_module_connection_and_curvature',()=>{
  for(let i=0;i<3;i++)for(const x of basis)for(const a of basis)
    equal(connection(G,G,i,mul(x,a)),add(mul(connection(G,G,i,x),a),mul(x,delta(G,i,a))),'right-module connection Leibniz');
  for(const A of trials)for(let i=0;i<3;i++)for(let j=0;j<3;j++){
    equal(F(G,c,A,i,j),compactF(G,c,A,i,j),'compact inner-frame curvature');
    for(const x of basis){
      const lhs=sub(sub(connection(G,A,i,connection(G,A,j,x)),connection(G,A,j,connection(G,A,i,x))),
                    sum(c[i][j].map((v,l)=>scale(connection(G,A,l,x),v))));
      equal(lhs,mul(F(G,c,A,i,j),x),'regular module curvature');
    }
  }
  const A=G;
  equal(F(G,c,A,0,1),scale(RK,4),'scaled connection curvature',true);
  return {connection_samples:27,module_curvature_instances:972,cut_turn_coefficients_retained:true};
});
const u=add(one,R),uInv=scale(sub(one,R),'1/2');
group('gauge_covariance_in_native_words',()=>{
  equal(u.times(uInv),one,'two-sided gauge inverse left',true);
  equal(uInv.times(u),one,'two-sided gauge inverse right',true);
  const gauge=A=>A.map((a,i)=>add(mul(mul(uInv,a),u),mul(uInv,delta(G,i,u))));
  for(const A of trials){
    const Ap=gauge(A);
    for(let i=0;i<3;i++){
      equal(add(G[i],Ap[i]),mul(mul(uInv,add(G[i],A[i])),u),'H gauge transform');
      for(const x of basis)
        equal(connection(G,Ap,i,x),mul(uInv,connection(G,A,i,mul(u,x))),'connection covariance');
      for(let j=0;j<3;j++)
        equal(F(G,c,Ap,i,j),mul(mul(uInv,F(G,c,A,i,j)),u),'curvature covariance');
    }
  }
  const pure=gauge([zero,zero,zero]);
  for(let i=0;i<3;i++)for(let j=0;j<3;j++)equal(F(G,c,pure,i,j),zero,'pure gauge flat');
  equal(F(G,c,pure,0,1),zero,'pure gauge flat word',true);
  return {connection_samples:27,gauge_curvature_instances:243,gauge_module_instances:324};
});
const affineG=[RK,add(K,R),one],affineC=Array.from({length:3},()=>Array.from({length:3},()=>[0,0,0]));
setBracket(affineC,0,1,1,2);
const affineA=[R,K,R];
group('Bianchi_and_noncommuting_frame_correction',()=>{
  for(const A of trials)equal(bianchi(G,c,A).total,zero,'EMK Bianchi');
  frameChecks(affineG,affineC);
  const b=bianchi(affineG,affineC,affineA);
  equal(b.derivative,scale(RK,-8),'affine derivative residue',true);
  equal(b.frame,scale(RK,8),'affine frame repair',true);
  equal(b.total,zero,'affine Bianchi closure',true);
  different(b.derivative,zero,'frame-free Bianchi');
  return {EMK_connection_samples:27,nonunimodular_frame_control:true};
});
group('mixed_channel_closure',()=>{
  const g=[zero,zero,zero],t=Array.from({length:3},()=>Array.from({length:3},()=>[0,0,0]));
  const a=[R,zero,zero],b=[zero,K,zero],total=a.map((x,i)=>add(x,b[i]));
  for(let i=0;i<3;i++)for(let j=0;j<3;j++){
    equal(F(g,t,a,i,j),zero,'single R channel flat');
    equal(F(g,t,b,i,j),zero,'single K channel flat');
    equal(F(g,t,total,i,j),sum([F(g,t,a,i,j),F(g,t,b,i,j),comm(a[i],b[j]),comm(b[i],a[j])]),'mixed channel distribution');
  }
  equal(F(g,t,total,0,1),scale(RK,2),'mixed curvature word',true);
  different(F(g,t,total,0,1),zero,'separate channels hide full curvature');
});
group('observer_descent_and_curvature_blindness',()=>{
  const words=[[],['R'],['K'],['R','K']];
  const coordinates=x=>words.map(w=>N(x).terms.get(JSON.stringify(w))||o.ZERO);
  const action=fn=>{const columns=basis.map(x=>coordinates(fn(x)));return o.matrix(words.map((_,i)=>columns.map(v=>v[i])));};
  const block=(a,b)=>o.matrix(Array.from({length:8},(_,i)=>Array.from({length:8},(_,j)=>i<4&&j<4?a[i][j]:i>=4&&j>=4?b[i-4][j-4]:o.ZERO)));
  const retained=o.matrix(Array.from({length:4},(_,i)=>Array.from({length:8},(_,j)=>j===i+4?1:0)));
  const curved=G.map((_,i)=>action(x=>connection(G,G,i,x)));
  const flat=G.map((_,i)=>action(x=>connection(G,[zero,zero,zero],i,x)));
  const joint=curved.map((x,i)=>block(x,flat[i]));
  for(let i=0;i<3;i++)requireTrue(o.equal(o.mul(retained,joint[i]),o.mul(flat[i],retained)),'Observer derivative descent');
  for(let i=0;i<3;i++)for(let j=0;j<3;j++){
    let native=o.sub(o.mul(joint[i],joint[j]),o.mul(joint[j],joint[i]));
    for(let l=0;l<3;l++)native=o.sub(native,o.scale(joint[l],c[i][j][l]));
    const expected=block(action(x=>mul(F(G,c,G,i,j),x)),o.zeros(4));
    requireTrue(o.equal(native,expected)&&o.isZero(o.mul(retained,native)),'Curvature descent/hidden copy');
    if(i===0&&j===1)requireTrue(!o.isZero(native),'Missing hidden native curvature');
  }
  const visible=x=>coordinates(x)[0];
  requireTrue(visible(R).zero()&&!visible(connection(G,G,0,R)).zero(),'Missing noninvariant observer kernel');
  return {connection_descent_instances:3,curvature_descent_instances:9,
          nonzero_native_curvature_with_flat_observation:true,noninvariant_kernel_control:true};
});
const negative={};
group('boundary_and_mutation_controls',()=>{
  different(comm(G[0],G[1]),zero,'native direction order defect');negative.noncommuting_frame_is_not_extra_connection_curvature=true;
  equal(F(G,c,[zero,zero,zero],0,1),zero,'zero connection curved-frame control',true);
  const A=G,withoutMixed=sub(sub(delta(G,0,A[1]),delta(G,1,A[0])),linear(c[0][1],A));
  different(withoutMixed,F(G,c,A,0,1),'omitted A commutator');negative.omitted_mixed_commutator_rejected=true;
  different(comm(K,RK),scale(R,2),'wrong frame sign');negative.wrong_EMK_frame_sign_rejected=true;
  different(mul(uInv,delta(G,1,u)),zero,'gauge derivative term');negative.gauge_conjugation_without_derivative_rejected=true;
  negative.frame_free_Bianchi_rejected=bianchi(affineG,affineC,affineA).derivative.terms.size!==0;
  const proof=witnesses.find(w=>w.certificate.steps.length);
  requireTrue(proof,'No nontrivial word proof');
  const damaged=JSON.parse(JSON.stringify(proof.certificate));damaged.steps[0].rule='__altered_rule__';
  try{system.replay(damaged);negative.altered_word_certificate_rejected=false;}catch{negative.altered_word_certificate_rejected=true;}
  try{checkFiles(rkf,[{...pins.canonical_engine.files[0],sha256:'0'.repeat(64)}]);negative.altered_source_pin_rejected=false;}
  catch{negative.altered_source_pin_rejected=true;}
  requireTrue(Object.values(negative).every(Boolean),'Failed mutation control');
});
const chapter=fs.readFileSync(path.join(root,'papers/emk-ugd-algebra/CONNECTION_CALCULUS.md'),'utf8');
const sections=[...chapter.matchAll(/## C([1-7])\. ([^\n]+)\n([\s\S]*?)(?=\n## |\s*$)/g)];
requireTrue(sections.length===7&&sections.every((x,i)=>Number(x[1])===i+1&&(x[3].includes('**Proof.**')||Number(x[1])===7)),'Missing chapter proof/section binding');
const result={schema:'publications.emkc1.v1',status:'PASS_EMKC1_NATIVE_CONNECTION_CALCULUS',
 source_pins_sha256:hash(pinsBytes),canonical_engine_commit:pins.canonical_engine.commit,
 presentation_sha256:system.hash(),exact_check_groups:checks.length,exact_identities:identityCount,checks,
 written_sections:sections.map(s=>({id:'C'+s[1],title:s[2],sha256:hash(Buffer.from(s[0]))})),
 symbolic_replay_count:witnesses.length,symbolic_replays:witnesses,rejected_false_alternatives:negative,
 scope:{written_general_connection_identities:true,finite_native_EMK_execution:true,
  full_Morphic_manuscript_certified:false,formal_proof_assistant_verified:false,
  physical_frame_or_metric_selected:false,physical_clock_selected:false,infinite_domain_control_proved:false,
  UGD_numeral_layer_replaced:false}};
const bytes=Buffer.from(JSON.stringify(p.stable(result),null,2)+'\n');
const resultPath=path.join(__dirname,'EMKC1_RESULT.json'),expectedPath=path.join(__dirname,'EXPECTED_EMKC1.sha256');
if(process.argv.includes('--write')){
  fs.writeFileSync(resultPath,bytes);fs.writeFileSync(expectedPath,hash(bytes)+'\n');
}else{
  requireTrue(process.argv.includes('--check'),'Use --check or --write');
  requireTrue(fs.readFileSync(resultPath).equals(bytes)&&fs.readFileSync(expectedPath,'utf8').trim()===hash(bytes),'Frozen certificate mismatch');
}
console.log(JSON.stringify({status:result.status,exact_check_groups:checks.length,exact_identities:identityCount,
 symbolic_replays:witnesses.length,negative_controls:Object.keys(negative).length,certificate_sha256:hash(bytes)}));
