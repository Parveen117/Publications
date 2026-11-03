'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const assert=require('node:assert/strict'),cp=require('node:child_process');
const here=__dirname,repo=path.resolve(here,'../..'),args=process.argv.slice(2),pos=args.indexOf('--rkf-root');
if(pos<0||!args[pos+1]||(args.includes('--write')&&args.includes('--check')))throw Error('Supply --rkf-root PATH and --check or --write');
const rkf=path.resolve(args[pos+1]),sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(here,'SOURCE_PINS.json'),'utf8'));
for(const [root,files] of [[repo,pins.publications_files],[rkf,pins.rkf_files]])
  for(const f of files)assert.equal(sha(fs.readFileSync(path.join(root,f.path))),f.sha256,'Source drift: '+f.path);
const p=require(path.join(rkf,'operator_foundation/core/paninian_operator.cjs'));
const {Cut,F}=require(path.join(rkf,'operator_foundation/core/native_operator.cjs'));
const specification=JSON.parse(fs.readFileSync(path.join(rkf,'operator_foundation/examples/emk_job.json'),'utf8')).presentation;
const s=new p.Presentation(specification),one=s.one(),R=s.word(['R']),K=s.word(['K']),L=R.times(K);
const star={R:R.scale(-1),K};assert.equal(p.daggerGate(s,star).status,'DAGGER_DESCENDS');
const dag=a=>p.star(a,star),C=x=>new Cut(String(x));
const equal=(a,b)=>{const e=p.proveEquality(a,b);assert.equal(e.status,'EQUAL_IN_DECLARED_QUOTIENT');assert(s.replay(e.certificate));return e.certificate;};
const distinct=(a,b)=>assert.equal(p.proveEquality(a,b).status,'DISTINCT_IN_DECLARED_QUOTIENT');
const aperture=N=>one.plus(N).scale(C('1/2')),P=aperture(K),Q=one.minus(P);
const counts={finite_transports:0,moving_memory_ledgers:0,transport_cocycles:0,native_generators:0};
const proofs={};
const rotorPairs=[['1','0'],['0','1'],['-1','0'],['3/5','4/5'],['4/5','-3/5'],['-3/5','-4/5']];
const rotors=rotorPairs.map(([a,b])=>one.scale(C(a)).plus(R.scale(C(b))));
for(let i=0;i<rotors.length;i++){
  const U=rotors[i],N=U.times(K).times(dag(U)),P2=aperture(N),Q2=one.minus(P2);
  equal(dag(U).times(U),one);equal(N.times(N),one);
  equal(P.times(P2).times(P),P.scale(C(new F(rotorPairs[i][0]).mul(new F(rotorPairs[i][0])))));
  equal(P.times(Q2).times(P),P.scale(C(new F(rotorPairs[i][1]).mul(new F(rotorPairs[i][1])))));
  equal(P2.times(U).times(Q),s.zero());
  equal(P2.times(R),P2.times(R).times(P).plus(P2.times(R).times(Q)));
  counts.moving_memory_ledgers++;
  for(const z of ['1/2','2/3','2']){
    const z2=new F(z).mul(new F(z)),B=U.scale(C(z)),D1=K.scale(3),D2=N.scale(C(z2.mul(new F(3))));
    equal(B.times(D1).times(dag(B)),D2);equal(dag(B).times(B),one.scale(C(z2)));
    equal(D2.times(D2),one.scale(C(z2.mul(z2).mul(new F(9)))));
    if(z!=='1')distinct(dag(B).times(B),one);
    counts.finite_transports++;
  }
}
for(const U of rotors.slice(0,4))for(const V of rotors.slice(0,4)){
  const A=U.scale(C('1/2')),B=V.scale(C('2/3')),AB=B.times(A);
  equal(B.times(A.times(K).times(dag(A))).times(dag(B)),AB.times(K).times(dag(AB)));
  equal(dag(AB).times(AB),one.scale(C('1/9')));counts.transport_cocycles++;
}
proofs.return=equal(P.times(R).times(Q).times(R).times(P),P.scale(-1));
proofs.full_turn=equal(R.times(R),one.scale(-1));
equal(R.times(R).times(K).times(dag(R.times(R))),K);distinct(R.times(R),one);
// All-axis unitary contraction refusal: a squared amplitude of 1 cannot become 1/16.
distinct(K.times(K),K.scale(C('1/4')).times(K.scale(C('1/4'))));
const run=cp.spawnSync('python',[path.join(here,'controls.py')],{encoding:'utf8',timeout:30000});
assert.equal(run.status,0,run.stderr);const result=JSON.parse(run.stdout);
for(const row of result.native_generator_rows){
  const D=K.scale(C(row.s)).plus(L.scale(C(row.v))),dD=K.scale(C(row.ds)).plus(L.scale(C(row.dv)));
  const G=one.scale(C(row.rho)).plus(R.scale(C(row.omega))).scale(C('1/2'));
  equal(G.times(D).plus(D.times(dag(G))),dD);
  if(row.rho!=='0')distinct(G.times(D).minus(D.times(G)),dD); // scalar loss cannot be a commutator
  counts.native_generators++;
}
delete result.native_generator_rows;
const sources=fs.readdirSync(here).filter(f=>/\.(md|py|cjs|json)$/.test(f)&&f!=='CERTIFICATE.json').sort();
const report=p.stable({protocol:'NATIVE_CUT_ENERGY_TRANSPORT_R1',status:'PASS_EXACT_CONDITIONAL_NATIVE_CONTROLS',
  native_counts:counts,proofs,chart_controls:result,rkf_commit:pins.rkf_commit,
  source_sha256:Object.fromEntries(sources.map(f=>[f,sha(fs.readFileSync(path.join(here,f)))])),
  scope:{response_attenuation_is_process_postulate:true,internal_clock_derived:true,
    physical_seconds_derived:false,spacetime_gravity_derived:false,unique_fundamental_action_derived:false,
    upstream_engine_recertified:false,proof_assistant:false}});
const text=JSON.stringify(report,null,2)+'\n',digest=sha(text);
if(args.includes('--write')){fs.writeFileSync(path.join(here,'CERTIFICATE.json'),text);fs.writeFileSync(path.join(here,'EXPECTED.sha256'),digest+'\n');}
else{assert.equal(fs.readFileSync(path.join(here,'CERTIFICATE.json'),'utf8'),text,'Certificate mismatch; no rewrite');assert.equal(fs.readFileSync(path.join(here,'EXPECTED.sha256'),'utf8').trim(),digest);}
console.log(JSON.stringify({status:report.status,sha256:digest,native_counts:counts,chart_counts:result.counts,runtime:process.version,mode:args.includes('--write')?'write':'check'}));
