'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const assert=require('node:assert/strict'),cp=require('node:child_process');
const here=__dirname,repo=path.resolve(here,'../..'),args=process.argv.slice(2);
const pos=args.indexOf('--rkf-root');
if(pos<0||!args[pos+1]||(args.includes('--write')&&args.includes('--check')))
  throw new Error('Supply --rkf-root PATH and either --check or --write');
const rkf=path.resolve(args[pos+1]),sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(here,'SOURCE_PINS.json'),'utf8'));
for(const [root,files] of [[repo,pins.publications_files],[rkf,pins.rkf_files]])
  for(const f of files)assert.equal(sha(fs.readFileSync(path.join(root,f.path))),f.sha256,'Source drift: '+f.path);
const engine=path.join(rkf,'operator_foundation');
const p=require(path.join(engine,'core/paninian_operator.cjs'));
const {Cut,F}=require(path.join(engine,'core/native_operator.cjs'));
const specification=JSON.parse(fs.readFileSync(path.join(engine,'examples/emk_job.json'),'utf8')).presentation;
const s=new p.Presentation(specification),one=s.one(),R=s.word(['R']),K=s.word(['K']),L=R.times(K);
const checks=[],proofs={};
function check(name,fn){fn();checks.push(name);}
function equal(a,b){const e=p.proveEquality(a,b);assert.equal(e.status,'EQUAL_IN_DECLARED_QUOTIENT');assert(s.replay(e.certificate));return e.certificate;}
function distinct(a,b){assert.equal(p.proveEquality(a,b).status,'DISTINCT_IN_DECLARED_QUOTIENT');}
check('native_distinction_plane',()=>{
  equal(L.times(L),one);equal(K.times(L).plus(L.times(K)),s.zero());
  const star={R:R.scale(-1),K};
  assert.equal(p.daggerGate(s,star).status,'DAGGER_DESCENDS');
  equal(p.star(L,star),L);
});
check('apertures_and_centre_direction_dependence',()=>{
  for(const [a,b] of [[1,0],[-1,0],[0,1],[0,-1],['3/5','4/5'],['4/5','-3/5']]){
    const n=K.scale(new Cut(a)).plus(L.scale(new Cut(b)));
    equal(n.times(n),one);
    const aperture=one.plus(n).scale(new Cut('1/2'));
    equal(aperture.times(aperture),aperture);
  }
  const plus=one.plus(K).scale(new Cut('1/2')),minus=one.minus(K).scale(new Cut('1/2'));
  distinct(plus,minus);const half=one.scale(new Cut('1/2'));
  distinct(half.times(half),half);
});
const P=one.plus(K).scale(new Cut('1/2')),Q=one.minus(P),PL=one.plus(L).scale(new Cut('1/2'));
check('formation_record_reconstruction_and_returning_memory',()=>{
  equal(P.plus(Q),one);equal(P.times(Q),s.zero());
  proofs.invisible=equal(P.times(R).times(P),s.zero());
  proofs.return=equal(P.times(R).times(Q).times(R).times(P),P.scale(-1));
  distinct(P.times(R).times(Q),s.zero());
  proofs.cut_order=equal(P.times(PL).minus(PL.times(P)),R.scale(new Cut('-1/2')));
});
const run=cp.spawnSync('python',[path.join(here,'controls.py')],{encoding:'utf8',timeout:90000});
assert.equal(run.status,0,run.stderr);const result=JSON.parse(run.stdout);
check('energy_derived_seam_uses_native_square_law',()=>{
  for(const row of result.native_seams){
    const a=new F(row.a),b=new F(row.b),D=K.scale(new Cut(a)).plus(L.scale(new Cut(b)));
    equal(D.times(D),one.scale(new Cut(new F(row.squared_norm))));
    assert(a.mul(a).add(b.mul(b)).eq(new F(row.squared_norm)));
  }
});
const sources=fs.readdirSync(here).filter(f=>/\.(md|py|cjs|json)$/.test(f)&&f!=='CERTIFICATE.json').sort();
const report=p.stable({protocol:'COHERENCE_FIRST_THERMO_R1',status:'PASS_EXACT_CONDITIONAL_NATIVE_CONTROLS',
  native_checks:checks,native_proofs:proofs,thermo_and_information:result,
  source_sha256:Object.fromEntries(sources.map(f=>[f,sha(fs.readFileSync(path.join(here,f)))])),
  rkf_commit:pins.rkf_commit,scope:{coherence_first_is_proposed_principle:true,constitutive_energy_declared:true,
    physical_clock_selected:false,spacetime_gravity_derived:false,infinite_Shannon_information_derived:false,
    upstream_engine_recertified:false,independent_proof_assistant:false}});
const text=JSON.stringify(report,null,2)+'\n',digest=sha(text);
if(args.includes('--write')){
  fs.writeFileSync(path.join(here,'CERTIFICATE.json'),text);
  fs.writeFileSync(path.join(here,'EXPECTED.sha256'),digest+'\n');
}else{
  assert.equal(fs.readFileSync(path.join(here,'CERTIFICATE.json'),'utf8'),text,'Certificate mismatch; no rewrite');
  assert.equal(fs.readFileSync(path.join(here,'EXPECTED.sha256'),'utf8').trim(),digest);
}
console.log(JSON.stringify({status:report.status,sha256:digest,native_checks:checks.length,
  counts:result.counts,runtime:process.version,mode:args.includes('--write')?'write':'check'}));
