'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const here=__dirname,repo=path.resolve(here,'../..'),args=process.argv.slice(2),pos=args.indexOf('--rkf-root');
if(pos<0||!args[pos+1]||(args.includes('--write')&&args.includes('--check')))throw Error('Supply --rkf-root PATH and --check or --write');
const rkf=path.resolve(args[pos+1]),sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(here,'SOURCE_PINS.json'),'utf8'));
for(const [root,files] of [[repo,pins.publications_files],[rkf,pins.rkf_files]])for(const f of files)
 assert.equal(sha(fs.readFileSync(path.join(root,f.path))),f.sha256,'Source drift: '+f.path);
const {F,Cut}=require(path.join(rkf,'operator_foundation/core/native_operator.cjs'));
const p=require(path.join(rkf,'operator_foundation/core/paninian_operator.cjs'));
const r1=require(path.join(rkf,'research/recognition_return/return_solver.cjs'));
const source=require(path.join(rkf,'research/recognition_return/r2/variable_return.cjs'));
const f=x=>new F(x),eq=(a,b)=>assert(F.of(a).eq(b)),lt=(a,b)=>F.of(a).le(b)&&!F.of(a).eq(b);
const s=r1.emk(),I=s.one(),R=s.word(['R']),K=s.word(['K']),L=K.times(R);
const P=I.plus(L).scale(new Cut('1/2')),Q=I.minus(P);
function prove(a,b){const e=p.proveEquality(a,b);assert.equal(e.status,'EQUAL_IN_DECLARED_QUOTIENT');assert(s.replay(e.certificate));return e.certificate;}
function distinct(a,b){assert.equal(p.proveEquality(a,b).status,'DISTINCT_IN_DECLARED_QUOTIENT');}
const proofs={plus:prove(P.times(P),P),minus:prove(Q.times(Q),Q),orthogonal:prove(P.times(Q),s.zero()),
 exchange:prove(Q.times(R),R.times(P))};
const counts={exact_native_inverses:0,source_observer_controls:0,pole_error_enclosures:0,third_probe_predictions:0,refusals:0};
const exactRows=[];
for(const bs of ['1/4','1/2','1','2','5'])for(const xs of ['1/4','1/2','3/4','7/8']){
 const b=f(bs),x=f(xs),a=b.add(1),z=x.div(a.sub(b.mul(x.pow(2)))),t=f(1).sub(z);
 eq(source.periodTwoResidual(a.mul(z),b.mul(z),x),0);assert(lt(0,z)&&lt(z,1));
 const Fz=I.plus(L.scale(new Cut(x))),G=P.scale(new Cut(f(1).div(f(1).add(x)))).plus(Q.scale(new Cut(f(1).div(f(1).sub(x)))));
 prove(Fz.times(G),I);prove(G.times(Fz),I);
 prove(G.times(P),P.scale(new Cut(f(1).div(f(1).add(x)))));
 prove(P.times(G).times(R).times(P),s.zero());
 prove(Q.times(G).times(R).times(P),R.times(P).scale(new Cut(f(1).div(f(1).sub(x)))));
 counts.source_observer_controls+=3;
 const D=f(1).add(b.mul(z).mul(f(1).add(x))),residue=b.mul(2).add(1);
 eq(t.div(f(1).sub(x)),D);assert(D.le(residue));assert(residue.sub(D).le(b.mul(t).mul(3)));
 exactRows.push({b,x,z,normalized_minus_response:D,pole_coefficient:residue});counts.exact_native_inverses++;
}
prove(I.plus(L).times(Q),s.zero());distinct(Q,s.zero());
const retainedInverse=P.scale(new Cut('1/2'));
prove(I.plus(L).times(retainedInverse),P);distinct(I.plus(L).times(retainedInverse),I);
const poleRows=[];
function inverseInterval(l,u){l=f(l);u=f(u);if(!l.le(u)||!f(0).le(l)||!lt(u,1))throw RangeError('Uncertified inverse denominator');return {lower:f(1).div(f(1).sub(l)),upper:f(1).div(f(1).sub(u))};}
for(const bs of ['1/4','1','2','5'])for(const ts of ['1/8','1/32','1/128']){
 const b=f(bs),a=b.add(1),t=f(ts),z=f(1).sub(t);
 const out=source.solvePeriodic({pattern:[a.mul(z),b.mul(z)],tolerance:'1/100000000000',maxCells:1024});
 assert.equal(out.status,'CERTIFIED_RECOGNITION_RESPONSE');
 const {lower:l,upper:u,midpoint:m,error:eps}=out.interval;
 const Dl=f(1).add(b.mul(z).mul(f(1).add(l))),Du=f(1).add(b.mul(z).mul(f(1).add(u)));
 const Dm=f(1).add(b.mul(z).mul(f(1).add(m))),residue=b.mul(2).add(1),bound=b.mul(t).mul(3).add(b.mul(z).mul(eps));
 eq(Du.sub(Dl),b.mul(z).mul(u.sub(l)));assert(residue.sub(Dm).abs().le(bound));
 const inv=inverseInterval(l,u);eq(inv.upper.sub(inv.lower),u.sub(l).div(f(1).sub(u).mul(f(1).sub(l))));
 // Two independent enclosures of the same normalized completed inverse overlap.
 assert(Dl.le(t.mul(inv.upper))&&t.mul(inv.lower).le(Du));
 poleRows.push({b,t,normalized_interval:{lower:Dl,upper:Du},pole_coefficient:residue,total_error_bound:bound,cells:out.cellsUsed});counts.pole_error_enclosures++;
}
const loose=source.enclose(['99/50','99/100']).interval;
assert(lt(loose.lower,1)&&f(1).le(loose.upper));assert.throws(()=>inverseInterval(loose.lower,loose.upper),/denominator/);counts.refusals++;
// Read the previously committed calibration; no third reading enters parameter reconstruction.
const prior=JSON.parse(fs.readFileSync(path.join(repo,'papers/native-return-identification/CERTIFICATE.json'),'utf8'));
const first=prior.records[0],b=f(first.expected.b),delta=f(first.expected.delta),u=f(first.readings.u),v=f(first.readings.v);
const A=f(1).sub(u.pow(2)),B=f(1).sub(v.pow(2));
eq(A.mul(B).mul(2).mul(b.pow(2)).add(A.add(B).mul(2).sub(u.mul(B)).sub(v.mul(A)).mul(b)).add(f(2).sub(u).sub(v)),0);
eq(u.div(f(1).add(b.mul(A))).sub(1),delta);eq(f(1).sub(v.div(f(1).add(b.mul(B)))),delta);
eq(b,'5/7');eq(delta,'3/4');
const z=f(1).add(delta.div(2)),a=b.add(1);eq(z,'11/8');
const prediction=source.solvePeriodic({pattern:[a.mul(z),b.mul(z)],tolerance:'1/10000000000000',maxCells:1024});
assert.equal(prediction.status,'CERTIFIED_RECOGNITION_RESPONSE');
function outward(x,up){const scale=1000000000000n,n=x.n*scale;return new F(up?(n+x.d-1n)/x.d:n/x.d,scale);}
const lower=outward(prediction.interval.lower,false),upper=outward(prediction.interval.upper,true);
const polynomial=x=>f(55).mul(x.pow(2)).add(f(56).mul(x)).sub(132);
assert(polynomial(lower).le(0)&&f(0).le(polynomial(upper)));
counts.third_probe_predictions++;
assert(!polynomial(f('6/5')).zero());assert(lt(upper,'6/5'));counts.refusals++;
const third={calibration_readings:{u,v},b,delta,z,polynomial:['55','56','-132'],lower,upper,
 pole_coefficient:b.mul(2).add(1),native_cells:prediction.cellsUsed,experimental_observation:false};
const names=fs.readdirSync(here).filter(x=>/\.(md|cjs|json)$/.test(x)&&x!=='CERTIFICATE.json').sort();
const report=p.stable(JSON.parse(JSON.stringify({protocol:'NATIVE_CRITICAL_RESPONSE_R1',status:'PASS_EXACT_NATIVE_CRITICAL_CHANNEL_CONTROLS',
 counts,proofs,exactRows,poleRows,third_probe:third,
 source_sha256:Object.fromEntries(names.map(x=>[x,sha(fs.readFileSync(path.join(here,x)))])),rkf_commit:pins.rkf_commit,
 scope:{inverse_of_boundary_return:true,inverse_identified_as_photon_propagator:false,
 b_predicted_from_primitives:false,physical_momentum_adapter:false,numerical_alpha:'NOT_DERIVED',
 experimental_third_probe:false,upstream_engine_recertified:false,proof_assistant:false}})));
const text=JSON.stringify(report,null,2)+'\n',digest=sha(text);
if(args.includes('--write')){fs.writeFileSync(path.join(here,'CERTIFICATE.json'),text);fs.writeFileSync(path.join(here,'EXPECTED.sha256'),digest+'\n');}
else{assert.equal(fs.readFileSync(path.join(here,'CERTIFICATE.json'),'utf8'),text,'Certificate mismatch; no rewrite');assert.equal(fs.readFileSync(path.join(here,'EXPECTED.sha256'),'utf8').trim(),digest);}
console.log(JSON.stringify({status:report.status,sha256:digest,counts,third_probe:third,alpha:'NOT_DERIVED',runtime:process.version,mode:args.includes('--write')?'write':'check'}));
