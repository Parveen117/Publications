'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const here=__dirname,repo=path.resolve(here,'../..'),args=process.argv.slice(2),pos=args.indexOf('--rkf-root');
if(pos<0||!args[pos+1]||(args.includes('--write')&&args.includes('--check')))throw Error('Supply --rkf-root PATH and --check or --write');
const rkf=path.resolve(args[pos+1]),sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(here,'SOURCE_PINS.json'),'utf8'));
for(const [root,files] of [[repo,pins.publications_files],[rkf,pins.rkf_files]])for(const f of files)
 assert.equal(sha(fs.readFileSync(path.join(root,f.path))),f.sha256,'Source drift: '+f.path);
const {F}=require(path.join(rkf,'operator_foundation/core/native_operator.cjs'));
const source=require(path.join(rkf,'research/recognition_return/r2/variable_return.cjs'));
const p=require(path.join(rkf,'operator_foundation/core/paninian_operator.cjs'));
const f=x=>new F(x),eq=(x,y)=>assert(F.of(x).eq(y)),lt=(x,y)=>F.of(x).le(y)&&!F.of(x).eq(y);
function polynomial(u,v){
 u=f(u);v=f(v);if(!lt(1,u)||!lt(0,v)||!lt(v,1)||!lt(u.add(v),2))throw RangeError('Inadmissible symmetric-probe readings');
 const A=f(1).sub(u.pow(2)),B=f(1).sub(v.pow(2));
 const coefficients=[f(2).sub(u).sub(v),A.add(B).mul(2).sub(u.mul(B)).sub(v.mul(A)),A.mul(B).mul(2)];
 return {coefficients,upper:f(1).div(u.pow(2).sub(1)),evaluate:b=>coefficients[2].mul(b).add(coefficients[1]).mul(b).add(coefficients[0])};
}
function recover(u,v,tolerance='1/100000000',budget=256){
 const q=polynomial(u,v),tol=f(tolerance);if(!lt(0,tol)||!Number.isSafeInteger(budget)||budget<1)throw RangeError('Invalid root budget');
 let lo=f(0),hi=q.upper;assert(lt(0,q.evaluate(lo))&&lt(q.evaluate(hi),0));
 for(let n=0;n<budget;n++){
  if(hi.sub(lo).le(tol))return {status:'CERTIFIED_PARAMETER_INTERVAL',lower:lo,upper:hi,steps:n};
  const mid=lo.add(hi).div(2),value=q.evaluate(mid);
  if(value.zero())return {status:'EXACT_PARAMETER',lower:mid,upper:mid,steps:n+1};
  if(lt(0,value))lo=mid;else hi=mid;
 }
 return {status:'REFINEMENT_BUDGET_EXHAUSTED',lower:lo,upper:hi,steps:budget};
}
function interval(u,v){
 if(!f(u.lower).le(u.upper)||!f(v.lower).le(v.upper))throw RangeError('Unordered readings');
 polynomial(u.lower,v.lower);polynomial(u.upper,v.upper);
 const low=recover(u.upper,v.upper),high=recover(u.lower,v.lower);
 if(low.status==='REFINEMENT_BUDGET_EXHAUSTED'||high.status==='REFINEMENT_BUDGET_EXHAUSTED')throw Error('Unclosed root budget');
 return {lower:low.lower,upper:high.upper};
}
const counts={exact_reconstructions:0,native_probe_enclosures:0,parameter_enclosures:0,
 gain_invariance:0,balanced_directional_profiles:0,refusals:0};
const fixtures=[['6/5','2/5','5/7','3/4'],['9/7','3/7','7/20','2/3'],['8/7','11/14','7/15','1/3'],
 ['14/11','2/3','3/25','3/8'],['4/3','4/9','3/13','5/8'],['13/11','3/5','5/8','4/7'],
 ['12/11','3/4','11/7','5/9'],['12/11','6/7','11/13','3/10']];
const records=[];
for(const [us,vs,bs,ds] of fixtures){
 const u=f(us),v=f(vs),b=f(bs),delta=f(ds),a=b.add(1);
 const q=polynomial(u,v);eq(q.evaluate(b),0);
 const rec=recover(u,v);assert(rec.status!=='REFINEMENT_BUDGET_EXHAUSTED');assert(rec.lower.le(b)&&b.le(rec.upper));
 const zp=u.div(f(1).add(b.mul(f(1).sub(u.pow(2))))),zm=v.div(f(1).add(b.mul(f(1).sub(v.pow(2)))));
 eq(zp.sub(1),delta);eq(f(1).sub(zm),delta);eq(zp.add(zm),2);
 eq(source.periodTwoResidual(a.mul(zp),b.mul(zp),u),0);
 eq(source.periodTwoResidual(a.mul(zm),b.mul(zm),v),0);
 const positive=source.solvePeriodic({pattern:[a.mul(zp),b.mul(zp)],tolerance:'1/100000000000',maxCells:512});
 const negative=source.solvePeriodic({pattern:[a.mul(zm),b.mul(zm)],tolerance:'1/100000000000',maxCells:512});
 for(const [out,x] of [[positive,u],[negative,v]]){assert.equal(out.status,'CERTIFIED_RECOGNITION_RESPONSE');assert(out.interval.lower.le(x)&&x.le(out.interval.upper));counts.native_probe_enclosures++;}
 const inferred=interval(positive.interval,negative.interval);assert(inferred.lower.le(b)&&b.le(inferred.upper));
 assert(inferred.upper.sub(inferred.lower).le('1/1000000'));counts.parameter_enclosures++;
 const slope=f(1).div(b.mul(2).add(1)),second=b.mul(2).mul(b.mul(4).add(3)).neg().div(b.mul(2).add(1).pow(3));
 const shape=second.neg().div(slope.pow(2));eq(shape,f(2).div(slope).sub(1).sub(slope));
 eq(b.pow(2).mul(8).add(f(6).sub(shape.mul(2)).mul(b)).sub(shape),0);
 for(const [A,g] of [['2','3'],['1/2','-2'],['3/5','1/4']]){
  eq(f(A).mul(f(A).mul(f(g).pow(2)).mul(second)).neg().div(f(A).mul(g).mul(slope).pow(2)),shape);
  eq(f(A).mul(u).div(A),u);eq(f(A).mul(v).div(A),v);counts.gain_invariance++;
 }
 records.push({readings:{u,v},reconstructed:rec,expected:{b,a,delta},native_enclosure_parameter:inferred,
 positiveCells:positive.cellsUsed,negativeCells:negative.cellsUsed,shape});counts.exact_reconstructions++;
}
for(const ts of ['2','3','4','5']){
 const t=f(ts),openingA=t.add(f(1).div(t)).div(2),openingB=t.sub(f(1).div(t)).div(2);
 const a=openingA.pow(2),b=openingB.pow(2);eq(a.sub(b),1);
 const bonds=[openingA,openingA,openingB,openingB].map(x=>({opening:x,closing:x}));
 const cells=source.pairedCellsFromBonds(bonds);eq(cells[0],a);eq(cells[1],b);
 assert.equal(source.synthesize(1,b).exactNativeCut,true);assert(!source.step(b,1).eq(1));counts.balanced_directional_profiles++;
}
for(const [u,v] of [[1,'1/2'],['6/5',0],['6/5',1],['3/2','1/2'],['3/2','3/4']]){assert.throws(()=>recover(u,v),/Inadmissible/);counts.refusals++;}
assert.throws(()=>interval({lower:f('6/5'),upper:f('3/2')},{lower:f('2/5'),upper:f('3/4')}),/Inadmissible/);counts.refusals++;
assert.equal(recover('6/5','2/5','1/100000000000000000',1).status,'REFINEMENT_BUDGET_EXHAUSTED');counts.refusals++;
assert.throws(()=>recover('6/5','2/5',0),/budget/);counts.refusals++;
const names=fs.readdirSync(here).filter(x=>/\.(md|cjs|json)$/.test(x)&&x!=='CERTIFICATE.json').sort();
const report=p.stable(JSON.parse(JSON.stringify({protocol:'NATIVE_RETURN_IDENTIFICATION_R1',status:'PASS_EXACT_IDENTIFICATION_AND_NATIVE_APERTURE_CONTROLS',
 counts,records,source_sha256:Object.fromEntries(names.map(x=>[x,sha(fs.readFileSync(path.join(here,x)))])),rkf_commit:pins.rkf_commit,
 scope:{parameter_identified_from_supplied_readings:true,parameter_predicted_from_primitives:false,physical_probe_selected:false,
 numerical_alpha:'NOT_DERIVED',experimental_data_used:false,upstream_engine_recertified:false,proof_assistant:false}})));
const text=JSON.stringify(report,null,2)+'\n',digest=sha(text);
if(args.includes('--write')){fs.writeFileSync(path.join(here,'CERTIFICATE.json'),text);fs.writeFileSync(path.join(here,'EXPECTED.sha256'),digest+'\n');}
else{assert.equal(fs.readFileSync(path.join(here,'CERTIFICATE.json'),'utf8'),text,'Certificate mismatch; no rewrite');assert.equal(fs.readFileSync(path.join(here,'EXPECTED.sha256'),'utf8').trim(),digest);}
console.log(JSON.stringify({status:report.status,sha256:digest,counts,alpha:report.scope.numerical_alpha,mode:args.includes('--write')?'write':'check',runtime:process.version}));
