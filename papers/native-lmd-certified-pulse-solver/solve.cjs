'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const R=require('./interval.cjs');
const {f,iv,point,add,neg,sub,mul,scale,div,mid,width,absmax,widen,overlap,contains,exp,sinh,PI,normalCDF,trig,decimal,pack,floor}=R;
const dir=process.argv[2],request=JSON.parse(fs.readFileSync(path.join(dir,'solver.request.json')));
assert.equal(request.schema,'lmd-certified-uniform-pulse-v1');
assert.deepEqual(Object.keys(request).sort(),['schema','h','times','k','central_bins','central_radius','maximum_probability_halfwidth','include_naive_comparison'].sort());
assert(request.include_naive_comparison===true);
assert(request.h.length>0&&request.h.every(h=>['1/8','1/4','3/8'].includes(h)));
assert(request.times.length>0&&request.times.every(t=>f(t).abs().le(1)));
assert(f(request.k).abs().le(4));assert.equal(request.central_radius,'2');
assert(request.central_bins.length>0&&request.central_bins.every(n=>[32,64,128].includes(n)));
const target=f(request.maximum_probability_halfwidth);assert(target.n>0n);
const source=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-lmd-pulse-characteristic/results/source_intervals.json')));
const checks=[],negatives=[],check=(name,v)=>{assert(v,name);checks.push(name);},reject=(name,v)=>{assert(v,name);negatives.push(name);};
const unpack=x=>iv(x.lower,x.upper),nodes=new Map();
for(let block=0;block<9;block++){
 const report=JSON.parse(fs.readFileSync(path.join(dir,`phase_${block}.result.json`)));
 assert.equal(report.status,'completed');assert.equal(report.backend,'cyclotomic64');
 for(const answer of report.answers){
  assert.equal(answer.algebraically_exact,true);assert.equal(answer.status,'completed');
  const [,label,axis]=answer.id.match(/^phase_(\d+)_(X|Y)$/),j=Number(label);
  const display=answer.value.decimal_enclosure,value=iv(decimal(display.lower),decimal(display.upper));
  if(!nodes.has(j))nodes.set(j,{});nodes.get(j)[axis==='X'?'cos':'sin']=value;
 }
}
assert.equal(nodes.size,17);
for(const [j,node] of nodes){
 const independent=trig(scale(PI,f(j).div(32)));
 check('NRQC_vs_native_cos:'+j,overlap(node.cos,independent.cos));
 check('NRQC_vs_native_sin:'+j,overlap(node.sin,independent.sin));
 check('NRQC_phase_normalization:'+j,contains(add(mul(node.cos,node.cos),mul(node.sin,node.sin)),1));
}
function node(j){
 j=((j%64)+64)%64;
 if(j<=16)return nodes.get(j);
 if(j<32){const a=nodes.get(32-j);return{cos:neg(a.cos),sin:a.sin};}
 if(j<48){const a=nodes.get(j-32);return{cos:neg(a.cos),sin:neg(a.sin)};}
 const a=nodes.get(64-j);return{cos:a.cos,sin:neg(a.sin)};
}
const phaseStep=scale(PI,'1/32'),phaseError=absmax(mul(mul(phaseStep,phaseStep),mul(phaseStep,phaseStep))).div(24);
function phase(a){
 if(a.lo.zero()&&a.hi.zero())return{cos:point(1),sin:point(0)};
 const z=div(a,phaseStep),j=Number(floor(mid(z))),u=sub(z,point(j));
 assert(Number.isSafeInteger(j)&&f('-1/100').le(u.lo)&&u.hi.le('101/100'));
 const um1=sub(u,point(1)),um2=sub(u,point(2)),up1=add(u,point(1));
 const weights=[scale(mul(mul(u,um1),um2),'-1/6'),scale(mul(mul(up1,um1),um2),'1/2'),
  scale(mul(mul(up1,u),um2),'-1/2'),scale(mul(mul(up1,u),um1),'1/6')];
 let cos=point(0),sin=point(0);
 for(let i=0;i<4;i++){const v=node(j+i-1);cos=add(cos,mul(weights[i],v.cos));sin=add(sin,mul(weights[i],v.sin));}
 return{cos:widen(cos,phaseError),sin:widen(sin,phaseError)};
}
const partitionCache=new Map();
function partition(gs,bins){
 const key=gs+':'+bins;if(partitionCache.has(key))return partitionCache.get(key);
 const g=f(gs),radius=f(request.central_radius),endpoints=[null];
 for(let j=0;j<=bins;j++)endpoints.push(radius.neg().add(radius.mul(2*j).div(bins)));endpoints.push(Infinity);
 const rows=[];let totalP=point(0),totalA=point(0),second=point(0);
 for(let j=0;j<endpoints.length-1;j++){
  const a=endpoints[j],b=endpoints[j+1],shift=x=>x===null||x===Infinity?x:x.sub(g);
  const p=sub(normalCDF(g,b),normalCDF(g,a));assert(p.lo.n>0n);
  const first=sub(normalCDF(g,shift(b)),normalCDF(g,shift(a)));assert(first.lo.n>0n);
  const mass=div(first,p);rows.push({p,first,mass});totalP=add(totalP,p);totalA=add(totalA,first);
  second=add(second,div(mul(first,first),p));
 }
 const I=unpack(source.integrals[gs]);assert(I);
 const residual=sub(I,second);assert(residual.lo.n>0n);
 check('partition_probability:'+key,contains(totalP,1));check('partition_mass_mean:'+key,contains(totalA,1));
 const baseResidual=sub(I,exp(point(g)));
 check('partition_residual_budget:'+key,!residual.hi.le(0)&&baseResidual.lo.le(residual.hi));
 const result={gamma:gs,bins,rows,second,residual,totalP};partitionCache.set(key,result);return result;
}
function solve(gs,a,damping,bins){
 const p=partition(gs,bins);let re=point(0),im=point(0),nativeRe=point(0),nativeIm=point(0);
 for(const row of p.rows){
  const angle=mul(a,row.mass),value=phase(angle),baseline=trig(angle);
  re=add(re,mul(row.p,value.cos));im=add(im,mul(row.p,value.sin));
  nativeRe=add(nativeRe,mul(row.p,baseline.cos));nativeIm=add(nativeIm,mul(row.p,baseline.sin));
 }
 const label=gs+':'+mid(a)+':'+bins;
 check('assembled_phase_cos:'+label,overlap(re,nativeRe));check('assembled_phase_sin:'+label,overlap(im,nativeIm));
 const compression=absmax(a).pow(2).mul(p.residual.hi).div(2);
 const actualRe=mul(damping,widen(re,compression)),actualIm=mul(damping,widen(im,compression));
 const clamp=x=>iv(x.lo.le(0)?0:x.lo,f(1).le(x.hi)?1:x.hi);
 const pc=clamp(scale(add(point(1),actualRe),'1/2')),ps=clamp(scale(add(point(1),actualIm),'1/2'));
 const certified=width(pc).div(2).le(target)&&width(ps).div(2).le(target);
 return{gamma:gs,bins,phi_cos:pack(actualRe),phi_sin:pack(actualIm),p_cos_plus:pack(pc),p_sin_plus:pack(ps),
  tolerance_status:certified?'certified_within_requested_halfwidth':'certified_enclosure_tolerance_not_met',
  probability_halfwidth:{cos:width(pc).div(2).toString(),sin:width(ps).div(2).toString()},
  error_budget:{source_and_partition_squared_L2:pack(p.residual),conditional_expectation_component_error:compression.toString(),
   interpolated_phase_component_error:phaseError.toString(),native_arithmetic_and_output_intervals:'included throughout'},
  discrete_native_baseline:{cos:pack(mul(damping,nativeRe)),sin:pack(mul(damping,nativeIm))}};
}
const cases=[],k=f(request.k);
for(const hs of request.h){
 const h=f(hs),damping=exp(point(h.mul(k.pow(2)).div(4).neg()));
 for(const ts of request.times){
  const t=f(ts),a=scale(sinh(point(h.mul(k).div(2))),t.mul(2).div(h)),b=point(t.mul(k));
  const refinements=[];
  for(const bins of request.central_bins){
   const actual=solve(hs,a,damping,bins),naive=solve(h.mul(2).toString(),b,damping,bins);
   const cosGap=sub(unpack(actual.p_cos_plus),unpack(naive.p_cos_plus)),sinGap=sub(unpack(actual.p_sin_plus),unpack(naive.p_sin_plus));
   const separated=x=>x.lo.n>0n||x.hi.n<0n;
   refinements.push({bins,actual,naive,probability_difference:{cos:pack(cosGap),sin:pack(sinGap)},
    separation:{cos:separated(cosGap),sin:separated(sinGap)}});
   if(t.zero()){
    check('zero_pulse_cos:'+hs+':'+bins,overlap(unpack(actual.p_cos_plus),scale(add(point(1),damping),'1/2')));
    check('zero_pulse_sin:'+hs+':'+bins,contains(unpack(actual.p_sin_plus),'1/2'));
   }
  }
  cases.push({h:hs,t:ts,k:request.k,refinements});
 }
}
for(const gs of new Set(request.h.flatMap(h=>[h,f(h).mul(2).toString()]))) {
 for(let j=1;j<request.central_bins.length;j++){
  const coarse=partition(gs,request.central_bins[j-1]),fine=partition(gs,request.central_bins[j]);
  check('nested_partition_progress:'+gs+':'+fine.bins,fine.residual.hi.le(coarse.residual.lo));
 }
}
for(const hs of request.h){
 const pos=cases.find(c=>c.h===hs&&c.t==='1/4'),negative=cases.find(c=>c.h===hs&&c.t==='-1/4');
 if(pos&&negative){
  for(let j=0;j<pos.refinements.length;j++){
   const p=pos.refinements[j].actual,n=negative.refinements[j].actual;
   check('time_reversal_cos:'+hs+':'+j,overlap(unpack(p.phi_cos),unpack(n.phi_cos)));
   check('time_reversal_sin:'+hs+':'+j,overlap(unpack(p.phi_sin),neg(unpack(n.phi_sin))));
  }
 }
}
const witness=cases.find(c=>c.h==='1/8'&&c.t==='1/4');
if(witness){
 const r=witness.refinements.at(-1);
 reject('naive_classical_density_shift',r.separation.cos||r.separation.sin);
 reject('reversing_pulse_sign_changes_sine',!overlap(unpack(r.actual.phi_sin),neg(unpack(r.actual.phi_sin))));
 const baseline=unpack(r.actual.discrete_native_baseline.cos),correct=unpack(r.actual.phi_cos);
 const h=f('1/8'),damping=exp(point(h.neg()));
 reject('omitting_detector_overlap_factor',!overlap(correct,div(baseline,damping)));
}
reject('phase_orientation_detected',!overlap(nodes.get(8).sin,neg(nodes.get(8).sin)));
const partitionSummary=[...partitionCache.values()].map(p=>({gamma:p.gamma,bins:p.bins,cells_including_tails:p.rows.length,
 mass_second_moment:pack(p.second),source_and_partition_squared_L2:pack(p.residual)}));
const summary=cases.map(c=>{const r=c.refinements.at(-1);return{h:c.h,t:c.t,k:c.k,
 actual_cos_plus:[r.actual.p_cos_plus.display_lower,r.actual.p_cos_plus.display_upper],
 actual_sin_plus:[r.actual.p_sin_plus.display_lower,r.actual.p_sin_plus.display_upper],
 naive_cos_plus:[r.naive.p_cos_plus.display_lower,r.naive.p_cos_plus.display_upper],
 separation:r.separation,tolerance:r.actual.tolerance_status};});
const report={status:'PASS',checks:checks.length,negative_controls:negatives.length,check_names:checks,negative_control_names:negatives,
 scope:'Certified continuum uniform-vacuum marginal expectations conditional on the frozen LMD/common-parent theorem. Native conditional-expectation solver with exact NRQC phase-node queries.',
 requested_probability_halfwidth:target.toString(),cases_certified_within_tolerance:summary.filter(x=>x.tolerance==='certified_within_requested_halfwidth').length,
 total_cases:summary.length,summary,phase_error_bound:phaseError.toString(),
 hardware_called:false,full_field_circuit_executed:false,advantage_claimed:false};
fs.writeFileSync(path.join(dir,'solver.results.json'),JSON.stringify({request,partitions:partitionSummary,cases},null,2)+'\n');
fs.writeFileSync(path.join(dir,'verification.json'),JSON.stringify(report,null,2)+'\n');
process.stdout.write(JSON.stringify({status:report.status,checks:checks.length,negative_controls:negatives.length,summary},null,2)+'\n');
