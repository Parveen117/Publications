'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const R=require('./source.cjs'),{phaseEngine}=require('../native-lmd-first-mode-pulse-solver/phase.cjs'),{verify}=require('./verify.cjs');
const {f,iv,point,add,sub,neg,mul,scale,div,exp,sinh,width,absmax,overlap,contains,widen,trig,pack,productPartition}=R;
const directory=process.argv[2];assert(directory,'output directory required');
const read=p=>JSON.parse(fs.readFileSync(p)),write=(name,value)=>fs.writeFileSync(path.join(directory,name),JSON.stringify(value,null,2)+'\n');
const request=read(path.join(directory,'request.json'));
assert.deepEqual(request,read(path.join(__dirname,'request.json')),'unsupported campaign; commit a new protocol before changing it');
assert.equal(request.schema,'lmd-exact-total-mass-marginal-v1');
assert.equal(request.q_radius,'2');assert.deepEqual(request.exponential_octaves,[-16,4]);
assert.equal(request.gamma_integral_endpoint,'64');assert.equal(request.gamma_series_max_index,384);
assert.equal(request.maximum_probability_halfwidth,'1/200');
assert.deepEqual(request.partitions,[{q_bins:32,subdivisions_per_octave:2},{q_bins:64,subdivisions_per_octave:4}]);
const pins=read(path.join(__dirname,'SOURCE_PINS.json'));
for(const [p,hash] of Object.entries(pins.sha256)){
 assert.equal(crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,p))).digest('hex'),hash,'source pin:'+p);
}
const inherited=path.join(__dirname,'../native-lmd-certified-pulse-solver/results');
const source=read(path.join(__dirname,'../native-lmd-pulse-characteristic/results/source_intervals.json'));
const original=read(path.join(inherited,'solver.results.json'));
const previous=read(path.join(__dirname,'../native-lmd-first-mode-pulse-solver/results/solver.results.json'));
const previousReport=read(path.join(__dirname,'../native-lmd-first-mode-pulse-solver/results/verification.json'));
const verified=verify(),unpack=x=>iv(x.lower,x.upper),checks=verified.checks,negativeControls=verified.negatives;
const check=(name,ok)=>{assert(ok,name);checks.push(name);};
const reject=(name,ok)=>{assert(ok,name);negativeControls.push(name);};
const engine=phaseEngine(inherited,check),target=f(request.maximum_probability_halfwidth),cache=new Map(),summaries=[];
const clamp=x=>iv(x.lo.n<0n?0:x.lo,f(1).le(x.hi)?1:x.hi);
const boundary=x=>x===null?'-infinity':String(x);
const cells=rows=>rows.map(r=>({left:boundary(r.left),right:boundary(r.right),p:pack(r.p),first:pack(r.first),mass:pack(r.mass)}));

function partition(h,bins){
 const key=h+':'+bins.q_bins+':'+bins.subdivisions_per_octave;if(cache.has(key))return cache.get(key);
 const p=productPartition(h,bins);
 check('source_probability:'+key,contains(mul(p.q.totalP,p.exponential.totalP),1));
 check('source_mean:'+key,contains(mul(p.q.totalA,p.exponential.totalA),1));
 check('independent_source_second_moment:'+key,overlap(p.fullSecond,unpack(source.integrals[h])));
 check('variance_decomposition:'+key,overlap(p.residual,Object.values(p.decomposition).reduce(add,point(0))));
 check('positive_source_cells:'+key,p.rows.every(r=>r.p.lo.n>0n&&r.mass.lo.n>0n));
 summaries.push({h,...bins,cells:p.rows.length,q_probability:pack(p.q.totalP),q_mass:pack(p.q.totalA),
  exponential_probability:pack(p.exponential.totalP),exponential_mass:pack(p.exponential.totalA),
  gamma_shape:p.exponential.shape.toString(),gamma_normalizer:pack(p.exponential.normalizer),
  q_second:pack(p.q.second),exponential_second:pack(p.exponential.second),
  full_second:pack(p.fullSecond),captured_second:pack(p.second),total_variance_budget:pack(p.residual),
  variance_decomposition:Object.fromEntries(Object.entries(p.decomposition).map(([k,v])=>[k,pack(v)])),
  omitted_fourier_mode_budget:'0 under the credited full total-mass law',
  q_cells:cells(p.q.rows),exponential_cells:cells(p.exponential.rows)});
 cache.set(key,p);console.log('Prepared full-mass partition',key,'cells',p.rows.length);return p;
}

function solve(c,bins){
 const h=f(c.h),t=f(c.t),k=f(c.k),p=partition(c.h,bins);
 const damping=exp(point(h.mul(k.pow(2)).div(4).neg())),a=scale(sinh(point(h.mul(k).div(2))),t.mul(2).div(h));
 let re=point(0),im=point(0),totalP=point(0);
 const sampleIndices=new Set([0,1,Math.floor(p.rows.length/4),Math.floor(p.rows.length/2),Math.floor(3*p.rows.length/4),p.rows.length-2,p.rows.length-1]);
 const sampled=[];
 for(let index=0;index<p.rows.length;index++){
  const row=p.rows[index],angle=mul(a,row.mass),value=engine.phase(angle);
  re=add(re,mul(row.p,value.cos));im=add(im,mul(row.p,value.sin));totalP=add(totalP,row.p);
  if(sampleIndices.has(index)){
   const native=trig(angle),label=[c.id,bins.q_bins,bins.subdivisions_per_octave,index].join(':');
   check('selected_phase_cos:'+label,overlap(value.cos,native.cos));
   check('selected_phase_sin:'+label,overlap(value.sin,native.sin));
   sampled.push({index,angle:pack(angle),interpolated:{cos:pack(value.cos),sin:pack(value.sin)},independent_taylor:{cos:pack(native.cos),sin:pack(native.sin)}});
  }
 }
 check('assembled_source_probability:'+c.id+':'+bins.q_bins,contains(totalP,1));
 const compression=absmax(a).pow(2).mul(p.residual.hi).div(2);
 const real=mul(damping,widen(re,compression)),imaginary=mul(damping,widen(im,compression));
 const pc=clamp(scale(add(point(1),real),'1/2')),ps=clamp(scale(add(point(1),imaginary),'1/2'));
 const halfCos=width(pc).div(2),halfSin=width(ps).div(2),passed=halfCos.le(target)&&halfSin.le(target);
 const transfer=scale(mul(damping,mul(a,a)),'1/4');
 console.log('Evaluated',c.id,'partition',bins.q_bins,bins.subdivisions_per_octave,'tolerance',passed?'met':'not met');
 return{...bins,cells:p.rows.length,p_cos_plus:pack(pc),p_sin_plus:pack(ps),phi_cos:pack(real),phi_sin:pack(imaginary),
  phase_argument:pack(a),damping:pack(damping),probability_halfwidth:{cos:halfCos.toString(),sin:halfSin.toString()},
  tolerance_status:passed?'certified_within_requested_halfwidth':'certified_enclosure_tolerance_not_met',
  error_budget:{source_variance:pack(p.residual),source_component_error:compression.toString(),
   phase_component_error:t.zero()?'0':engine.error.toString(),
   source_probability_components:Object.fromEntries(Object.entries(p.decomposition).map(([k,v])=>[k,pack(mul(transfer,v))])),
   omitted_fourier_mode_budget:'0 conditional on Remy theorem and the audited source match',
   arithmetic_and_phase_node_enclosures:'propagated in final endpoints',
   interpretation:'conditional-variance upper certificate; not an actual-error lower bound'},
  selected_independent_phase_checks:sampled};
}

const cases=[];
for(const c of request.cases){
 const refinements=request.partitions.map(b=>solve(c,b));
 let old,oldSource;
 if(c.id==='remaining'){
  old=previous.cases.find(x=>x.h===c.h&&x.t===c.t&&x.k===c.k).refinements.at(-1);oldSource=pins.frozen_parent;
 }else if(c.id==='overlap'||c.id==='zero'){
  old=original.cases.find(x=>x.h===c.h&&x.t===c.t&&x.k===c.k).refinements.at(-1).actual;
  oldSource=pins.original_solver_commit;
 }
 if(old)for(const r of refinements){
  check('predecessor_cos_overlap:'+c.id+':'+r.q_bins,overlap(unpack(old.p_cos_plus),unpack(r.p_cos_plus)));
  check('predecessor_sin_overlap:'+c.id+':'+r.q_bins,overlap(unpack(old.p_sin_plus),unpack(r.p_sin_plus)));
 }
 cases.push({...c,baseline:old?{source_commit:oldSource,p_cos_plus:old.p_cos_plus,p_sin_plus:old.p_sin_plus,
  probability_halfwidth:old.probability_halfwidth,tolerance_status:old.tolerance_status}:null,refinements});
}
for(const h of [...new Set(request.cases.map(c=>c.h))]){
 const coarse=partition(h,request.partitions[0]),fine=partition(h,request.partitions[1]);
 check('nested_variance_improves:'+h,fine.residual.hi.le(coarse.residual.lo));
}
const positive=cases.find(c=>c.id==='remaining'),negative=cases.find(c=>c.id==='reverse'),zero=cases.find(c=>c.id==='zero');
for(let j=0;j<request.partitions.length;j++){
 const p=positive.refinements[j],n=negative.refinements[j],z=zero.refinements[j];
 check('time_reversal_cos:'+j,overlap(unpack(p.phi_cos),unpack(n.phi_cos)));
 check('time_reversal_sin:'+j,overlap(unpack(p.phi_sin),neg(unpack(n.phi_sin))));
 check('zero_pulse_cos:'+j,overlap(unpack(z.phi_cos),unpack(z.damping)));
 check('zero_pulse_sin:'+j,contains(unpack(z.phi_sin),0));
}
reject('pulse_direction_is_detected',!overlap(unpack(positive.refinements.at(-1).phi_sin),unpack(negative.refinements.at(-1).phi_sin)));
const last=positive.refinements.at(-1),I=partition(positive.h,request.partitions.at(-1)).fullSecond;
const weight=div(point(1),I),twoPointPhase=engine.phase(mul(unpack(last.phase_argument),I));
const twoPointCos=scale(add(point(1),mul(unpack(last.damping),add(sub(point(1),weight),mul(weight,twoPointPhase.cos)))),'1/2');
check('two_point_law_mean',contains(mul(weight,I),1));
check('two_point_law_second',overlap(mul(weight,mul(I,I)),I));
reject('two_moments_do_not_determine_target',!overlap(twoPointCos,unpack(last.p_cos_plus)));
const retainedOriginal=original.cases.filter(c=>c.refinements.at(-1).actual.tolerance_status==='certified_within_requested_halfwidth');
const retainedFirst=previous.cases.filter(c=>c.refinements.at(-1).tolerance_status==='certified_within_requested_halfwidth');
assert.equal(retainedOriginal.length,8);assert.equal(retainedFirst.length,3);
assert.equal(previousReport.combined_cases_within_tolerance,11);assert.equal(original.cases.length,12);
const added=last.tolerance_status==='certified_within_requested_halfwidth'?1:0;
const report={status:'PASS',checks:checks.length,negative_controls:negativeControls.length,check_names:checks,negative_control_names:negativeControls,
 target_halfwidth:target.toString(),new_cases:cases.length,
 new_cases_within_tolerance:cases.filter(c=>c.refinements.at(-1).tolerance_status==='certified_within_requested_halfwidth').length,
 retained_previous_passes:11,new_original_case_passes:added,combined_cases_within_tolerance:11+added,total_original_cases:12,
 combined_count_note:'eight original plus three first-mode predecessor certifications plus the remaining case here; four new cases including controls, not twelve new executions',
 fresh_nrqc_requests:0,reused_nrqc_requests:9,source_model_changed:false,hardware_called:false,
 new_external_dependency:request.new_external_dependency,gravitational_dynamics_derived:false,advantage_claimed:false,
 two_moment_counterexample:{probability_at_I:pack(weight),cos_plus:pack(twoPointCos),disjoint_from_full_mass_certificate:true},
 summary:cases.map(c=>{const r=c.refinements.at(-1);return{id:c.id,h:c.h,t:c.t,k:c.k,
  actual_cos_plus:[r.p_cos_plus.display_lower,r.p_cos_plus.display_upper],actual_sin_plus:[r.p_sin_plus.display_lower,r.p_sin_plus.display_upper],
  tolerance:r.tolerance_status};})};
write('solver.results.json',{request,partitions:summaries,cases});write('verification.json',report);
console.log(JSON.stringify({status:report.status,checks:report.checks,negative_controls:report.negative_controls,
 combined_passes:report.combined_cases_within_tolerance,summary:report.summary},null,2));
