'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const R=require('./source.cjs'),{phaseEngine}=require('./phase.cjs'),{verify}=require('./verify.cjs');
const {f,iv,point,add,sub,neg,mul,scale,exp,sinh,width,absmax,overlap,contains,widen,trig,pack,productPartition}=R;
const directory=process.argv[2];assert(directory,'output directory required');
const read=p=>JSON.parse(fs.readFileSync(p)),write=(name,value)=>fs.writeFileSync(path.join(directory,name),JSON.stringify(value,null,2)+'\n');
const request=read(path.join(directory,'request.json'));
assert.deepEqual(request,read(path.join(__dirname,'request.json')),'unsupported campaign; commit a new protocol before changing it');
assert.equal(request.schema,'lmd-first-spatial-mode-pulse-v1');
assert.equal(request.q_radius,'2');assert.equal(request.radial_endpoint,'16');
assert.equal(request.poisson_max_index,32);assert.equal(request.bessel_max_index,32);
assert.equal(request.maximum_probability_halfwidth,'1/200');
assert.deepEqual(request.partitions,[{q_bins:64,radial_bins:16},{q_bins:128,radial_bins:32}]);
const pins=read(path.join(__dirname,'SOURCE_PINS.json'));
for(const [p,hash] of Object.entries(pins.sha256)){
 assert.equal(crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,p))).digest('hex'),hash,'source pin:'+p);
}
const inherited=path.join(__dirname,'../native-lmd-certified-pulse-solver/results');
const source=read(path.join(__dirname,'../native-lmd-pulse-characteristic/results/source_intervals.json'));
const baseline=read(path.join(inherited,'solver.results.json'));
const verified=verify(),unpack=x=>iv(x.lower,x.upper),checks=verified.checks,negativeControls=verified.negatives;
const check=(name,ok)=>{assert(ok,name);checks.push(name);};
const reject=(name,ok)=>{assert(ok,name);negativeControls.push(name);};
const engine=phaseEngine(inherited,check),target=f(request.maximum_probability_halfwidth),cache=new Map();
const rawPartitionSummaries=[];

function partition(h,bins){
 const key=h+':'+bins.q_bins+':'+bins.radial_bins;
 if(cache.has(key))return cache.get(key);
 const p=productPartition(h,bins,unpack(source.integrals[h]));
 check('source_probability:'+key,contains(mul(p.q.totalP,p.radial.totalP),1));
 check('source_mean:'+key,contains(mul(p.q.totalA,p.radial.totalA),1));
 check('improves_same_q_partition:'+key,f(1).le(p.radial.second.lo));
 check('variance_decomposition:'+key,overlap(p.residual,Object.values(p.decomposition).reduce(add,point(0))));
 check('finite_above_ideal_floor:'+key,p.floors.first.hi.le(p.residual.lo));
 const oldN1=source.cutoff_moments.find(c=>c.h===h).rows.find(c=>c.N===1);
 check('Bessel_vs_frozen_Laurent_moment:'+key,overlap(p.idealSecond,unpack(oldN1.I_N)));
 rawPartitionSummaries.push({h,...bins,cells:p.rows.length,
   q_probability:pack(p.q.totalP),q_mass:pack(p.q.totalA),
   radial_probability:pack(p.radial.totalP),radial_mass:pack(p.radial.totalA),
   q_second:pack(p.q.second),radial_second:pack(p.radial.second),
   total_second:pack(p.second),total_variance_budget:pack(p.residual),
   zero_mode_ideal_floor:pack(p.floors.zero),first_mode_ideal_floor:pack(p.floors.first),
   variance_decomposition:Object.fromEntries(Object.entries(p.decomposition).map(([k,v])=>[k,pack(v)])),
   poisson_remainder:p.radial.poissonRemainder.toString(),
   q_cells:p.q.rows.map(r=>({left:r.left===null?'-infinity':String(r.left),right:String(r.right),p:pack(r.p),first:pack(r.first),mass:pack(r.mass)})),
   radial_cells:p.radial.rows.map(r=>({left:String(r.left),right:String(r.right),p:pack(r.p),first:pack(r.first),mass:pack(r.mass)}))});
 cache.set(key,p);console.log('Prepared exact source partition',key,'cells',p.rows.length);return p;
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
   const native=trig(angle),label=[c.h,c.t,bins.q_bins,bins.radial_bins,index].join(':');
   check('selected_phase_cos:'+label,overlap(value.cos,native.cos));
   check('selected_phase_sin:'+label,overlap(value.sin,native.sin));
   sampled.push({index,angle:pack(angle),interpolated:{cos:pack(value.cos),sin:pack(value.sin)},independent_taylor:{cos:pack(native.cos),sin:pack(native.sin)}});
  }
 }
 check('assembled_source_probability:'+c.h+':'+c.t+':'+bins.q_bins,contains(totalP,1));
 const compression=absmax(a).pow(2).mul(p.residual.hi).div(2);
 const real=mul(damping,widen(re,compression)),imaginary=mul(damping,widen(im,compression));
 const clamp=x=>iv(x.lo.n<0n?0:x.lo,f(1).le(x.hi)?1:x.hi);
 const pc=clamp(scale(add(point(1),real),'1/2')),ps=clamp(scale(add(point(1),imaginary),'1/2'));
 const halfCos=width(pc).div(2),halfSin=width(ps).div(2);
 const passed=halfCos.le(target)&&halfSin.le(target);
 const transfer=scale(mul(damping,mul(a,a)),'1/4');
 const zeroFloor=mul(transfer,p.floors.zero),firstFloor=mul(transfer,p.floors.first);
 console.log('Evaluated',c.h,c.t,'partition',bins.q_bins,bins.radial_bins,'tolerance',passed?'met':'not met');
 return{...bins,cells:p.rows.length,p_cos_plus:pack(pc),p_sin_plus:pack(ps),
  phi_cos:pack(real),phi_sin:pack(imaginary),phase_argument:pack(a),damping:pack(damping),
  probability_halfwidth:{cos:halfCos.toString(),sin:halfSin.toString()},
  tolerance_status:passed?'certified_within_requested_halfwidth':'certified_enclosure_tolerance_not_met',
  error_budget:{source_variance:pack(p.residual),source_component_error:compression.toString(),
    phase_component_error:engine.error.toString(),
    zero_mode_ideal_probability_budget:pack(zeroFloor),first_mode_ideal_probability_budget:pack(firstFloor),
    finer_first_mode_cells_can_remove_this_budget_obstacle:firstFloor.hi.le(target),
    interpretation:'ideal probability budgets are terms in this variance-based certificate, not lower bounds on actual probability error'},
  selected_independent_phase_checks:sampled};
}

const cases=[];
for(const c of request.cases){
 const old=baseline.cases.find(x=>x.h===c.h&&x.t===c.t&&x.k===c.k)?.refinements.at(-1)?.actual;
 assert(old&&old.bins===128&&old.tolerance_status==='certified_enclosure_tolerance_not_met');
 const refinements=request.partitions.map(bins=>solve(c,bins));
 for(const value of refinements){
  check('old_new_cos_overlap:'+c.h+':'+c.t+':'+value.q_bins,overlap(unpack(old.p_cos_plus),unpack(value.p_cos_plus)));
  check('old_new_sin_overlap:'+c.h+':'+c.t+':'+value.q_bins,overlap(unpack(old.p_sin_plus),unpack(value.p_sin_plus)));
 }
 cases.push({...c,baseline:{source_commit:pins.frozen_parent,bins:old.bins,p_cos_plus:old.p_cos_plus,p_sin_plus:old.p_sin_plus,probability_halfwidth:old.probability_halfwidth,tolerance_status:old.tolerance_status},refinements});
}
for(const h of [...new Set(request.cases.map(c=>c.h))]){
 const coarse=partition(h,request.partitions[0]),fine=partition(h,request.partitions[1]);
 check('nested_product_variance_improves:'+h,fine.residual.hi.le(coarse.residual.lo));
 reject('zero_mode_budget_is_not_first_mode_budget:'+h,fine.floors.first.hi.le(fine.floors.zero.lo));
}
const positive=cases.find(c=>c.h==='3/8'&&c.t==='1/4'),negative=cases.find(c=>c.h==='3/8'&&c.t==='-1/4');
for(let j=0;j<request.partitions.length;j++){
 check('time_reversal_cos:'+j,overlap(unpack(positive.refinements[j].phi_cos),unpack(negative.refinements[j].phi_cos)));
 check('time_reversal_sin:'+j,overlap(unpack(positive.refinements[j].phi_sin),neg(unpack(negative.refinements[j].phi_sin))));
}
reject('pulse_direction_is_detected',!overlap(unpack(positive.refinements.at(-1).phi_sin),unpack(negative.refinements.at(-1).phi_sin)));
const previousPasses=baseline.cases.filter(c=>c.refinements.at(-1).actual.tolerance_status==='certified_within_requested_halfwidth');
assert.equal(previousPasses.length,8);assert.equal(baseline.cases.length,12);
const newPasses=cases.filter(c=>c.refinements.at(-1).tolerance_status==='certified_within_requested_halfwidth').length;
const report={status:'PASS',checks:checks.length,negative_controls:negativeControls.length,check_names:checks,negative_control_names:negativeControls,
 target_halfwidth:target.toString(),new_cases:cases.length,new_cases_within_tolerance:newPasses,
 retained_previous_passes:previousPasses.length,combined_cases_within_tolerance:previousPasses.length+newPasses,total_original_cases:12,
 combined_count_note:'eight unchanged predecessor certifications plus this four-case refinement; not twelve new executions',
 fresh_nrqc_requests:0,reused_nrqc_requests:9,source_model_changed:false,hardware_called:false,
 gravitational_dynamics_derived:false,advantage_claimed:false,
 summary:cases.map(c=>{const r=c.refinements.at(-1);return{h:c.h,t:c.t,k:c.k,
  actual_cos_plus:[r.p_cos_plus.display_lower,r.p_cos_plus.display_upper],actual_sin_plus:[r.p_sin_plus.display_lower,r.p_sin_plus.display_upper],
  tolerance:r.tolerance_status,ideal_first_mode_probability_budget:r.error_budget.first_mode_ideal_probability_budget};})};
write('solver.results.json',{request,partitions:rawPartitionSummaries,cases});write('verification.json',report);
console.log(JSON.stringify({status:report.status,checks:report.checks,negative_controls:report.negative_controls,
 new_passes:newPasses,combined_passes:report.combined_cases_within_tolerance,summary:report.summary},null,2));
