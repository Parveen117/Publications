'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const R=require('./source.cjs'),{verify}=require('./verify.cjs'),{phaseEngine}=require('../native-lmd-first-mode-pulse-solver/phase.cjs');
const {f,iv,point,add,sub,neg,mul,scale,exp,sinh,absmax,width,contains,overlap,widen,trig,pack,productPartition,angular0}=R;
const directory=process.argv[2];assert(directory,'output directory required');
const read=p=>JSON.parse(fs.readFileSync(p)),write=(name,value)=>fs.writeFileSync(path.join(directory,name),JSON.stringify(value,null,2)+'\n');
const request=read(path.join(directory,'request.json'));
assert.deepEqual(request,read(path.join(__dirname,'request.json')),'commit a new protocol before changing the campaign');
assert.equal(request.schema,'lmd-nonuniform-first-mode-pulse-v1');
assert.equal(request.q_radius,'2');assert.equal(request.radial_endpoint,'16');assert.equal(request.gamma_mixture_max_index,32);
assert.equal(request.maximum_probability_halfwidth,'1/200');
assert.deepEqual(request.partitions,[{q_bins:64,radial_bins:16},{q_bins:128,radial_bins:32}]);
assert.deepEqual(request.angular_series_rules,[['1',16],['4',24],['16',40],['32',64]]);
const pins=read(path.join(__dirname,'SOURCE_PINS.json'));
for(const [p,hash] of Object.entries(pins.sha256))assert.equal(crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,p))).digest('hex'),hash,'source pin:'+p);
const inherited=path.join(__dirname,'../native-lmd-certified-pulse-solver/results');
const original=read(path.join(inherited,'solver.results.json'));
const source=read(path.join(__dirname,'../native-lmd-pulse-characteristic/results/source_intervals.json'));
const verified=verify(),checks=verified.checks,negativeControls=verified.negatives;
const check=(n,ok)=>{assert(ok,n);checks.push(n);},reject=(n,ok)=>{assert(ok,n);negativeControls.push(n);};
const unpack=x=>iv(x.lower,x.upper),engine=phaseEngine(inherited,check),target=f(request.maximum_probability_halfwidth);
const cache=new Map(),sources=new Map(),summaries=[];
const boundary=x=>x===null?'-infinity':String(x);

function partition(h,e,bins){
 const eabs=f(e).abs().toString(),sourcekey=[h,bins.q_bins,bins.radial_bins].join(':'),key=sourcekey+':'+eabs;
 if(cache.has(key))return cache.get(key);
 const p=productPartition(h,eabs,bins,unpack(source.integrals[h]));
 check('source_probability:'+key,contains(mul(p.q.totalP,p.radial.totalP),1));
 check('angle_averaged_source_mean:'+key,contains(mul(p.q.totalA,p.radial.totalA),1));
 check('variance_decomposition:'+key,overlap(p.residual,Object.values(p.decomposition).reduce(add,point(0))));
 check('positive_conditional_source:'+key,p.rows.every(r=>r.mass1.hi.mul(eabs).le(r.mass0.lo)));
 if(!sources.has(sourcekey))sources.set(sourcekey,{h,...bins,cells:p.rows.length,
  q_probability:pack(p.q.totalP),q_mean:pack(p.q.totalA),radial_probability:pack(p.radial.totalP),radial_uniform_mean:pack(p.radial.totalA),
  radial_first_harmonic_integral:pack(p.radial.first1Total),mixture_remainder:p.radial.mixtureRemainder.toString(),
  q_cells:p.q.rows.map(r=>({left:boundary(r.left),right:boundary(r.right),p:pack(r.p),first:pack(r.first),mass:pack(r.mass)})),
  radial_cells:p.radial.rows.map(r=>({left:boundary(r.left),right:boundary(r.right),p:pack(r.p),first0:pack(r.first),first1:pack(r.first1),mass0:pack(r.mass),mass1:pack(r.mass1)}))});
 summaries.push({h,e_abs:eabs,...bins,cells:p.rows.length,full_second:pack(p.fullSecond),ideal_first_mode_second:pack(p.idealSecond),
  captured_second:pack(p.second),source_variance:pack(p.residual),ideal_first_mode_variance:pack(p.floor),
  variance_decomposition:Object.fromEntries(Object.entries(p.decomposition).map(([k,v])=>[k,pack(v)]))});
 cache.set(key,p);console.log('Prepared weighted source',key,'cells',p.rows.length);return p;
}

function solve(c,bins){
 const h=f(c.h),t=f(c.t),k=f(c.k),e=f(c.e),p=partition(c.h,c.e,bins);
 const a=scale(sinh(point(h.mul(k).div(2))),t.mul(2).div(h)),damping=exp(point(h.mul(k.pow(2)).div(4).neg()));
 let re=point(0),im=point(0),totalP=point(0),fallbackProbability=point(0),seriesError=point(0),fallbackCells=0;
 const sampled=[],indices=new Set([0,Math.floor(p.rows.length/2),p.rows.length-1]);
 for(let index=0;index<p.rows.length;index++){
  const row=p.rows[index],angle=mul(a,row.mass0),phase=engine.phase(angle),arg=scale(mul(a,row.mass1),e),angular=angular0(arg);
  const weight=mul(row.p,angular.value);
  re=add(re,mul(weight,phase.cos));im=add(im,mul(weight,phase.sin));totalP=add(totalP,row.p);
  if(angular.fallback){fallbackCells++;fallbackProbability=add(fallbackProbability,row.p);}
  else seriesError=add(seriesError,scale(row.p,angular.tail));
  if(indices.has(index)){
   const native=trig(angle),label=[c.id,bins.q_bins,bins.radial_bins,index].join(':');
   check('selected_phase_cos:'+label,overlap(phase.cos,native.cos));check('selected_phase_sin:'+label,overlap(phase.sin,native.sin));
   sampled.push({index,angle:pack(angle),interpolated:{cos:pack(phase.cos),sin:pack(phase.sin)},independent_taylor:{cos:pack(native.cos),sin:pack(native.sin)},
    angular_argument:pack(arg),angular_factor:pack(angular.value),angular_terms:angular.terms,angular_fallback:angular.fallback});
  }
 }
 check('assembled_probability:'+c.id+':'+bins.q_bins,contains(totalP,1));
 const compression=absmax(a).pow(2).mul(p.residual.hi).div(2),real=mul(damping,widen(re,compression)),imaginary=mul(damping,widen(im,compression));
 const clamp=x=>iv(x.lo.n<0n?0:x.lo,f(1).le(x.hi)?1:x.hi);
 const pc=clamp(scale(add(point(1),real),'1/2')),ps=clamp(scale(add(point(1),imaginary),'1/2'));
 const hc=width(pc).div(2),hs=width(ps).div(2),passed=hc.le(target)&&hs.le(target);
 const transfer=scale(mul(damping,mul(a,a)),'1/4'),ideal=mul(transfer,p.floor);
 console.log('Evaluated',c.id,bins.q_bins,bins.radial_bins,'tolerance',passed?'met':'not met');
 return{...bins,cells:p.rows.length,p_cos_plus:pack(pc),p_sin_plus:pack(ps),phi_cos:pack(real),phi_sin:pack(imaginary),
  phase_argument:pack(a),damping:pack(damping),probability_halfwidth:{cos:hc.toString(),sin:hs.toString()},
  tolerance_status:passed?'certified_within_requested_halfwidth':'certified_enclosure_tolerance_not_met',
  error_budget:{source_variance:pack(p.residual),source_component_error:compression.toString(),
   source_probability_components:Object.fromEntries(Object.entries(p.decomposition).map(([k,v])=>[k,pack(mul(transfer,v))])),
   ideal_first_mode_probability_budget:pack(ideal),ideal_budget_exceeds_target:!ideal.lo.le(target),
   phase_component_error:t.zero()?'0':engine.error.toString(),angular_series_weighted_error:pack(seriesError),
   angular_fallback_cells:fallbackCells,angular_fallback_probability:pack(fallbackProbability),
   arithmetic_and_special_function_intervals:'propagated through final endpoints',
   interpretation:'variance-based certificate; ideal budget is not an actual-error lower bound'},selected_independent_phase_checks:sampled};
}

const cases=request.cases.map(c=>({...c,refinements:request.partitions.map(b=>solve(c,b))}));
const byId=id=>cases.find(c=>c.id===id);
for(const c of cases){
 const coarse=partition(c.h,c.e,request.partitions[0]),fine=partition(c.h,c.e,request.partitions[1]);
 check('nested_variance:'+c.id,fine.residual.hi.le(coarse.residual.lo));
}
for(let j=0;j<request.partitions.length;j++){
 const p=byId('small_h').refinements[j],s=byId('amplitude_sign').refinements[j],n=byId('time_reverse').refinements[j],z=byId('zero_pulse').refinements[j];
 check('amplitude_sign_cos:'+j,overlap(unpack(p.phi_cos),unpack(s.phi_cos)));
 check('amplitude_sign_sin:'+j,overlap(unpack(p.phi_sin),unpack(s.phi_sin)));
 check('time_reversal_cos:'+j,overlap(unpack(p.phi_cos),unpack(n.phi_cos)));
 check('time_reversal_sin:'+j,overlap(unpack(p.phi_sin),neg(unpack(n.phi_sin))));
 check('zero_pulse_cos:'+j,overlap(unpack(z.phi_cos),unpack(z.damping)));
 check('zero_pulse_sin:'+j,contains(unpack(z.phi_sin),0));
 const old=original.cases.find(c=>c.h==='1/8'&&c.t==='1/2'&&c.k==='2').refinements.at(-1).actual;
 const u=byId('uniform_reference').refinements[j];
 check('uniform_predecessor_cos_overlap:'+j,overlap(unpack(old.p_cos_plus),unpack(u.p_cos_plus)));
 check('uniform_predecessor_sin_overlap:'+j,overlap(unpack(old.p_sin_plus),unpack(u.p_sin_plus)));
}
reject('pulse_direction_detected',!overlap(unpack(byId('small_h').refinements.at(-1).phi_sin),unpack(byId('time_reverse').refinements.at(-1).phi_sin)));
reject('weighted_second_rejects_uniform_source',!overlap(partition('3/8','1/2',request.partitions.at(-1)).fullSecond,unpack(source.integrals['3/8'])));
const uniform=byId('uniform_reference').refinements.at(-1),comparisons=[];
for(const id of ['small_h','strong_modulation']){
 const c=byId(id),r=c.refinements.at(-1),dc=sub(unpack(r.p_cos_plus),unpack(uniform.p_cos_plus)),ds=sub(unpack(r.p_sin_plus),unpack(uniform.p_sin_plus));
 comparisons.push({id,e:c.e,cos_difference:pack(dc),sin_difference:pack(ds),cos_disjoint:!contains(dc,0),sin_disjoint:!contains(ds,0)});
}
const primaryIds=['small_h','middle_h','large_h','difficult','strong_modulation'];
const pass=c=>c.refinements.at(-1).tolerance_status==='certified_within_requested_halfwidth';
const report={status:'PASS',checks:checks.length,negative_controls:negativeControls.length,check_names:checks,negative_control_names:negativeControls,
 target_halfwidth:target.toString(),cases:cases.length,cases_within_tolerance:cases.filter(pass).length,
 primary_nonuniform_cases:primaryIds.length,primary_nonuniform_cases_within_tolerance:cases.filter(c=>primaryIds.includes(c.id)&&pass(c)).length,
 unresolved_cases:cases.filter(c=>!pass(c)).map(c=>c.id),uniform_comparisons:comparisons,
 fresh_nrqc_requests:0,reused_nrqc_requests:9,hardware_called:false,weighted_full_mass_law_assumed:false,
 relative_phase_integrated_with_modulation:true,gravitational_dynamics_derived:false,advantage_claimed:false,
 summary:cases.map(c=>{const r=c.refinements.at(-1);return{id:c.id,h:c.h,t:c.t,k:c.k,e:c.e,
  cos_plus:[r.p_cos_plus.display_lower,r.p_cos_plus.display_upper],sin_plus:[r.p_sin_plus.display_lower,r.p_sin_plus.display_upper],
  tolerance:r.tolerance_status,ideal_budget_exceeds_target:r.error_budget.ideal_budget_exceeds_target};})};
write('solver.results.json',{request,sources:[...sources.values()],partitions:summaries,cases});write('verification.json',report);
console.log(JSON.stringify({status:report.status,checks:report.checks,negative_controls:report.negative_controls,passes:report.cases_within_tolerance,
 primary_passes:report.primary_nonuniform_cases_within_tolerance,unresolved:report.unresolved_cases,comparisons,summary:report.summary},null,2));
