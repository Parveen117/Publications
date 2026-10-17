'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const {F}=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const {Interval:I}=require('../native-prime-observation-closure/observation_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const W=require('../native-layer-wave-propagation/propagation_adapter.cjs');
const A=require('./analytic_adapter.cjs');
const args=process.argv.slice(2);assert.ok(args.length===1&&['--write','--check'].includes(args[0]));
const checks=[],negative_controls=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}};
const overlap=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.hi)&&b.lo.le(a.hi);};
const subset=(a,b)=>{a=I.of(a);b=I.of(b);return b.lo.le(a.lo)&&a.hi.le(b.hi);};
const sum=xs=>xs.reduce((s,x)=>s.add(x),new F(0));

// A saturating signed sequence checks every binomial sign and majorant power.
const bound=A.majorant(3),signed=Array.from({length:42},(_,n)=>new F(n%2?-3:3).div(A.R.pow(n)));
ck('majorant:exact_growth_constant',bound.A.eq('1984/225'));
for(let j=0;j<=20;j++){
 const d=L.exactDifference(signed,j+1,j);
 ck('majorant:signed_saturation:'+j,d.abs().eq(A.fullBound(j,bound)));
 const K=19,unknown=sum(Array.from({length:j+1},(_,ell)=>j+1+ell>K?new F(L.choose(j,ell)).mul(signed[j+1+ell]).mul(ell%2?-1:1):new F(0)));
 ck('majorant:missing_saturation:'+j,unknown.abs().eq(A.missing(K,j,bound)));
}
ck('majorant:first_missing_order',A.missing(19,9,bound).zero()&&!A.missing(19,10,bound).zero());
for(const t of [0,'1/2',1,2,4]){
 const rest=A.tail(20,t,bound),first=bound.C.mul(rest.x.pow(20)).div(W.factorial(40));
 ck('tail:geometric_identity:'+t,rest.radius.mul(new F(1).sub(rest.ratio)).eq(first));
 const nextRatio=rest.x.div(43*44);ck('tail:decreasing_ratios:'+t,nextRatio.le(rest.ratio));
}
fixtures.push({kind:'universal_majorant_fixtures',majorant:bound,first_missing_j:10});

// Positive rational atoms provide an exact source with a known global budget.
const atoms=['1/5','1/2','4/5'].map(F.of),c=Array.from({length:82},(_,n)=>sum(atoms.map(y=>y.pow(n))).div(3)),toy=[];
for(const K of [3,7,19])for(const t of ['1/2',1,2]){
 const raw=c.slice(0,K+1),b=A.majorant(1),read=A.enclosure(raw,t,b,24),p=A.known(raw,t),grouped=sum(raw.slice(1).map((v,n)=>v.mul(A.kernel(n+1,t))));
 ck('kernel:regrouping:'+K+':'+t,p.value.contains(grouped));
 ck('kernel:evenness:'+K+':'+t,A.known(raw,F.of(t).neg()).value.lo.eq(p.value.lo)&&A.known(raw,F.of(t).neg()).value.hi.eq(p.value.hi));
 // Independent finite cosine expansion at the three atoms, with a certified
 // positive-measure Taylor remainder. This is not a numerical transcendental call.
 const direct=sum(atoms.map(y=>sum(Array.from({length:40},(_,j)=>A.coefficient(j,t).mul(y.pow(j+1)).mul(new F(1).sub(y).pow(j)))))).div(3);
 const radius=c[1].mul(F.of(t).pow(80)).div(W.factorial(80)),oracle=new I(direct.sub(radius),direct.add(radius));
 ck('toy:positive_atom_enclosure:'+K+':'+t,subset(oracle,read.enclosure));
 toy.push({K,t,reading:read,cosineOracle:oracle});
}
fixtures.push({kind:'positive_atom_return',atoms,cases:toy});

// Entire return and first finite positivity do not imply positive realization.
const bad=[new F(1),new F(1,2),new F(-1,4)],badAtOne=A.known(bad,1).value;
ck('counterexample:explicit_entire_polynomial',badAtOne.contains(1)&&A.kernel(1,1).eq(1)&&A.kernel(2,1).eq(-2));
ck('counterexample:first_pair_positive',!bad[1].le(0)&&!bad[0].sub(bad[1]).le(0));
no('inference:entire_implies_unitary',!badAtOne.lo.le(bad[1]));
no('inference:first_pair_implies_all_orders',bad[2].sub(bad[1].pow(2)).eq('-1/2'));
fixtures.push({kind:'entire_nonunitary_counterexample',moments:bad,all_later_moments:'zero',return_polynomial:['1/2',0,'1/2'],at_one:badAtOne,hankel_determinant:'-1/2'});

const parent=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-layer-depth-compatibility/CERTIFICATE.json'),'utf8'));
const retained=parent.fixtures.find(x=>x.kind==='actual_depth_compatibility');assert.ok(retained);
const lwp=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-layer-wave-propagation/CERTIFICATE.json'),'utf8'));
const wave=lwp.fixtures.find(x=>x.kind==='actual_wave_propagation');assert.ok(wave);
const sources=[];
for(const source of retained.sources){
 const b=A.sourceMajorant(source),raw=source.moments,K=raw.length-1,readings=[];
 ck('actual:retained_order:'+source.s,K===19);
 ck('actual:source_budget:'+source.s,b.A.eq('1984/225')&&b.q0.positive());
 for(const t of ['1/2',1,2]){
  const coarse=A.coarse(raw,t,b),refined=A.enclosure(raw,t,b,24),known=A.known(raw,t),grouped=raw.slice(1).reduce((v,x,n)=>v.add(I.of(x).mul(A.kernel(n+1,t))),new I(0));
  ck('actual:regrouping:'+source.s+':'+t,overlap(known.value,grouped));
  ck('actual:radius_improves:'+source.s+':'+t,refined.radius.le(coarse.radius));
  ck('actual:coarse_refined_consistency:'+source.s+':'+t,overlap(coarse.enclosure,refined.enclosure));
  ck('actual:source_uncertainty_retained:'+source.s+':'+t,refined.enclosure.lo.le(refined.partial.lo.sub(refined.radius))&&refined.partial.hi.add(refined.radius).le(refined.enclosure.hi));
  if(t!== '1/2'){
   const prior=wave.sources.find(x=>String(x.s)===String(source.s)).exponential.find(x=>F.of(x.t).eq(t));assert.ok(prior);
   ck('actual:LWP_prefix_matches:'+source.s+':'+t,coarse.partial.lo.eq(I.of(prior.partial).lo)&&coarse.partial.hi.eq(I.of(prior.partial).hi));
   ck('actual:conditional_interval_consistent:'+source.s+':'+t,overlap(refined.enclosure,prior.enclosure));
  }
  readings.push({t,coarse,refined,normalized:refined.enclosure.div(raw[1])});
 }
 const zero=A.enclosure(raw,0,b,24);
 ck('actual:initial_value:'+source.s,zero.radius.zero()&&zero.enclosure.lo.eq(I.of(raw[1]).lo)&&zero.enclosure.hi.eq(I.of(raw[1]).hi));
 sources.push({s:source.s,a:source.a,majorant:b,highest_derivative:K,polynomial_degree:2*(K-1),readings});
}
fixtures.push({kind:'actual_entire_source_return',sources});
no('guard:nonpositive_majorant',rejects(()=>A.majorant(0)));
no('guard:zero_coefficient_radius',rejects(()=>A.majorant(1,0)));
no('guard:source_outside_domain',rejects(()=>A.sourceMajorant({...retained.sources[0],s:4})));
no('guard:inconsistent_anchor',rejects(()=>A.sourceMajorant({...retained.sources[0],a:1})));
no('guard:uncertified_q0',rejects(()=>A.sourceMajorant({...retained.sources[0],q0:[-1,1]})));
no('guard:invalid_tail_ratio',rejects(()=>A.tail(1,2,bound)));
no('guard:insufficient_tail_cutoff',rejects(()=>A.enclosure(c.slice(0,20),1,A.majorant(1),18)));
no('guard:coefficient_budget_mismatch',rejects(()=>A.enclosure([1,20],1,A.majorant(1),4)));
no('guard:missing_normalization',rejects(()=>A.enclosure([2,'1/2'],1,A.majorant(3),4)));
no('guard:empty_jet',rejects(()=>A.enclosure([],1,bound,4)));
no('inference:LWP_tiny_radius_unconditional',lwp.scope.all_order_actual_source_positivity===false);
const payload={protocol:'LAR1_LAYER_ANALYTIC_RETURN_V1',status:'PASS_LAR1_LAYER_ANALYTIC_RETURN',source_parent:'3a843da9377f1a03e54869e5f3ea22f82ac42499',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,fixtures,scope:{actual_entire_return:'unconditional using LAA source Cauchy bound; written analytic proof',all_known_coefficient_polynomial:true,explicit_unknown_tail:true,all_order_positive_realization:'conditional on the inherited all-order source criterion',rh_proved:false,physical_interpretation:'open; tau is dimensionless',new_source_derivatives:0,highest_derivative:19,arithmetic:'frozen native F and outward intervals; Python3.12 I/O only',execution:'GitHub Actions; no local scientific replay'},summary:{exact_checks:checks.length,symbolic_certificates:0,negative_controls:negative_controls.length,actual_source_anchors:sources.length,actual_source_enclosures:sources.reduce((s,x)=>s+x.readings.length,0)}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(args[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
