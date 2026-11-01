'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const S=require('../native-infinite-layer-cut/source_layers.cjs');
const Q=require('../native-prime-innovation-window/full_source.cjs');
const LAT=require('../native-layer-anchor-transport/anchor_adapter.cjs');
const A=require('./atlas_adapter.cjs');
const {F}=N,{Interval:I}=O,argv=process.argv.slice(2);assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});},no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const pos=x=>!F.of(x).le(0),rejects=f=>{try{f();return false;}catch{return true;}},encloses=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.lo)&&b.hi.le(a.hi);},overlap=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.hi)&&b.lo.le(a.hi);};
const json=x=>JSON.stringify(x,(_,v)=>typeof v==='bigint'?v.toString():v),sum=xs=>xs.reduce((a,b)=>a.add(b),new F(0));

const tailRecords=[];
for(const [K,j,z] of [[0,0,'1/3'],[3,0,'1/2'],[3,2,'1/4'],[7,7,'1/2'],[19,5,'1/5'],[19,19,'1/2']]){
 const v=A.exactTail(K,j,z),end=K+7,partial=sum(Array.from({length:7},(_,i)=>{const n=K+1+i;return new F(L.choose(n,j)).mul(F.of(z).pow(n-j));})),rest=A.exactTail(end,j,z);
 ck('tail:finite_split_identity:'+K+':'+j,v.eq(partial.add(rest)));
 ck('tail:positive_remainder:'+K+':'+j,pos(v)&&pos(rest));
 tailRecords.push({K,j,z,tail:v,split_at:end,partial,remainder:rest});
}
for(const [K,j,z] of [[3,0,'1/4'],[3,2,'1/4'],[19,5,'4/25']])ck('tail:improves_old_geometric_bound:'+K+':'+j,A.exactTail(K,j,z).le(LAT.tail(K,j,z)));
ck('tail:zero_displacement',Array.from({length:20},(_,j)=>A.exactTail(19,j,0)).every(x=>x.zero()));
ck('tail:top_layer_remains_bounded_inside_disk',pos(A.exactTail(19,19,'1/2'))&&rejects(()=>LAT.tail(19,19,'1/2')));
no('guard:tail_at_boundary',rejects(()=>A.exactTail(19,5,1)));
no('guard:missing_tail_layer',rejects(()=>A.exactTail(19,20,'1/2')));
no('guard:negative_tail_radius',rejects(()=>A.exactTail(19,1,'-1/2')));
fixtures.push({kind:'exact_positive_tail',records:tailRecords});

for(const K of [0,1,4,19])for(const [rho,sigma] of [['5/4','7/6'],['7/8','8/7']]){
 const r=F.of(rho),s=F.of(sigma),one=A.matrix(K,1),R=A.matrix(K,r),T=A.matrix(K,s),total=A.matrix(K,r.mul(s)),inverse=A.matrix(K,new F(1).div(r));
 ck('polynomial:composition:'+K+':'+rho,N.equal(N.mul(T,R),N.matrix(total)));
 ck('polynomial:identity_and_inverse:'+K+':'+rho,N.equal(one,N.identity(K+1))&&N.equal(N.mul(inverse,R),N.identity(K+1)));
}
const toy=Array.from({length:4},(_,j)=>new F('1/2').pow(j)),rho=new F('5/4'),delta=rho.sub(1),matrix=A.matrix(3,rho),moved=matrix.map(row=>sum(row.map((v,j)=>v.mul(toy[j])))),actual0=new F(1).div(new F(1).add(delta.div(2))),actual1=rho.div(new F(2).add(delta));
const residual=delta.div(2).pow(4).div(new F(1).add(delta.div(2)));
ck('toy:exact_value_remainder',actual0.sub(moved[0]).eq(residual)&&pos(residual));
no('restart:shifted_polynomial_equals_actual_value',!moved[0].eq(actual0));
no('restart:normalization_removes_source_error',!moved[1].div(moved[0]).eq(actual1));
const shortened=moved.map((v,j)=>j===3?new F(0):v),back=A.matrix(3,'4/5').map(row=>sum(row.map((v,j)=>v.mul(shortened[j]))));
no('restart:discarded_polynomial_layer_is_free',!back[0].eq(toy[0]));
no('guard:nonpositive_polynomial_ratio',rejects(()=>A.matrix(3,0)));
fixtures.push({kind:'exact_polynomial_and_source_distinction',toy,rho,moved,actual_value:actual0,actual_first_moment:actual1,residual,shortened_inverse:back});

const T=new P.Presentation({tokens:['r','s','x'],rules:[{id:'s_r',lhs:['s','r'],rhs:[[['r','s'],1]],source:'commuting real anchor ratios'},{id:'x_r',lhs:['x','r'],rhs:[[['r','x'],1]],source:'commuting displacement and ratio'},{id:'x_s',lhs:['x','s'],rhs:[[['s','x'],1]],source:'commuting displacement and ratio'}]});
const audit=T.audit();ck('words:affine_transport_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name:'source_anchor_affine_composition',spec:T.spec(),audit});
const one=T.one(),r=T.word(['r']),s=T.word(['s']),x=T.word(['x']);
function prove(name,left,right){const z=T.reduce(left.minus(right),{witness:true});ck('words:'+name,z.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(z.certificate));ck('words:replay:'+name,T.replay(certificate));proofs.push({name,presentation_hash:T.hash(),lhs:left.toJSON(),rhs:right.toJSON(),certificate});}
prove('composed_anchor_variable',r.minus(one).plus(r.times(s.minus(one).plus(s.times(x)))),r.times(s).minus(one).plus(r.times(s).times(x)));
prove('geometric_remainder_K3',one.minus(x).times(one.plus(x).plus(x.times(x)).plus(x.times(x).times(x))),one.minus(x.times(x).times(x).times(x)));

const parent=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-layer-anchor-transport/CERTIFICATE.json'),'utf8')),retained=parent.fixtures.find(x=>x.kind==='actual_anchor_transport');assert.ok(retained);
const jets=[retained.jet,S.sourceLayers(31,19),S.sourceLayers(33,19)],charts=[],sources=[];
const piw=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-prime-innovation-window/CERTIFICATE.json'),'utf8')).fixtures.find(x=>x.kind==='complete_source_prime_update');assert.ok(piw);
ck('source:reference_majorant_exactly_retained',A.sourceContract(jets[0]).majorant.eq('221/228'));
ck('source:new_majorants',A.sourceContract(jets[1]).majorant.eq('225/244')&&A.sourceContract(jets[2]).majorant.eq('229/260'));
for(const jet of jets){
 const base=A.sourceContract(jet),forms=L.forms(jet.moments,10),triangle=[];
 for(let r=0;r<=19;r++)for(let k=0;k+r<=19;k++)triangle.push({r,k,value:L.difference(jet.moments,r,k)});
 ck('actual:complete_twenty_layer_jet:'+jet.s,jet.degree===19&&jet.moments.length===20&&base.a.eq(new F(2*jet.s-1,2).pow(2)));
 ck('actual:tenth_order_forms:'+jet.s,forms.factors.A.positive&&forms.factors.B.positive&&forms.factors.H.positive);
 ck('actual:210_positive_binomial_readings:'+jet.s,triangle.length===210&&triangle.every(v=>v.value.positive()));
 ck('actual:retained_complete_tails:'+jet.s,jet.components.hurwitz.length===19&&jet.components.hurwitz.every(v=>pos(v.remainder))&&jet.components.euler.tails.length===20&&jet.components.euler.tails.every(pos));
 ck('actual:source_interval_precision:'+jet.s,jet.moments.every(v=>{v=I.of(v);return v.hi.sub(v.lo).le(O.decimal('1e-62'));}));
 const local=[A.corridor(jet,3,'3/20'),A.corridor(jet,4,'1/12')];
 ck('actual:continuous_order3_chart:'+jet.s,local[0].certified);
 ck('actual:continuous_order4_chart:'+jet.s,local[1].certified);
 charts.push(...local);sources.push({s:jet.s,jet,forms,binomial_triangle:triangle});
}
ck('actual:fresh31_independent_value_overlap',overlap(jets[1].q0,piw.held_out_jet.q));
ck('actual:fresh31_independent_derivative_overlap',overlap(I.of(jets[1].moments[1]).mul(jets[1].q0).div(F.of(jets[1].a).neg()),piw.held_out_jet.derivative));
const atlas3=A.glue(charts.filter(v=>v.m===3),3),atlas4=A.glue(charts.filter(v=>v.m===4),4);
ck('atlas:exact_order3_range',atlas3.range[0].eq('55233/80')&&atlas3.range[1].eq('19435/16'));
ck('atlas:exact_order4_range',atlas4.range[0].eq('11913/16')&&atlas4.range[1].eq('54925/48'));
ck('atlas:strict_overlaps',atlas3.overlaps.concat(atlas4.overlaps).every(v=>pos(v[1].sub(v[0]))));
const covers=(atlas,s1,s2)=>atlas.range[0].le(new F(2*s1-1,2).pow(2))&&new F(2*s2-1,2).pow(2).le(atlas.range[1]);
ck('atlas:every_real_s27_through35',covers(atlas3,27,35));
ck('atlas:every_real_s28_through34',covers(atlas4,28,34));
no('atlas:bridge_missing_order4_chart',rejects(()=>A.glue([charts[1],charts[5]],4)));
no('atlas:promote_order3_chart_to_order4',rejects(()=>A.glue([charts[0]],4)));
no('atlas:use_failed_certificate',rejects(()=>A.glue([{...charts[0],certified:false}],3)));
no('actual:positivity_requires_prime_anchor',jets[2].s===33&&33%3===0&&sources[2].forms.factors.A.positive&&sources[2].forms.factors.B.positive);

const direct=A.forecast(jets[0],['4225/3249'],19),chained=A.forecast(jets[0],['3721/3249','4225/3721'],19),fresh=A.forecast(jets[1],['4225/3721'],19),roundTrip=A.forecast(jets[0],['3721/3249','3249/3721'],19);
ck('ledger:direct_and_chained_same_reading',json(direct.raw)===json(chained.raw)&&json(direct.budget)===json(chained.budget));
ck('ledger:retains_original_source',chained.source_s===29&&chained.source_degree===19&&chained.source_a.eq('3249/4')&&chained.target_a.eq('4225/4'));
ck('ledger:round_trip_has_zero_remainder',roundTrip.rho.eq(1)&&roundTrip.budget.values.every(v=>v.zero())&&roundTrip.raw.every((v,j)=>encloses(v,jets[0].moments[j])));
ck('ledger:all_twenty_actual_target_layers_enclosed',direct.normalized.every((v,j)=>encloses(v,jets[2].moments[j]))&&fresh.normalized.every((v,j)=>encloses(v,jets[2].moments[j])));
ck('ledger:fresh_source_restart_reduces_every_budget',fresh.budget.values.every((v,j)=>pos(direct.budget.values[j].sub(v))));
ck('ledger:reported_fifth_layer_budgets',direct.budget.values[5].le(O.decimal('0.356'))&&fresh.budget.values[5].le(O.decimal('8.251e-7')));
no('ledger:highest_retained_layer_is_precise',pos(direct.budget.values[19].sub(1))&&pos(fresh.budget.values[19].sub(1)));
no('ledger:actual_polynomial_error_can_be_reset',!overlap(direct.polynomial[0].mul(jets[0].q0),jets[2].q0));
no('guard:cumulative_disk_escape_by_small_steps',rejects(()=>A.forecast(jets[0],Array(4).fill('5/4'),1)));
no('guard:lost_source_anchor_metadata',rejects(()=>A.sourceContract({...jets[0],a:'3721/4'})));
no('guard:unavailable_twenty_first_layer',rejects(()=>A.forecast(jets[0],['3721/3249'],20)));
no('guard:source_outside_retained_domain',rejects(()=>A.sourceContract({...jets[0],s:3,w:'5/2',a:'25/4'})));

const held=Q.sourceJet(35),forecast35=A.forecast(jets[2],['4761/4225'],5),heldRaw1=I.of(held.derivative).mul(F.of(held.t).neg()).div(jets[2].q0),responses=[];
ck('heldout:independent_s35_value_inside_forecast',held.s===35&&encloses(forecast35.q,held.q));
ck('heldout:independent_s35_derivative_inside_forecast',encloses(forecast35.raw[1],heldRaw1));
for(const jet of jets){
 const b=A.sourceContract(jet),lambda=F.of(held.t).div(b.a).sub(1),response=L.response(jet.moments,10,lambda).value.mul(jet.q0),error=I.of(held.q).sub(response),wh=F.of(held.w),bound=wh.add(b.w).pow(2).div(wh).mul(wh.sub(b.w).div(wh.add(b.w).sub(2)).pow(20));
 ck('heldout:source_scale_contract:'+jet.s,b.q0.hi.le('1/2')&&held.q.positive());
 ck('heldout:unconditional_error_bound:'+jet.s,error.positive()&&error.absUpper().le(bound));
 responses.push({source_s:jet.s,response,error,analytic_bound:bound});
}
ck('heldout:successive_fresh_models_improve',responses[1].error.hi.le(responses[0].error.lo)&&responses[2].error.hi.le(responses[1].error.lo));
ck('heldout:reported_s33_error',responses[2].error.hi.le(O.decimal('1.971e-35')));
no('atlas:rebuilt_models_are_one_transported_response',!overlap(responses[0].response,responses[1].response)&&!overlap(responses[1].response,responses[2].response));
no('atlas:finite_charts_establish_all_depth',parent.scope.all_order_source_positivity===false&&parent.scope.rh_proved===false);
fixtures.push({kind:'actual_source_anchor_atlas',sources,charts,atlases:[atlas3,atlas4],transport:{direct,chained,fresh_restart:fresh,round_trip:roundTrip},held_out:{jet:held,forecast:forecast35,responses},scope:'Three finite source jets and six whole-interval certificates; no all-order or physical propagation conclusion'});
const payload={protocol:'LAA1_LAYER_ANCHOR_ATLAS_V1',status:'PASS_LAA1_LAYER_ANCHOR_ATLAS',source_parent:'08698bf8ad428ddd67f308ffc7299376e1beab0f',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,scope:{exact_differentiated_tail:true,source_preserving_polynomial_transport:true,connected_actual_source_atlas:true,original_source_anchors:[29,31,33],derivative_degree:19,confluent_directions_at_each_anchor:10,continuous_order3_a:['55233/80','19435/16'],continuous_order4_a:['11913/16','54925/48'],continuous_order3_includes_s:[27,35],continuous_order4_includes_s:[28,34],fresh_source_anchors:[31,33],all_order_source_positivity:false,rh_proved:false,physical_mass_Maxwell_or_quantum_gravity_closed:false,global_eta_improved:false,arithmetic:'frozen native rational/rewrite engines and inherited outward intervals; Python 3.12 I/O only'},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length,source_anchors:3,fresh_derivative_layers:40,retained_derivative_layers_per_anchor:20,positive_binomial_readings:630,continuous_source_charts:6,connected_atlases:2}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
