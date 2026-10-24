'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const S=require('../native-infinite-layer-cut/source_layers.cjs');
const C=require('../native-layer-continuation-window/continuation_adapter.cjs');
const A=require('./anchor_adapter.cjs');
const {F}=N,{Interval:I}=O,argv=process.argv.slice(2);assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});},no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const pos=x=>!F.of(x).le(0),overlap=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.hi)&&b.lo.le(a.hi);},rejects=f=>{try{f();return false;}catch{return true;}};
const sum=xs=>xs.reduce((a,b)=>a.add(b),new F(0)),moments=(ys,ws,count)=>Array.from({length:count},(_,j)=>sum(ws.map((w,k)=>w.mul(ys[k].pow(j)))));
const atoms=['1/5','4/5','9/25','16/25'].map(F.of),weights=['1/8','2/8','3/8','2/8'].map(F.of),c=moments(atoms,weights,22),records=[];
for(const rho of ['3/4','5/4','37/32']){
 const moved=A.exactMove(atoms,weights,rho),back=A.exactMove(moved.atoms,moved.weights,new F(1).div(rho)),twice=A.exactMove(moved.atoms,moved.weights,'7/6'),once=A.exactMove(atoms,weights,F.of(rho).mul('7/6'));
 ck('toy:probability_and_interior:'+rho,sum(moved.weights).eq(1)&&moved.atoms.every(y=>pos(y)&&pos(new F(1).sub(y))));
 ck('toy:inverse_anchor:'+rho,back.atoms.every((y,k)=>y.eq(atoms[k]))&&back.weights.every((w,k)=>w.eq(weights[k])));
 ck('toy:anchor_semigroup:'+rho,twice.atoms.every((y,k)=>y.eq(once.atoms[k]))&&twice.weights.every((w,k)=>w.eq(once.weights[k]))&&moved.normalization.mul(twice.normalization).eq(once.normalization));
 ck('toy:invariant_poles_and_residues:'+rho,moved.atoms.every((y,k)=>F.of(rho).mul(new F(1).sub(y)).div(y).eq(new F(1).sub(atoms[k]).div(atoms[k]))&&moved.normalization.mul(moved.weights[k]).mul(rho).div(y).eq(weights[k].div(atoms[k]))));
 const target=new F('7/3'),before=sum(weights.map((w,k)=>w.div(new F(1).add(target.sub(1).mul(atoms[k]))))),after=moved.normalization.mul(sum(moved.weights.map((w,k)=>w.div(new F(1).add(target.div(rho).sub(1).mul(moved.atoms[k]))))));
 ck('toy:invariant_rational_response:'+rho,before.eq(after));
 for(const m of [1,2,4]){
  const T=A.polynomialTransport(m,rho),inv=A.polynomialTransport(m,new F(1).div(rho)),next=A.polynomialTransport(m,'7/6'),total=A.polynomialTransport(m,F.of(rho).mul('7/6'));
  ck('toy:polynomial_inverse_and_composition:'+rho+':'+m,N.equal(N.mul(T,inv),N.identity(m))&&N.equal(N.mul(T,next),N.matrix(total)));
 }
 records.push(moved);
}
const dots=A.hierarchy(c),direct=c.slice(0,-1).map((_,j)=>sum(weights.map((w,k)=>w.mul(atoms[k].pow(j)).mul(c[1].sub(atoms[k]).add(new F(j).mul(new F(1).sub(atoms[k])))))));
ck('toy:all_moment_flow_readings',dots.every((z,j)=>z.contains(direct[j]))&&direct[0].zero());
for(const [r,k] of [[0,0],[1,0],[0,3],[2,4],[7,5]]){
 const D=L.exactDifference(c,r,k),next=L.exactDifference(c,r+1,k),dot=sum(weights.map((w,j)=>w.mul(atoms[j].pow(r)).mul(new F(1).sub(atoms[j]).pow(k)).mul(c[1].sub(atoms[j]).add(new F(r).mul(new F(1).sub(atoms[j]))).sub(atoms[j].mul(k)))));
 ck('toy:binomial_flow:'+r+':'+k,dot.eq(new F(r).add(c[1]).mul(D).sub(next.mul(r+k+1))));
 const rho=new F('5/4'),moved=A.exactMove(atoms,weights,rho),cm=moments(moved.atoms,moved.weights,14),formula=sum(weights.map((w,j)=>w.mul(atoms[j].pow(r)).mul(new F(1).sub(atoms[j]).pow(k)).div(new F(1).add(rho.sub(1).mul(atoms[j])).pow(r+k+1)))).mul(rho.pow(r)).div(moved.normalization);
 ck('toy:binomial_anchor_transform:'+r+':'+k,L.exactDifference(cm,r,k).eq(formula));
}
fixtures.push({kind:'exact_anchor_semigroup_and_hierarchy',atoms,weights,moments:c,transports:records,flow:direct});
const unweighted=A.exactMove(atoms,weights,'5/4');
no('anchor:unchanged_weights',!moments(unweighted.atoms,weights,2)[1].eq(moments(unweighted.atoms,unweighted.weights,2)[1]));
no('anchor:unchanged_spectral_coordinates',!moments(atoms,unweighted.weights,2)[1].eq(moments(unweighted.atoms,unweighted.weights,2)[1]));
const changed=c.slice(0,21);changed[20]=changed[20].add(1);
no('finite_hierarchy:closed_without_next_layer',A.hierarchy(c.slice(0,21))[19].hi.le(A.hierarchy(changed)[19].lo)===false&&!overlap(A.hierarchy(c.slice(0,21))[19],A.hierarchy(changed)[19]));
no('guard:nonpositive_anchor',rejects(()=>A.exactMove(atoms,weights,0)));
no('guard:unnormalized_weights',rejects(()=>A.exactMove(atoms,[1,1,1,1],2)));
no('guard:noncontractive_tail_ratio',rejects(()=>A.tail(19,19,'1/2')));
no('guard:unsupported_source_disk',rejects(()=>A.errorBounds('1/100',19,'1/10','1/2')));
no('guard:incomplete_surrogate',rejects(()=>A.corridor([1,'1/2'],1,19,1,'1/10')));
no('guard:nonpositive_surrogate',rejects(()=>A.corridor([1,2],1,1,1,'1/10')));
no('guard:nonpositive_forecast_anchor',rejects(()=>A.forecast(c,1,0,1)));
no('guard:missing_forecast_depth',rejects(()=>A.forecast(c,1,'5/4',22)));

const T=new P.Presentation({tokens:['y','r','s'],rules:[{id:'r_y',lhs:['r','y'],rhs:[[['y','r'],1]],source:'commuting real spectral coordinate and anchor ratio'},{id:'s_y',lhs:['s','y'],rhs:[[['y','s'],1]],source:'commuting scalar ratios'},{id:'s_r',lhs:['s','r'],rhs:[[['r','s'],1]],source:'commuting scalar ratios'}]});
const audit=T.audit();ck('words:scalar_anchor_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name:'anchor_mobius_identities',spec:T.spec(),audit});
const one=T.one(),y=T.word(['y']),r=T.word(['r']),s=T.word(['s']),den=one.plus(r.minus(one).times(y));
function prove(name,left,right){const z=T.reduce(left.minus(right),{witness:true});ck('words:'+name,z.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(z.certificate));ck('words:replay:'+name,T.replay(certificate));proofs.push({name,presentation_hash:T.hash(),lhs:left.toJSON(),rhs:right.toJSON(),certificate});}
prove('mobius_composed_denominator',one.plus(s.minus(one).times(y)).plus(r.minus(one).times(s).times(y)),one.plus(r.times(s).minus(one).times(y)));
prove('complementary_branch',den.minus(r.times(y)),one.minus(y));
prove('logistic_flow_cross_multiplied',r.times(y).times(den).minus(r.times(r).times(y).times(y)),r.times(y).times(one.minus(y)));

const parent=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-layer-continuation-window/CERTIFICATE.json'),'utf8')),old=parent.fixtures.find(x=>x.kind==='actual_two_layer_continuation');assert.ok(old);
const jet=S.sourceLayers(29,19),a=C.acceptance(jet.moments,9,jet.moments[18],jet.moments[19]),forms=L.forms(jet.moments,10);
ck('actual:complete_20_source_layers',jet.s===29&&jet.a.eq('3249/4')&&jet.moments.length===20);
ck('actual:old_18_layers_overlap',jet.moments.slice(0,18).every((z,j)=>overlap(z,old.jet.moments[j])));
ck('actual:predicted_c18_window_overlap',overlap(a.even.lower,old.next_unexecuted_even_window.lower)&&overlap(a.even.upper,old.next_unexecuted_even_window.upper));
ck('actual:next_pair_passes_sharp_windows',a.accepted);
const inside=(z,lo,hi)=>pos(z.lo.sub(O.decimal(lo)))&&pos(O.decimal(hi).sub(z.hi));
ck('actual:canonical_positions',inside(a.evenFraction,'0.137','0.138')&&inside(a.oddFraction,'0.136','0.137'));
ck('actual:new_pivot_floors',pos(a.oddLower.lo.sub(O.decimal('3.479e-15')))&&pos(a.oddUpper.lo.sub(O.decimal('2.196e-14'))));
ck('actual:source_interval_precision',jet.moments.every(z=>z.hi.sub(z.lo).le(O.decimal('1e-62'))));
ck('actual:both_tenth_order_forms_positive',forms.factors.A.positive&&forms.factors.B.positive&&forms.factors.H.positive);
ck('actual:all_new_source_tails_retained',jet.components.hurwitz.length===19&&jet.components.hurwitz.every(z=>pos(z.remainder))&&jet.components.euler.tails.length===20&&jet.components.euler.tails.every(pos));
const triangle=[];for(let r=0;r<=19;r++)for(let k=0;k+r<=19;k++)triangle.push({r,k,value:L.difference(jet.moments,r,k)});
ck('actual:210_positive_binomial_readings',triangle.length===210&&triangle.every(z=>z.value.positive()));
const corridors=[A.corridor(jet.moments,jet.q0,19,3,'3/20'),A.corridor(jet.moments,jet.q0,19,4,'1/12')];
ck('actual:continuous_order3_corridor',corridors[0].certified);
ck('actual:continuous_order4_corridor',corridors[1].certified);
const wider=A.corridor(jet.moments,jet.q0,19,3,'1/6');
no('actual:larger_corridor_automatically_certified',!wider.certified);
const rho=new F('3721/3249'),forecast=A.forecast(jet.moments,jet.q0,rho,5),primeParent=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-prime-innovation-window/CERTIFICATE.json'),'utf8')),held=primeParent.fixtures.find(x=>x.kind==='complete_source_prime_update').held_out_jet;
ck('actual:source31_inside_continuous_corridor',rho.sub(1).abs().le('3/20'));
ck('actual:source27_inside_continuous_corridor',new F('2809/3249').sub(1).abs().le('3/20'));
ck('actual:held_out_value_inside_forecast',forecast.q.lo.le(I.of(held.q).lo)&&I.of(held.q).hi.le(forecast.q.hi));
const heldRaw1=I.of(held.derivative).mul(new F('3721/4').neg()).div(jet.q0);
ck('actual:held_out_derivative_inside_forecast',forecast.raw[1].lo.le(heldRaw1.lo)&&heldRaw1.hi.le(forecast.raw[1].hi));
const fit=L.response(jet.moments,10,rho.sub(1)).value.mul(jet.q0),error=I.of(held.q).sub(fit),oldfit=I.of(old.held_out.layer_response),wh=new F('61/2'),analytic=wh.add(jet.w).pow(2).div(wh).mul(wh.sub(jet.w).div(wh.add(jet.w).sub(2)).pow(20));
ck('actual:held_out_response_improves',error.positive()&&fit.sub(oldfit).positive()&&error.hi.le(I.of(old.held_out.error).lo));
ck('actual:unconditional_analytic_error_bound',error.absUpper().le(analytic));
ck('actual:reported_response_bounds',error.hi.le(O.decimal('4.160e-35'))&&analytic.le(O.decimal('9.131e-28')));
no('actual:finite_surrogate_equals_whole_source',error.positive());
no('actual:finite_corridor_proves_all_depths',parent.scope.rh_proved===false&&parent.scope.all_order_actual_source_window_membership_proved===false);
fixtures.push({kind:'actual_anchor_transport',jet,acceptance:a,forms,binomial_triangle:triangle,corridors,wider_corridor_attempt:{m:3,d:'1/6',certified:false,meaning:'this sufficient bound is inconclusive, not a negative actual source form'},forecast31:forecast,held_out:{s:31,q:held.q,derivative:held.derivative,layer_response:fit,error,analytic_bound:analytic},scope:'20 actual layers at s29; continuum positivity only at orders3 and4 on printed intervals'});
const payload={protocol:'LAT1_LAYER_ANCHOR_TRANSPORT_V1',status:'PASS_LAT1_LAYER_ANCHOR_TRANSPORT',source_parent:'d1654412a0bd6f572d688eb7986bbb5ebdb74ed1',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,scope:{exact_source_anchor_hierarchy:true,model_mobius_semigroup:true,unconditional_continuous_actual_finite_positivity:true,anchor_a:'3249/4',source_derivatives:20,anchor_confluent_directions:10,continuous_order3_rho:['17/20','23/20'],continuous_order4_rho:['11/12','13/12'],all_order_source_positivity:false,rh_proved:false,physical_mass_Maxwell_or_quantum_gravity_closed:false,global_eta_improved:false,arithmetic:'frozen native rational/rewrite engines and inherited outward intervals; Python 3.12 I/O only'},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length,actual_derivative_layers:20,actual_confluent_directions:10,positive_binomial_readings:210,continuous_source_corridors:2}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
