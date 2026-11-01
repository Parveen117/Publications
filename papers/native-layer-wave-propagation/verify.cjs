'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs'),O=require('../native-prime-observation-closure/observation_adapter.cjs'),T=require('../native-layer-observer-transport/transport_adapter.cjs'),L=require('../native-infinite-layer-cut/layer_adapter.cjs'),A=require('./propagation_adapter.cjs');
const {F}=N,{Interval:I}=O,args=process.argv.slice(2);assert.ok(args.length===1&&['--write','--check'].includes(args[0]));
const checks=[],negative_controls=[],fixtures=[],presentations=[],proofs=[],ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});},no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const eq=N.equal,mul=N.mul,adj=N.dagger,positive=x=>!F.of(x).le(0),sum=a=>a.reduce((s,v)=>s.add(v),new F(0)),rejects=f=>{try{f();return false;}catch{return true;}},overlap=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.hi)&&b.lo.le(a.hi);};
const atoms=['1/10','1/5','2/5','3/5','4/5','9/10'].map(F.of),c=Array.from({length:24},(_,j)=>sum(atoms.map(v=>v.pow(j))).div(6)),toy=[];
const embed=(n,m)=>{const E=T.inclusion(n,m),zero=N.zeros(m,n);return L.block(E,zero,zero,E);},quadratic=(left,metric,right)=>mul(adj(left),mul(metric,right))[0][0].rad;
for(let n=1;n<=3;n++)for(let m=n;m<=4;m++){
 const old=T.weightedCut(c,m),deep=T.weightedCut(c,m+1),G=A.generator(old),Gn=A.generator(deep),V=embed(m,m+1),u=[Array.from({length:n},(_,j)=>j+1),Array.from({length:n},(_,j)=>1-j)],z=[Array.from({length:n},(_,j)=>j===n-1?1:0),Array.from({length:n},(_,j)=>j%2?-1:2)],uc=A.pairColumn(u,m),zc=A.pairColumn(z,m),ud=mul(V,uc),zd=mul(V,zc),margin=m-n;
 ck('toy:skew_adjoint_generator:'+n+':'+m,eq(T.adjoint(G,old.G,old.G),N.scale(G,-1))&&eq(G,N.scale(N.sub(old.wave,T.adjoint(old.wave,old.G,old.G)),'1/2')));
 for(const kind of ['wave','inverse','generator']){
  const B=kind==='wave'?old.wave:kind==='inverse'?T.adjoint(old.wave,old.G,old.G):G,Bn=kind==='wave'?deep.wave:kind==='inverse'?T.adjoint(deep.wave,deep.G,deep.G):Gn;
  for(let k=0;k<=2*margin;k++){
   const v=mul(N.power(B,k),zc),vd=mul(N.power(Bn,k),zd),poly=A.action(z,kind,k),reading=quadratic(uc,old.G,v),deepReading=quadratic(ud,deep.G,vd),expected=A.functional(c,A.correlation(u,poly));
   ck('toy:exact_correlation_band:'+n+':'+m+':'+kind+':'+k,reading.eq(deepReading)&&expected.contains(reading));
   if(k<=margin)ck('toy:exact_state_band:'+n+':'+m+':'+kind+':'+k,eq(mul(V,v),vd)&&eq(v,A.pairColumn(poly,m)));
  }
 }
 toy.push({n,m,margin,generator:G,state_steps:margin,correlation_steps:2*margin});
}
fixtures.push({kind:'exact_state_and_observer_bands',atoms,moments:c,cases:toy});

const first=[];
for(let m=1;m<=4;m++){
 const k=2*m-1,old=T.weightedCut(c,m),next=T.weightedCut(c,m+1),z=A.pairColumn([[1],[0]],m),zn=A.pairColumn([[1],[0]],m+1),before=quadratic(z,old.G,mul(N.power(old.wave,k),z)),after=quadratic(zn,next.G,mul(N.power(next.wave,k),zn)),norm=T.exactTransport(c,m).norm,formula=new F(2).pow(4*m-3).mul(norm),poly=A.returnPolynomial(k);
 ck('toy:first_unresolved_return:'+m,after.sub(before).eq(formula)&&positive(formula)&&poly.length===2*m+1&&poly[2*m].eq(new F(2).pow(4*m-3)));
 no('toy:unrestricted_correlation_stability:'+m,!before.eq(after));
 const seed=A.pairColumn([[1],[0]],m),deeperSeed=mul(embed(m,m+1),seed),oldState=mul(N.power(old.wave,m),seed),newState=mul(N.power(next.wave,m),deeperSeed);
 no('toy:state_stability_beyond_margin:'+m,!eq(mul(embed(m,m+1),oldState),newState));
 const native=A.finiteReturns(c,m,k);ck('toy:native_interval_return:'+m,native.values[k].contains(before));
 first.push({m,k,before,after,boundaryNorm:norm,formula});
}
fixtures.push({kind:'first_boundary_return',transitions:first});
const gtoy=A.finiteReturns(c,5,8,'generator');
ck('toy:generator_return_coefficients',gtoy.values.every((v,k)=>overlap(v,A.functional(c,A.correlation([[1],[0]],A.action([[1],[0]],'generator',k))))));
const half=T.weightedCut([1,'1/2'],1),halfExp=A.exponentialReading([1,'1/2'],1,1),halfReturn=quadratic(A.pairColumn([[1],[0]],1),half.G,mul(half.wave,A.pairColumn([[1],[0]],1)));
no('exponential:unit_parameter_equals_discrete_wave',halfExp.enclosure.positive()&&halfReturn.zero());
const signed=A.functional(c,A.returnPolynomial(2));no('return:correlation_is_a_probability',signed.hi.le(0)&&!signed.hi.zero());

const presentation=new P.Presentation({tokens:['P','Z'],rules:[{id:'cut_projection',lhs:['P','P'],rhs:[[['P'],1]],source:'weighted source cut'},{id:'channel_reflection',lhs:['Z','Z'],rhs:[[[],1]],source:'channel involution'}]}),audit=presentation.audit();ck('words:confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name:'source_reflections',spec:presentation.spec(),audit});
const p=presentation.word(['P']),z=presentation.word(['Z']),one=presentation.one(),reflection=p.plus(p).minus(one),wave=reflection.times(z),back=z.times(reflection),generator=p.times(z).minus(z.times(p));
function prove(name,left,right){const r=presentation.reduce(left.minus(right),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,presentation.replay(certificate));proofs.push({name,presentation_hash:presentation.hash(),lhs:left.toJSON(),rhs:right.toJSON(),certificate});}
prove('wave_inverse',wave.times(back),one);prove('generator_is_skew_part',wave.minus(back),generator.plus(generator));prove('generator_flips_channel',z.times(generator).plus(generator.times(z)),one.minus(one));

const budgets=[];
for(const tau of [0,1,2,4]){
 const b=A.budgets(1,10,tau),r=10,firstTerm=new F(tau).pow(r).div(A.factorial(r)),ratio=new F(tau,r+1);
 ck('tail:positive_geometric_majorant:'+tau,b.relativeStateBound.eq(firstTerm.div(new F(1).sub(ratio)).mul(2))&&new F(0).le(b.relativeCorrelationBound));
 budgets.push(b);
}
ck('tail:twentieth_factorial_exact',A.factorial(20).eq('2432902008176640000'));

const parent=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-layer-depth-compatibility/CERTIFICATE.json'),'utf8')),retained=parent.fixtures.find(x=>x.kind==='actual_depth_compatibility');assert.ok(retained);
const sources=[];
for(const source of retained.sources){
 const raw=source.moments,models=Array.from({length:10},(_,j)=>A.finiteReturns(raw,j+1,j===9?18:2*j+1)),returns=[],innovations=[];
 for(let m=1;m<=10;m++){
  const model=models[m-1],rows=[];
  for(let k=0;k<=2*m-2;k++){const direct=A.functional(raw,A.returnPolynomial(k)),finite=model.values[k];ck('actual:resolved_wave_return:'+source.s+':'+m+':'+k,overlap(direct,finite));rows.push({k,direct,finite});}
  returns.push({m,rows});
  if(m<=9){const k=2*m-1,before=model.values[k],after=models[m].values[k],difference=after.sub(before),formula=I.of(source.transitions[m-1].boundary.norm).mul(new F(2).pow(4*m-3));ck('actual:first_boundary_return:'+source.s+':'+m,difference.positive()&&formula.positive()&&overlap(difference,formula));innovations.push({m,k,before,after,difference,formula});}
 }
 const generatorReturns=A.finiteReturns(raw,10,18,'generator');
 ck('actual:generator_coefficients_match_source:'+source.s,generatorReturns.values.every((v,k)=>overlap(v,A.functional(raw,A.correlation([[1],[0]],A.action([[1],[0]],'generator',k))))));
 ck('actual:odd_seed_coefficients_vanish:'+source.s,generatorReturns.values.every((v,k)=>k%2===0||v.contains(0)));
 const exponential=[1,2,4].map(t=>A.exponentialReading(raw,10,t));
 for(const e of exponential){ck('actual:positive_binomial_inputs:'+source.s+':'+e.t,e.terms.every(t=>t.moment.positive()));ck('actual:finite_exponential_radius:'+source.s+':'+e.t,e.normalizedRadius.eq(e.t.pow(20).div(A.factorial(20)))&&e.enclosure.lo.le(e.enclosure.hi));}
 ck('actual:unit_parameter_radius:'+source.s,exponential[0].normalizedRadius.eq(new F(1).div(A.factorial(20))));
 sources.push({s:source.s,a:source.a,seedNormSquared:raw[1],returns,innovations,generatorReturns,exponential});
}
fixtures.push({kind:'actual_wave_propagation',sources,budgets,scope:'All finite source conclusions use only retained c_0 through c_19; deeper extension comparisons remain conditional on existence of those certified extensions'});
no('guard:missing_source_coefficient',rejects(()=>A.functional([1,'1/2'],[0,0,1])));
no('guard:nonpositive_source_pairings',rejects(()=>A.finiteReturns([1,2],1,2)));
no('guard:unavailable_carrier',rejects(()=>A.exponentialReading(retained.sources[0].moments,11,1)));
no('guard:invalid_tail_ratio',rejects(()=>A.tailBound(2,3)));
no('guard:zero_tail_order',rejects(()=>A.tailBound(0,1)));
no('guard:invalid_carrier_budget',rejects(()=>A.budgets(3,2,0)));
no('guard:unknown_action',rejects(()=>A.action([[1],[0]],'untyped',1)));
no('guard:polynomial_outside_carrier',rejects(()=>A.pairColumn([[1,1],[0]],1)));
no('actual:finite_execution_proves_all_order_source_signs',parent.scope.all_order_actual_source_positivity===false);
const payload={protocol:'LWP1_LAYER_WAVE_PROPAGATION_V1',status:'PASS_LWP1_LAYER_WAVE_PROPAGATION',source_parent:'a5f3a6234bd75aceb6cf7f2e3ebe2c35c5af2fa5',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,scope:{exact_state_depth_band:true,exact_doubled_correlation_band:true,first_return_innovation:true,unitary_exponential_of_skew_part:true,continuous_interpolation_is_discrete_wave_logarithm:false,finite_source_exponential_enclosures:true,conditional_infinite_limit:true,all_order_actual_source_positivity:false,rh_proved:false,physical_interpretation:'open; parameter has no assigned physical time or propagation speed',highest_derivative:19,new_source_derivatives:0,largest_carrier:10,arithmetic:'frozen native F/Cut and inherited outward intervals; Python3.12 I/O only'},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length,actual_source_anchors:3,resolved_wave_return_readings:300,first_boundary_return_innovations:27,exponential_enclosures:9}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(args[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
