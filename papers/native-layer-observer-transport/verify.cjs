'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs'),L=require('../native-infinite-layer-cut/layer_adapter.cjs'),C=require('../native-layer-continuation-window/continuation_adapter.cjs'),T=require('./transport_adapter.cjs');
const {F}=N,{Interval:I}=O,args=process.argv.slice(2);assert.ok(args.length===1&&['--write','--check'].includes(args[0]));
const checks=[],negative_controls=[],fixtures=[],presentations=[],proofs=[],ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});},no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const eq=N.equal,mul=N.mul,sub=N.sub,adj=N.dagger,scale=N.scale,id=N.identity,positive=v=>!F.of(v).le(0),sum=a=>a.reduce((s,v)=>s.add(v),new F(0));
const overlap=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.hi)&&b.lo.le(a.hi);},rejects=f=>{try{f();return false;}catch{return true;}};
const atoms=['1/10','1/5','2/5','3/5','4/5'].map(F.of),weights=['1/10','2/10','3/10','1/10','3/10'].map(F.of),c=Array.from({length:20},(_,j)=>sum(weights.map((v,k)=>v.mul(atoms[k].pow(j))))),toy=[];
for(let n=1;n<=4;n++){
 const a=T.exactTransport(c,n),{old,next,V,star,K,p,w,norm,boundary,compressed,compressedWave}=a,unit=id(2*n),nextUnit=id(2*(n+1)),projected=mul(V,mul(old.cut,star)),leak=sub(mul(next.cut,V),mul(V,old.cut)),leakStar=T.adjoint(leak,old.G,next.G),cutDefect=sub(compressed,N.power(compressed,2)),oldComplement=sub(unit,old.cut),ZKZ=mul(old.Z,mul(K,old.Z)),waveStar=T.adjoint(compressedWave,old.G,old.G),waveLeak=sub(mul(next.wave,V),mul(V,old.wave)),waveLeakStar=T.adjoint(waveLeak,old.G,next.G),window=C.exactWindow(c,n,c[2*n]);
 ck('toy:weighted_exact_projection:'+n,eq(N.power(old.cut,2),old.cut)&&eq(mul(adj(old.cut),old.G),mul(old.G,old.cut)));
 ck('toy:isometry_and_channel_intertwining:'+n,eq(mul(star,V),unit)&&eq(mul(next.Z,V),mul(V,old.Z)));
 ck('toy:weighted_wave_unitarity:'+n,eq(mul(T.adjoint(next.wave,next.G,next.G),next.wave),nextUnit));
 ck('toy:positive_cut_sector_retained:'+n,eq(mul(next.cut,mul(V,old.cut)),mul(V,old.cut)));
 ck('toy:exact_rank_one_projection_update:'+n,eq(sub(next.cut,projected),boundary)&&eq(N.power(boundary,2),boundary)&&N.rank(boundary)===1);
 ck('toy:orthogonal_new_boundary:'+n,N.isZero(mul(projected,boundary))&&eq(mul(adj(w),mul(next.G,w)),[[norm]]));
 ck('toy:compressed_new_boundary:'+n,eq(mul(star,mul(boundary,V)),K)&&eq(compressed,N.add(old.cut,K)));
 ck('toy:leakage_norm_and_complement_support:'+n,eq(mul(leakStar,leak),K)&&N.isZero(mul(old.cut,K))&&eq(mul(oldComplement,K),K));
 ck('toy:canonical_position_is_only_leakage_eigenvalue:'+n,eq(N.power(K,2),scale(K,p))&&N.rank(K)===1&&p.eq(norm.div(window.width))&&positive(p)&&positive(new F(1).sub(p)));
 ck('toy:projection_compression_defect:'+n,eq(cutDefect,scale(K,new F(1).sub(p))));
 ck('toy:compressed_wave_norm_loss:'+n,eq(sub(unit,mul(waveStar,compressedWave)),scale(ZKZ,new F(4).mul(new F(1).sub(p)))));
 ck('toy:wave_transport_discrepancy:'+n,eq(mul(waveLeakStar,waveLeak),scale(ZKZ,4)));
 const z=N.matrix(Array.from({length:2*n},(_,j)=>[j%3-1])),carried=mul(V,z),observer=N.add(old.Z,scale(old.cut,'2/3')),carriedObserver=mul(V,mul(observer,star)),before=mul(adj(z),mul(old.G,mul(observer,z))),after=mul(adj(carried),mul(next.G,mul(carriedObserver,carried))),boundaryReading=mul(adj(carried),mul(next.G,mul(boundary,carried)))[0][0].rad;
 ck('toy:co_carried_observer_reading:'+n,eq(before,after));
 ck('toy:boundary_readout_formula:'+n,boundaryReading.eq(norm.mul(z[n-1][0].rad.sub(z[2*n-1][0].rad).pow(2))));
 no('toy:all_observers_unchanged_without_carrying:'+n,!eq(compressed,old.cut));
 no('toy:exact_full_cut_intertwining:'+n,!N.isZero(leak));
 no('toy:compressed_wave_is_unitary:'+n,!eq(mul(waveStar,compressedWave),unit));
 const interval=T.sourceTransport(c,n);
 ck('toy:native_interval_contains_exact_position:'+n,interval.position.contains(p)&&interval.windowPosition.contains(p));
 if(n<=3){
  const deep=T.exactTransport(c,n,n+2),tail=T.exactTransport(c,n+1,n+2);
  ck('toy:composition_of_observer_transport:'+n,eq(deep.V,mul(tail.V,V))&&eq(deep.star,mul(star,tail.star)));
  ck('toy:one_layer_saturation:'+n,eq(deep.compressed,compressed)&&eq(mul(deep.next.cut,deep.V),mul(tail.V,mul(next.cut,V))));
 }
 toy.push({n,norm,canonical_position:p,window_width:window.width,V,adjoint:star,cut:old.cut,next_cut:next.cut,boundary_projection:boundary,compressed_cut:compressed,leakage:K,compressed_wave:compressedWave,boundary_readout:boundaryReading});
}
fixtures.push({kind:'exact_weighted_transport',atoms,weights,moments:c,transitions:toy});

// A rational physical square-root witness makes the change of metric explicit.
const single=[1,'9/25'],physical=L.cutWave([['9/25']],[['12/25']]),weighted=T.weightedCut(single,1),R=N.matrix([['3/5',0],[0,'4/5']]);
ck('coordinates:finite_physical_cut_equivalence',eq(mul(adj(R),R),weighted.G)&&eq(mul(physical.cut,R),mul(R,weighted.cut))&&eq(mul(physical.wave,R),mul(R,weighted.wave)));
no('coordinates:weighted_pairing_is_identity',!eq(weighted.G,id(2)));

const nested=new P.Presentation({tokens:['P','Q'],rules:[{id:'P_squared',lhs:['P','P'],rhs:[[['P'],1]],source:'old embedded orthogonal cut'},{id:'Q_squared',lhs:['Q','Q'],rhs:[[['Q'],1]],source:'new orthogonal cut'},{id:'nested_PQ',lhs:['P','Q'],rhs:[[['P'],1]],source:'nested cut range'},{id:'nested_QP',lhs:['Q','P'],rhs:[[['P'],1]],source:'nested cut range'}]}),wave=new P.Presentation({tokens:['C','Z'],rules:[{id:'channel_involution',lhs:['Z','Z'],rhs:[[[],1]],source:'fixed channel reflection'}]});
function prove(presentation,name,left,right){const audit=presentation.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:presentation.spec(),audit});const result=presentation.reduce(left.minus(right),{witness:true});ck('words:identity:'+name,result.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(result.certificate));ck('words:replay:'+name,presentation.replay(certificate));proofs.push({name,presentation_hash:presentation.hash(),lhs:left.toJSON(),rhs:right.toJSON(),certificate});}
const pword=nested.word(['P']),qword=nested.word(['Q']),delta=qword.minus(pword);prove(nested,'nested_boundary_is_projection',delta.times(delta),delta);
const cw=wave.word(['C']),zw=wave.word(['Z']),one=wave.one(),reflection=cw.plus(cw).minus(one),gap=cw.minus(cw.times(cw)),four=gap.plus(gap).plus(gap).plus(gap);prove(wave,'wave_compression_loss',one.minus(zw.times(reflection).times(reflection).times(zw)),zw.times(four).times(zw));

const parent=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-layer-depth-compatibility/CERTIFICATE.json'),'utf8')),retained=parent.fixtures.find(x=>x.kind==='actual_depth_compatibility');assert.ok(retained);
const sources=[];
for(const s of retained.sources){
 const transitions=[];
 for(let n=1;n<=9;n++){
  const t=T.sourceTransport(s.moments,n),old=s.transitions[n-1];
  ck('actual:strict_channel_pairings:'+s.s+':'+n,t.position.positive()&&new I(1).sub(t.position).positive());
  ck('actual:canonical_position_identity:'+s.s+':'+n,overlap(t.position,t.windowPosition)&&overlap(t.inverseWidth,new I(1).div(t.windowWidth)));
  ck('actual:retained_boundary_reproduced:'+s.s+':'+n,overlap(t.boundaryNorm,old.boundary.norm));
  ck('actual:nonzero_cut_and_wave_defects:'+s.s+':'+n,t.cutDefect.positive()&&t.waveLoss.positive()&&t.waveDiscrepancy.positive()&&t.waveLoss.hi.le(1)&&t.waveDiscrepancy.hi.le(4));
  transitions.push(t);
 }
 sources.push({s:s.s,a:s.a,transitions});
}
fixtures.push({kind:'actual_observer_transport',source_certificate:'../native-layer-depth-compatibility/CERTIFICATE.json',sources,scope:'Universal root-free identities specialized to 27 certified finite source transitions; no approximate matrix square roots or new derivatives'});
no('guard:uncertified_source_pairing',rejects(()=>T.weightedCut([1,2],1)));
no('guard:insufficient_next_source_layers',rejects(()=>T.sourceTransport(c.slice(0,18),9)));
no('guard:nonunit_source',rejects(()=>T.weightedCut([2,1],1)));
no('guard:finite_saturated_measure_extends_strictly',rejects(()=>T.sourceTransport(c,5)));
no('guard:invalid_carrier_dimensions',rejects(()=>T.inclusion(2,1)));
no('guard:no_depth_increase',rejects(()=>T.exactTransport(c,2,2)));
no('actual:finite_certification_supplies_infinite_premise',parent.scope.all_order_actual_source_positivity===false);
const payload={protocol:'LOT1_LAYER_OBSERVER_TRANSPORT_V1',status:'PASS_LOT1_LAYER_OBSERVER_TRANSPORT',source_parent:'96b9903491ff3dbd71835e4c944fd3306ef333bf',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,scope:{physical_interpretation:'open; source structure determines subsequent interpretation',weighted_channel_isometric_transport:true,nested_cut_ranges:true,rank_one_projection_update:true,boundary_leakage_equals_even_canonical_position:true,co_carried_observer_readings_preserved:true,fixed_cut_reading_preserved_for_every_state:false,natural_H_doubled_inclusion_intertwining:false,conditional_monotone_strong_cut_limit:true,all_order_actual_source_positivity:false,rh_proved:false,source_anchors:[29,31,33],highest_derivative:19,new_source_derivatives:0,largest_carrier:10,arithmetic:'frozen native F/Cut and inherited outward intervals; Python3.12 I/O only'},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length,actual_source_anchors:3,certified_depth_transitions:27}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(args[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
