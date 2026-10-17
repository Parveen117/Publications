'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs'),W=require('../native-prime-innovation-window/window_adapter.cjs'),L=require('../native-infinite-layer-cut/layer_adapter.cjs'),D=require('./depth_adapter.cjs');
const {F}=N,{Interval:I}=O,args=process.argv.slice(2);assert.ok(args.length===1&&['--write','--check'].includes(args[0]));
const checks=[],negative_controls=[],fixtures=[],proofs=[],presentations=[],ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});},no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const positive=v=>!F.of(v).le(0),rejects=f=>{try{f();return false;}catch{return true;}},overlap=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.hi)&&b.lo.le(a.hi);};
const eq=N.equal,mul=N.mul,adj=N.dagger,det=W.determinant,sum=a=>a.reduce((s,v)=>s.add(v),new F(0)),quad=(M,u)=>mul(adj(u),mul(N.inverse(M),u))[0][0].rad;
const atoms=['1/5','4/5','9/25','16/25'].map(F.of),weights=['1/8','2/8','3/8','2/8'].map(F.of),c=Array.from({length:20},(_,j)=>sum(weights.map((v,k)=>v.mul(atoms[k].pow(j))))),toy=[];
for(let n=1;n<=3;n++){
 const old=L.exactForms(c,n),next=L.exactForms(c,n+1),E=N.matrix(Array.from({length:n+1},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0))),star=mul(N.inverse(old.H),mul(adj(E),next.H)),h=N.matrix(c.slice(n,2*n).map(v=>[v])),alpha=mul(N.inverse(old.H),h),v=N.matrix([...alpha.map(row=>[row[0].neg()]),[1]]),ell=N.matrix(Array.from({length:n},(_,j)=>[j===n-1?1:0])),norm=det(next.H).div(det(old.H)),previous=n===1?new F(1):det(old.H).div(det(L.exactForms(c,n-1).H)),beta=norm.div(previous),leak=N.sub(mul(next.Y,E),mul(E,old.Y)),R=N.sub(mul(star,mul(N.power(next.Y,2),E)),N.power(old.Y,2));
 ck('toy:isometry_and_linear_compression:'+n,eq(mul(star,E),N.identity(n))&&eq(mul(star,mul(next.Y,E)),old.Y));
 ck('toy:one_boundary_direction:'+n,eq(leak,mul(v,adj(ell)))&&N.rank(leak)===1);
 ck('toy:orthogonal_boundary_norm:'+n,N.isZero(mul(adj(E),mul(next.H,v)))&&mul(adj(v),mul(next.H,v))[0][0].rad.eq(norm));
 ck('toy:rank_one_square_defect:'+n,eq(R,N.scale(mul(N.inverse(old.H),mul(ell,adj(ell))),norm))&&N.rank(R)===1);
 ck('toy:coupling_eigenvalue_and_native_interval:'+n,sum(R.map((row,j)=>row[j].rad)).eq(beta)&&positive(beta)&&positive(new F('1/4').sub(beta))&&D.boundary(c,n).beta.contains(beta));
 const u=N.matrix(c.slice(0,n).map(v=>[v])),un=N.matrix(c.slice(0,n+1).map(v=>[v])),innovations=[];
 for(const raw of ['-1/3','0','1/7','2']){
  const lambda=F.of(raw),p=N.add(old.H,N.scale(old.A,lambda)),pn=N.add(next.H,N.scale(next.A,lambda)),difference=quad(pn,un).sub(quad(p,u)),formula=lambda.pow(2*n).mul(det(next.H).pow(2)).div(det(p).mul(det(pn))),native=D.innovation(c,n,lambda);
  ck('toy:determinant_innovation:'+n+':'+raw,difference.eq(formula)&&native.formula.contains(formula)&&native.difference.contains(formula));
  ck('toy:strict_increase_or_anchor_equality:'+n+':'+raw,lambda.zero()?difference.zero():positive(difference));
  innovations.push({lambda,difference,formula});
 }
 let power=N.identity(n);const model=[];for(let j=0;j<=2*n;j++){model.push(mul(adj(N.matrix(Array.from({length:n},(_,i)=>[i===0?1:0]))),mul(old.H,mul(power,N.matrix(Array.from({length:n},(_,i)=>[i===0?1:0])))))[0][0].rad);power=mul(power,old.Y);}
 ck('toy:stable_readings_through_2N_minus1:'+n,model.slice(0,2*n).every((v,j)=>v.eq(c[j])));
 ck('toy:first_new_even_reading_is_boundary_norm:'+n,c[2*n].sub(model[2*n]).eq(norm));
 no('toy:all_readings_stabilize_at_finite_depth:'+n,!model[2*n].eq(c[2*n]));
 toy.push({n,norm,previous,beta,E,adjoint:star,boundary_vector:v,leak,square_defect:R,model_moments:model,innovations});
}
fixtures.push({kind:'exact_nested_carriers',atoms,weights,moments:c,transitions:toy});

const symmetric=Array.from({length:6},(_,j)=>new F('1/5').pow(j).add(new F('4/5').pow(j)).div(2)),small=L.exactForms(symmetric,1),large=L.exactForms(symmetric,2),e=N.matrix([[1],[0]]),star=mul(N.inverse(small.H),mul(adj(e),large.H)),v=L.block(e,N.zeros(2,1),N.zeros(2,1),e),vs=L.block(star,N.zeros(1,2),N.zeros(1,2),star),cut1=L.cutWave(small.Y,[[new F('1/2')]],small.H),cut2=L.cutWave(large.Y,N.scale(N.identity(2),'2/5'),large.H),compressed=mul(vs,mul(cut2.cut,v)),defect=N.sub(compressed,N.power(compressed,2)),compressedWave=mul(vs,mul(cut2.wave,v));
ck('cut:exact_compression_defect',eq(compressed,[['1/2','2/5'],['2/5','1/2']])&&eq(defect,N.scale(N.identity(2),'9/100')));
ck('cut:full_projection_identity',eq(N.power(cut1.cut,2),cut1.cut)&&eq(N.power(cut2.cut,2),cut2.cut));
no('cut:natural_inclusion_intertwines',!eq(mul(cut2.cut,v),mul(v,cut1.cut))&&!eq(compressed,cut1.cut));
no('wave:compression_remains_unitary',eq(mul(adj(compressedWave),compressedWave),N.scale(N.identity(2),'16/25'))&&!eq(mul(adj(compressedWave),compressedWave),N.identity(2)));
fixtures.push({kind:'exact_cut_compression_witness',moments:symmetric,compressed,defect,compressed_wave:compressedWave});

for(const m of [1,2,5,18]){
 const pair=D.adjacent(c,m);
 ck('toy:deletion_and_squared_error:'+m,pair.coarse.every((v,j)=>v.contains(pair.lower.weights[j].lo)&&overlap(v,pair.lower.weights[j]))&&overlap(pair.meanSquare,pair.expected));
 ck('toy:all_factorial_observers:'+m,Array.from({length:m+1},(_,r)=>D.reading(pair.lower.weights,r).contains(c[r])).every(Boolean));
 ck('toy:conditional_mean_identity:'+m,Array.from({length:m+2},(_,j)=>{const mean=new F(m+1-j,m+1).mul(new F(j,m)).add(new F(j,m+1).mul(new F(j-1,m)));return mean.eq(new F(j,m+1));}).every(Boolean));
}
no('guard:missing_next_source_pair',rejects(()=>D.boundary(c.slice(0,18),9)));
no('guard:record_beyond_source',rejects(()=>D.record(c,20)));
no('guard:record_without_positive_words',rejects(()=>D.record([1,2],1)));
no('guard:nonunit_source',rejects(()=>D.record([2,1],1)));
no('guard:zero_adjacent_record_length',rejects(()=>D.adjacent(c,0)));
no('guard:nonpositive_response_domain',rejects(()=>D.innovation(c,2,-1)));
no('guard:saturated_finite_source_extends_strictly',rejects(()=>D.boundary(c,4)));

const T=new P.Presentation({tokens:['P','Q'],rules:[{id:'P_projection',lhs:['P','P'],rhs:[[['P'],1]],source:'orthogonal carrier projection'},{id:'Q_projection',lhs:['Q','Q'],rhs:[[['Q'],1]],source:'deeper exact cut'}]}),audit=T.audit();ck('words:compression_presentation',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name:'cut_compression',spec:T.spec(),audit});
const p=T.word(['P']),q=T.word(['Q']),one=T.one(),corner=p.times(q).times(p);
function prove(name,left,right){const z=T.reduce(left.minus(right),{witness:true});ck('words:'+name,z.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(z.certificate));ck('words:replay:'+name,T.replay(certificate));proofs.push({name,presentation_hash:T.hash(),lhs:left.toJSON(),rhs:right.toJSON(),certificate});}
prove('compression_defect_is_discarded_channel',corner.minus(corner.times(corner)),p.times(q).times(one.minus(p)).times(q).times(p));
prove('retained_plus_discarded_norm',p.times(q).times(p).times(q).times(p).plus(p.times(q).times(one.minus(p)).times(q).times(p)),corner);

const parent=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-layer-anchor-atlas/CERTIFICATE.json'),'utf8')),atlas=parent.fixtures.find(x=>x.kind==='actual_source_anchor_atlas');assert.ok(atlas);
const sources=[];
for(const source of atlas.sources){
 const jet=source.jet,raw=jet.moments,records=Array.from({length:20},(_,m)=>D.record(raw,m)),transitions=[],couplings=[],lambda=F.of(atlas.held_out.jet.t).div(jet.a).sub(1);
 ck('actual:retained_twenty_layers:'+jet.s,jet.degree===19&&raw.length===20&&[29,31,33].includes(jet.s));
 ck('actual:all_record_rows_normalized:'+jet.s,records.every(row=>row.total.contains(1)));
 ck('actual:all_factorial_readings_retained:'+jet.s,records.every(row=>Array.from({length:row.m+1},(_,r)=>overlap(D.reading(row.weights,r),raw[r])).every(Boolean)));
 ck('actual:every_available_deletion:'+jet.s,records.slice(1).every((row,m)=>D.coarsen(row.weights).every((v,j)=>overlap(v,records[m].weights[j]))));
 for(let m=1;m<=18;m++){
  const pair=D.adjacent(raw,m);
  ck('actual:adjacent_record_error:'+jet.s+':'+m,pair.meanSquare.positive()&&overlap(pair.meanSquare,pair.expected));
  couplings.push({m,meanSquare:pair.meanSquare,expected:pair.expected});
 }
 for(let n=1;n<=9;n++){
  const boundary=D.boundary(raw,n),innovation=D.innovation(raw,n,lambda);
  ck('actual:positive_rank_one_boundary:'+jet.s+':'+n,boundary.norm.positive()&&boundary.beta.positive()&&positive(new F('1/4').sub(boundary.beta.hi))&&overlap(boundary.norm,boundary.schurNorm));
  ck('actual:monic_orthogonality:'+jet.s+':'+n,boundary.orthogonality.every(v=>v.contains(0)));
  ck('actual:determinant_schur_innovation:'+jet.s+':'+n,innovation.formula.positive()&&innovation.difference.positive()&&innovation.schur.positive()&&overlap(innovation.difference,innovation.formula)&&overlap(innovation.residual,innovation.residualFormula)&&overlap(innovation.schur,innovation.schurRatio));
  transitions.push({boundary,innovation});
 }
 const last=transitions[8],retained=atlas.held_out.responses.find(x=>x.source_s===jet.s),newResponse=last.innovation.after.mul(jet.q0),sourceError=I.of(atlas.held_out.jet.q).sub(newResponse),scale=I.of(raw[1]).sub(raw[2]);
 ck('actual:depth10_heldout_response_reproduced:'+jet.s,overlap(newResponse,retained.response)&&overlap(sourceError,retained.error)&&sourceError.positive());
 ck('actual:record_noise_scale_positive:'+jet.s,scale.positive()&&scale.hi.le('1/4'));
 sources.push({s:jet.s,a:jet.a,q0:jet.q0,moments:raw,records,adjacent_record_couplings:couplings,transitions,held_out:{s:35,lambda,response:newResponse,error:sourceError},conditional_L2_scale:scale});
}
no('actual:finite_records_supply_all_order_premise',parent.scope.all_order_source_positivity===false&&parent.scope.rh_proved===false);
no('actual:strict_growth_allows_exact_cut_intertwining',sources.every(s=>s.transitions.every(t=>t.boundary.norm.positive())));
fixtures.push({kind:'actual_depth_compatibility',source_certificate:'../native-layer-anchor-atlas/CERTIFICATE.json',sources,scope:'Nine strict transitions and records through length19 at each printed anchor; no new derivative or unbounded source premise supplied'});
const payload={protocol:'LDC1_LAYER_DEPTH_COMPATIBILITY_V1',status:'PASS_LDC1_LAYER_DEPTH_COMPATIBILITY',source_parent:'e05800d68fbbcce229785cb303e6bf6b59cef086',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,scope:{physical_interpretation:'open; no predetermined field identification',stable_word_probability_readings:true,rank_one_depth_boundary:true,exact_positive_response_innovation:true,natural_cut_and_wave_intertwining_obstructed_under_strict_growth:true,conditional_record_L2_limit:true,all_order_actual_source_positivity:false,source_anchors:[29,31,33],highest_derivative:19,new_source_derivatives:0,largest_carrier:10,largest_exact_source_record:19,arithmetic:'frozen native F/Cut and inherited outward intervals; Python3.12 I/O only'},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length,actual_source_anchors:3,certified_depth_transitions:27,positive_response_innovations:27,adjacent_record_couplings:54,largest_record_length:19}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(args[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
