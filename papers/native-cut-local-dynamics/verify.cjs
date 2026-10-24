'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,mm=(...a)=>a.reduce(N.mul),sc=N.scale,eq=N.equal,zero=N.isZero,cm=N.commutator;
const checks=[],negative_controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const tr=a=>a.reduce((s,row,i)=>s.add(row[i]),N.ZERO),real=z=>{assert.ok(z.turn.zero());return z.rad;},ex=(rho,a)=>real(tr(mm(rho,a)));
const id=N.identity(2),p=N.matrix([[1,0],[0,0]]),q=N.sub(id,p),R=N.matrix([[0,-1],[1,0]]),W=N.matrix([['3/5','-4/5'],['4/5','3/5']]),pd=mm(W,p,N.dagger(W)),qd=N.sub(id,pd),K=sc(R,'3/5'),dp=cm(K,pd),Z=N.zeros(2);
ck('moving:projection',eq(mm(pd,pd),pd)&&eq(N.dagger(pd),pd));
ck('moving:differentiated_projection',eq(N.add(mm(dp,pd),mm(pd,dp)),dp));
ck('moving:connection',eq(cm(dp,pd),K));
ck('moving:antihermitian_connection',eq(N.dagger(K),sc(K,-1)));
ck('moving:tangent_offdiagonal',zero(mm(pd,dp,pd))&&zero(mm(qd,dp,qd)));
for(const [label,h]of [['stationary_source',Z],['comoving_source',sc(K,N.IOTA)],['other_source',N.matrix([[2,1],[1,-1]])]]){
 const defect=N.add(dp,sc(cm(h,pd),N.IOTA)),ht=mm(N.dagger(W),N.sub(h,sc(K,N.IOTA)),W),b=mm(q,ht,p);
 ck(label+':defect_hermitian',eq(defect,N.dagger(defect)));
 ck(label+':defect_offdiagonal',zero(mm(pd,defect,pd))&&zero(mm(qd,defect,qd)));
 ck(label+':moving_h_hermitian',eq(ht,N.dagger(ht)));
 ck(label+':connection_block_identity',eq(b,sc(mm(N.dagger(W),qd,defect,pd,W),N.IOTA.neg())));
 for(const [stateName,rho]of [['basis0',p],['basis1',q],['superposition',N.matrix([['1/2','1/2'],['1/2','1/2']])]]){
  const rhoPrime=sc(cm(h,rho),N.IOTA.neg());
  ck(label+':weight_derivative:'+stateName,tr(N.add(mm(dp,rho),mm(pd,rhoPrime))).eq(tr(mm(defect,rho))));
 }
 if(label==='comoving_source')ck(label+':closed',zero(defect)&&zero(ht));
 else no(label+':arbitrary_moving_cut_closes',!zero(b));
 const gauge=N.matrix([[N.Cut.of(['3/5','4/5']),0],[0,N.Cut.of(['5/13','12/13'])]]),diag=N.matrix([[2,0],[0,-3]]),hg=N.add(mm(N.dagger(gauge),ht,gauge),diag),bg=mm(q,hg,p);
 ck(label+':gauge_unitary',eq(mm(N.dagger(gauge),gauge),id));
 ck(label+':gauge_block_covariance',eq(bg,mm(N.dagger(gauge),b,gauge)));
 ck(label+':gauge_squared_norm',eq(mm(N.dagger(bg),bg),mm(N.dagger(gauge),N.dagger(b),b,gauge)));
}
const projected=mm(p,pd,p),returned=mm(p,qd,p);
ck('arrows:retained_amplitude',eq(projected,sc(p,'9/25')));
ck('arrows:return_amplitude',eq(returned,sc(p,'16/25')));
ck('arrows:coherent_sum',eq(N.add(projected,returned),p));
ck('arrows:selected_joint_probability',ex(p,mm(N.dagger(projected),projected)).eq('81/625'));
ck('arrows:no_intervention_survival',ex(p,p).eq(1));
no('arrows:reset_preserves_endpoint_probability',!ex(p,mm(N.dagger(projected),projected)).eq(1));
no('arrows:zero_initial_hidden_means_no_return',!zero(returned));
ck('arrows:comoving_step0',zero(mm(qd,W,p)));
ck('arrows:comoving_step1',zero(mm(q,N.dagger(W),pd)));
ck('arrows:comoving_composition',eq(mm(p,N.dagger(W),pd,W,p),p));
ck('detector:stationary_source_weight',ex(p,pd).eq('9/25'));
ck('detector:comoving_source_weight',ex(pd,pd).eq(1));
no('detector:omit_connection',!ex(p,pd).eq(1));
ck('detector:stationary_weight_without_closure',ex(sc(id,'1/2'),dp).zero()&&!zero(dp));
no('detector:one_state_conservation_implies_all_state',!zero(dp));
results.two_step={full_endpoint_survival:'1',retained_path_amplitude:'9/25',return_path_amplitude:'16/25',selected_joint_probability:'81/625'};
results.detector={stationary_source_moving_detector:'9/25',comoving_source_detector:'1'};
const budget=f('1/10'),leak=budget.pow(2),survival=f(1).sub(leak);
ck('budget:survival',survival.eq('99/100'));ck('budget:amplitude',leak.div(2).eq('1/200'));
ck('budget:unconditional',leak.eq('1/100'));ck('budget:conditional',leak.mul(2).div(survival).eq('2/99'));
no('budget:conditioning_at_zero_survival',(()=>{try{f(1).div(0);return false;}catch{return true;}})());
results.budget={Lambda:'1/10',survival_lower:'99/100',pure_amplitude_error_upper:'1/200',unconditioned_unit_observable_error_upper:'1/100',conditioned_unit_observable_error_upper:'2/99'};
function rotor(radius){const n=2*radius+1,id=N.identity(n),u=N.zeros(n),l=N.zeros(n),s=N.zeros(n);
 for(let j=0;j<n;j++){l[j][j]=N.Cut.of(j-radius);s[n-1-j][j]=N.ONE;if(j+1<n)u[j+1][j]=N.ONE;}
 const ud=N.dagger(u),sin=sc(N.sub(u,ud),N.Cut.of([0,'-1/2'])),sin2=sc(N.sub(mm(u,u),mm(ud,ud)),N.Cut.of([0,'-1/2'])),hg=N.sub(N.add(mm(l,l),sc(id,8)),N.add(mm(u,u),mm(ud,ud))),v=N.sub(sc(id,'13/10'),sc(N.add(u,ud),'3/20'));
 return{n,radius,id,u,ud,l,s,sin,sin2,hg,v,plus:sc(N.add(id,s),'1/2'),minus:sc(N.sub(id,s),'1/2')};}
function cpow(z,n){let a=N.ONE,base=n<0?z.inv():z;for(let j=0;j<Math.abs(n);j++)a=a.mul(base);return a;}
const fixtures=[];
for(const radius of [2,3]){
 const r=rotor(radius);
 for(const [c,s]of [['1','0'],['-1','0'],['0','1'],['3/5','4/5']]){
  const frame=N.zeros(r.n),phase=N.Cut.of([c,f(s).neg()]);for(let j=0;j<r.n;j++)frame[j][j]=cpow(phase,j-radius);
  ck('rotor:unitary:'+radius+':'+c,eq(mm(N.dagger(frame),frame),r.id));
  for(const k of [0,2])for(const omega of ['0','3/5']){
   const label=[radius,c,s,k,omega].join(':'),h=N.sub(r.hg,sc(r.v,k)),ht=N.sub(mm(N.dagger(frame),h,frame),sc(r.l,omega)),b=mm(r.minus,ht,r.plus);
   const expected=mm(r.minus,N.sub(N.sub(sc(r.sin2,f(c).mul(s).mul(4)),sc(r.sin,f(k).mul('3/10').mul(s))),sc(r.l,omega)),r.plus);
   ck('rotor:odd_block:'+label,eq(b,expected));
   ck('rotor:hermitian:'+label,eq(ht,N.dagger(ht)));
   if(s==='0'&&omega==='0')ck('rotor:stationary_closed:'+label,zero(b));
   if(c==='0'&&k===2&&omega==='0')no('rotor:geometry_only_closure_survives_matter:'+radius,!zero(b));
   if(c==='1'&&omega!=='0')no('rotor:instantaneous_aligned_cut_suffices:'+radius,!zero(b));
   if(radius===3)fixtures.push({cos_alpha:c,sin_alpha:s,matter_wave_number:k,omega,offblock_zero:zero(b)});
  }
 }
 const ket=N.zeros(r.n,1);ket[radius-1][0]=N.ONE;ket[radius+1][0]=N.ONE;
 const rho=sc(mm(ket,N.dagger(ket)),'1/2'),b=sc(mm(r.minus,r.l,r.plus),'-3/5');
 ck('rotor:initial_cosine_retained:'+radius,ex(rho,r.plus).eq(1));
 ck('rotor:initial_leakage_coefficient:'+radius,ex(rho,mm(N.dagger(b),b)).eq('9/25'));
}
results.rotor={scope:'Symmetric finite Fourier fixtures of the common-core identity; not proof of unbounded evolution or spatial constraints.',fixtures,initial_cosine1_omega_three_fifths_leakage_coefficient:'9/25'};
function register(name,spec){const s=new M.Presentation(spec),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(spec.tokens.map(t=>{const n=typeof t==='string'?t:t.name;return[n,s.word([n])];}))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{p,d,h}}=register('differentiated_projection',{tokens:['p','d','h'],rules:[{id:'pp',lhs:['p','p'],rhs:[[['p'],1]],source:'Orthogonal projection idempotence'},{id:'dp',lhs:['d','p'],rhs:[[['d'],1],[['p','d'],-1]],source:'Derivative of projection identity'}]});
 const q=s.one().minus(p),K=d.times(p).minus(p.times(d)),D=d.plus(h.times(p).minus(p.times(h)).scale(N.IOTA));
 prove('projection:retained_tangent_zero',p.times(d).times(p),s.zero());
 prove('projection:hidden_tangent_zero',q.times(d).times(q),s.zero());
 prove('projection:connection_transports_cut',K.times(p).minus(p.times(K)),d);
 prove('projection:connection_retained_zero',p.times(K).times(p),s.zero());
 prove('projection:defect_retained_zero',p.times(D).times(p),s.zero());
 prove('projection:defect_hidden_zero',q.times(D).times(q),s.zero());
 prove('projection:defect_block_sign',q.times(h.minus(K.scale(N.IOTA))).times(p),q.times(D).times(p).scale(N.IOTA.neg()));
 prove('projection:connection_equals_negative_NJC',s.one().scale(2).times(p).minus(s.one()).times(d),K.scale(-1));
}
{
 const source=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-invariant-ratio-transport/CERTIFICATE.json'),'utf8')),spec=source.presentations.find(x=>x.name==='circle_with_reflection').spec;
 const{s:sys,w:{u,v,p,s}}=register('inherited_circle_reflection',spec),one=sys.one(),plus=one.plus(s).scale('1/2'),minus=one.minus(s).scale('1/2'),sin=u.minus(v).scale(N.Cut.of([0,'-1/2'])),sin2=u.times(u).minus(v.times(v)).scale(N.Cut.of([0,'-1/2']));
 for(const [name,odd]of [['sine',sin],['double_sine',sin2],['momentum',p]]){
  prove('reflection:odd_'+name,s.times(odd).plus(odd.times(s)),sys.zero());
  prove('reflection:maps_'+name,minus.times(odd).times(plus),odd.times(plus));
 }
 prove('reflection:even_kinetic_block',minus.times(p).times(p).times(plus),sys.zero());
}
const chosen=proofs.find(x=>x.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_word_proof',(()=>{try{systems.get(chosen.certificate.presentation).replay(bad);return false;}catch{return true;}})());
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),payload={protocol:'CLD1_MOVING_CUT_LOCAL_DYNAMICS_V1',status:'PASS_NATIVE_MOVING_CUT_CONNECTION_RETURN_AND_REFLECTION',scope:'Native finite algebra and fixtures; written finite analytic theorem and separately proved infinite AGD core identity. No unbounded memory norm theorem or physical quantum-gravity completion.',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results,summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
