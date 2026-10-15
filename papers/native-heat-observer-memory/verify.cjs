'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,m=N.matrix,mm=(...xs)=>xs.reduce(N.mul),dg=N.dagger,eq=N.equal,sc=N.scale;
const checks=[],negative_controls=[],proofs=[],presentations=[],systems=new Map();
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
const trace=A=>A.reduce((s,row,i)=>s.add(row[i]),N.ZERO),I2=N.identity(2),I3=N.identity(3);
const parent=JSON.parse(fs.readFileSync(path.join(__dirname,'../yang-mills-certified-benchmark/certificates/YM54_RESULT.json'),'utf8')).thermo_fixture;
const H=m([[2,1],[1,5]]),B=m([[2,1],[0,3]]),ds=[m([[0,1],[1,0]]),m([[1,0],[0,0]])],Hi=N.inverse(H),Bi=N.inverse(B);
ck('source:published_H',eq(H,m(parent.H)));
ck('source:factor',eq(mm(dg(B),B),sc(H,2)));
const vectors=ds.map(D=>{const X=mm(Hi,D),shape=N.sub(X,sc(I2,trace(X).div(2))),Y=mm(B,shape,Bi);ck('source:shape_self_adjoint:'+vectorsIndex(D),eq(Y,dg(Y)));return [Y[0][0],Y[0][1],N.ZERO];});
function vectorsIndex(D){return ds.indexOf(D);}
ck('source:published_vectors',eq(vectors,m(parent.vectors)));
const C=sc(mm(dg(vectors),vectors),'1/2'),Lambda=N.sub(sc(I3,trace(C)),C);
ck('source:published_protocol_tensor',eq(C,m(parent.C)));
ck('source:actual_heat_block',eq(Lambda,m([['5/72',0,0],[0,'5/162',0],[0,0,'65/648']])));
// Direct differentiation of the actual quaternion-coordinate quadratic functions.
const J=[m([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]]),m([[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]]),m([[0,0,0,-1],[0,0,-1,0],[0,1,0,0],[1,0,0,0]])];
const deriv=(F,a)=>{const T=sc(J[a],'-1/2');return N.add(mm(dg(T),F),mm(F,T));};
const forms=[m([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,-1]]),m([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]]),m([[0,0,-1,0],[0,0,0,1],[-1,0,0,0],[0,1,0,0]])];
for(let i=0;i<3;i++){
 ck('basis:centered:'+i,trace(forms[i]).zero());
 for(let j=0;j<3;j++)ck('basis:reference_Gram:'+i+':'+j,trace(mm(forms[i],forms[j])).div(12).eq(i===j?'1/3':0));
}
for(const [k,tensor]of [C,m([[2,1,0],[1,3,0],[0,0,1]])].entries()){
 const block=N.sub(sc(I3,trace(tensor)),tensor);
 for(let r=0;r<3;r++){
  let direct=N.zeros(4),expected=N.zeros(4);
  for(let a=0;a<3;a++)for(let b=0;b<3;b++)direct=N.sub(direct,sc(deriv(deriv(forms[r],b),a),tensor[a][b]));
  for(let j=0;j<3;j++)expected=N.add(expected,sc(forms[j],block[r][j]));
  ck('block:source_differential_intertwiner:'+k+':'+r,eq(direct,expected));
 }
}
const u=m([['3/5','4/5',0]]),v=m([['-4/5','3/5',0]]),w=m([[0,0,1]]),O=u.concat(v),G=sc(Lambda,-1),P=mm(dg(u),u),Q=N.sub(I3,P);
const Hobs=mm(O,Lambda,dg(O)),closure=N.rowClosure(u,[G]),descended=N.descendedAction(G,O);
ck('observer:orthonormal_rows',eq(mm(O,dg(O)),I2));
ck('observer:source_rates_in_retained_frame',eq(Hobs,m([['29/648','-1/54'],['-1/54','1/18']])));
ck('observer:minimal_rank_two',closure.rank===2&&N.rank(u)===1);
ck('observer:blind_third_direction',N.isZero(mm(closure.observer,dg(w))));
ck('observer:complete_repair_span',N.rank(O.concat(closure.observer))===2);
ck('observer:repaired_generator',eq(descended,sc(Hobs,-1)));
ck('observer:exact_intertwiner',eq(mm(O,G),mm(descended,O)));
no('observer:single_channel_autonomous',N.descendedAction(G,u)===null);
no('observer:projector_is_spectral',!N.isZero(N.commutator(P,Lambda)));
ck('leakage:initial_zero',N.isZero(mm(u,dg(v))));
ck('leakage:first_derivative',mm(u,G,dg(v))[0][0].eq('1/54'));
ck('defect:leading_coefficient',mm(u,Lambda,Q,Lambda,dg(u))[0][0].eq('1/2916'));
const l1=f('5/72'),l2=f('5/162'),delta=l1.sub(l2),weight=f('144/625');
ck('defect:source_rate_difference',delta.eq('25/648'));
ck('defect:variance_coefficient',weight.mul(delta.pow(2)).eq('1/2916'));
ck('memory:kernel_amplitude',Hobs[0][1].mul(Hobs[1][0]).eq('1/2916'));
ck('memory:decay_rate',Hobs[1][1].eq('1/18'));
ck('memory:integrated_tail_constant',f('1/2916').div('1/18').eq('1/162'));
ck('memory:hidden_initial_coefficient',Hobs[0][1].neg().eq('1/54'));
ck('memory:strict_positive_heat_determinant',Hobs[0][0].mul(Hobs[1][1]).sub(Hobs[0][1].mul(Hobs[0][1])).eq('25/11664'));
for(const z0 of ['1/2','1','3']){
 const z=f(z0),R=N.inverse(N.add(sc(I3,z),Lambda)),den=z.add('29/648').sub(f('1/2916').div(z.add('1/18')));
 ck('memory:full_source_resolvent_matches_kernel:'+z0,mm(u,R,dg(u))[0][0].eq(f(1).div(den)));
 ck('memory:hidden_source_transfer_matches:'+z0,mm(u,R,dg(v))[0][0].eq(f('1/54').div(z.add('1/18').mul(den))));
 no('memory:omit_kernel_changes_actual_resolvent:'+z0,!mm(u,R,dg(u))[0][0].eq(f(1).div(z.add('29/648'))));
}
// Alternating factorial-series enclosures; all coefficients are native exact rationals.
function expMinus(r){assert.ok(f(0).le(r)&&r.le(1));let term=f(1),sum=f(1),hi;for(let k=1;k<=21;k++){term=term.mul(r.neg()).div(k);sum=sum.add(term);if(k===20)hi=sum;}return {lo:sum,hi};}
const bounds=[];
for(const t0 of ['1','6','12']){
 const t=f(t0),fast=expMinus(l1.mul(t)),slow=expMinus(l2.mul(t)),gapLo=slow.lo.sub(fast.hi),gapHi=slow.hi.sub(fast.lo);
 const defectLo=weight.mul(gapLo.pow(2)),defectHi=weight.mul(gapHi.pow(2)),lower=t.pow(2).div(2916).mul(f(1).sub(l1.mul(t)).pow(2)),upper=t.pow(2).div(2916);
 ck('bounds:positive_gap:'+t0,f(0).le(gapLo)&&!gapLo.zero());
 ck('bounds:lower_budget:'+t0,lower.le(defectLo));
 ck('bounds:upper_budget:'+t0,defectHi.le(upper));
 ck('bounds:increasing_enclosure:'+t0,defectLo.le(defectHi));
 bounds.push({heat_ticks:t0,defect_lower:lower,defect_upper:upper,series_lower:defectLo,series_upper:defectHi});
}
ck('bounds:one_tick_lower',bounds[0].defect_lower.eq('4489/15116544'));
for(const row of [m([[1,0,0]]),m([[0,1,0]]),w])ck('controls:eigenchannel_closed:'+row.map(r=>r.map(z=>z.toJSON()).join(',')).join(''),N.descendedAction(G,row)!==null);
ck('controls:isotropic_channel_closed',N.descendedAction(sc(I3,-2),u)!==null);
no('controls:positive_heat_implies_every_observer_closes',!N.isZero(mm(Q,Lambda,dg(u))));
no('controls:zero_initial_readout_means_no_future_signal',!N.isZero(mm(u,G,dg(v))));
no('controls:discard_hidden_initial_data',!Hobs[0][1].zero());
no('controls:zero_memory_tail',!f('1/162').zero());
// HCM-2: exact phase-aware adapter; historical UGD operations are not replaced.
const alpha=N.Cut.of('-65/1296'),beta=new N.Cut('7/1296','24/1296');
const R=m([[0,-1],[1,0]]),K=m([[1,0],[0,-1]]),RK=mm(R,K),S=sc(N.add(sc(K,7),sc(RK,24)),'1/25');
const Ps=sc(N.add(I2,S),'1/2'),Pf=sc(N.sub(I2,S),'1/2'),sourceG=sc(Hobs,-1);
ck('phase:exact_EMK_generator',eq(sourceG,N.add(sc(I2,alpha),N.add(sc(K,'7/1296'),sc(RK,'24/1296')))));
ck('phase:reflection',eq(mm(S,S),I2));
ck('phase:complementary_spectral_sectors',N.isZero(mm(Ps,Pf))&&eq(N.add(Ps,Pf),I2));
ck('phase:slow_sector_rate',eq(mm(sourceG,Ps),sc(Ps,'-5/162')));
ck('phase:fast_sector_rate',eq(mm(sourceG,Pf),sc(Pf,'-5/72')));
no('phase:single_complex_scalar_generator',!N.isZero(N.commutator(sourceG,R)));
const phase=a=>[N.ONE,N.IOTA,N.ONE.neg(),N.IOTA.neg()][((a%4)+4)%4];
const pow2=n=>n>=0?f(2).pow(n):f(1).div(f(2).pow(-n));
const encode=z=>{
 if(z.zero())return {zero:true};
 const q=z.rad,h=z.turn,positive=v=>!v.zero()&&f(0).le(v);
 const a=positive(q)&&f(0).le(h)?0:positive(h)&&q.le(0)?1:positive(q.neg())&&h.le(0)?2:3;
 let w=z.mul(phase(-a)),n=0;
 while(w.norm2().le(1)&&!w.norm2().eq(1)){w=w.mul(2);n--;}
 while(f(4).le(w.norm2())){w=w.div(2);n++;}
 return {zero:false,a,n,w};
};
const decode=r=>r.zero?N.ZERO:phase(r.a).mul(pow2(r.n)).mul(r.w);
const complexOf=v=>new N.Cut(v[0][0].rad,v[1][0].rad);
const column=z=>m([[z.rad],[z.turn]]);
const samples=[[0,0],[1,0],[0,1],[-1,0],[0,-1],[1,1],[-1,2],[-3,-1],[2,-3],['1/16','1/32'],[16,24],[2,0]];
const phase_records=[];
for(const [j,pair]of samples.entries()){
 const z=new N.Cut(...pair),record=encode(z),dz=alpha.mul(z).add(beta.mul(z.dagger()));
 ck('record:round_trip:'+j,decode(record).eq(z));
 ck('phase:actual_source_intertwiner:'+j,complexOf(mm(sourceG,column(z))).eq(dz));
 if(!z.zero()){
  ck('record:canonical_sector:'+j,!record.w.rad.zero()&&f(0).le(record.w.rad)&&f(0).le(record.w.turn));
  ck('record:canonical_scale:'+j,f(1).le(record.w.norm2())&&record.w.norm2().le(4)&&!record.w.norm2().eq(4));
  const rate=dz.div(z),norm=z.norm2(),c2=z.rad.pow(2).sub(z.turn.pow(2)).div(norm),s2=z.rad.mul(z.turn).mul(2).div(norm);
  ck('phase:radial_rate:'+j,rate.rad.eq(f(-65).add(c2.mul(7)).add(s2.mul(24)).div(1296)));
  ck('phase:angular_rate:'+j,rate.turn.eq(c2.mul(24).sub(s2.mul(7)).div(1296)));
 }
 phase_records.push({input:z,record});
}
ck('record:zero_has_no_phase',encode(N.ZERO).zero&&!('a' in encode(N.ZERO)));
ck('record:scale_boundary_half_open',encode(N.Cut.of(2)).n===1&&encode(N.Cut.of(2)).w.eq(1));
const sameLabel1=encode(new N.Cut(1,0)),sameLabel2=encode(new N.Cut(1,1));
ck('record:equal_discrete_labels',sameLabel1.a===sameLabel2.a&&sameLabel1.n===sameLabel2.n);
no('record:labels_recover_residual',!sameLabel1.w.eq(sameLabel2.w));
no('record:visible_reading_recovers_hidden_phase',!complexOf(mm(sourceG,column(new N.Cut(1,0)))).rad.eq(complexOf(mm(sourceG,column(new N.Cut(1,1)))).rad));
ck('phase:nonzero_initial_angular_velocity',alpha.add(beta).turn.eq('1/54'));
for(let a=0;a<4;a++)for(const n of [-2,0,3]){
 const w0=new N.Cut('3/2','1/4'),frame=phase(a).mul(pow2(n)),z=frame.mul(w0),chartDerivative=alpha.mul(z).add(beta.mul(z.dagger())).div(frame);
 ck('chart:generator_intertwiner:'+a+':'+n,chartDerivative.eq(alpha.mul(w0).add(beta.mul(a%2?-1:1).mul(w0.dagger()))));
 const aNew=(a+1)%4,nNew=n+1,wNew=pow2(n-nNew),residual=w0.mul(wNew).mul(phase(a-aNew));
 ck('chart:transition_preserves_decoding:'+a+':'+n,phase(aNew).mul(pow2(nNew)).mul(residual).eq(z));
 if(a%2)no('chart:omit_conjugate_phase_factor:'+a+':'+n,!chartDerivative.eq(alpha.mul(w0).add(beta.mul(w0.dagger()))));
}
for(let a=0;a<4;a++)for(let b=0;b<4;b++)ck('UGD_label:cyclic_composition:'+a+':'+b,phase(a+b).eq(phase(a).mul(phase(b))));
// Source-pinned operation audit: mul_digits declares a*b mod K, not a+b mod K.
no('UGD_interface:numeral_product_equals_phase_composition',!phase(1*1).eq(phase(1).mul(phase(1))));
no('UGD_interface:neutral_digit_is_complex_zero',!phase(0).zero());
const zeroSource=new N.Cut(0,0);ck('phase:source_keeps_zero',alpha.mul(zeroSource).add(beta.mul(zeroSource.dagger())).zero());
function register(name,tokens){const rules=[];for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Central scalar coefficients in the retained heat block'});const s=new M.Presentation({tokens,rules}),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return s;}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
const ring=register('heat_phase_variables',['x','y']),x=ring.word(['x']),y=ring.word(['y']),a=x.scale('9/25').plus(y.scale('16/25'));
prove('words:exact_two_time_variance',x.times(x).scale('9/25').plus(y.times(y).scale('16/25')).minus(a.times(a)),x.times(x).minus(x.times(y).scale(2)).plus(y.times(y)).scale('144/625'));
const emk=new M.Presentation({tokens:['r','k'],rules:[{id:'rr',lhs:['r','r'],rhs:[[[],-1]],source:'Native circular generator'},{id:'kk',lhs:['k','k'],rhs:[[[],1]],source:'Native reflection'},{id:'kr',lhs:['k','r'],rhs:[[['r','k'],-1]],source:'Native anti-commutation'}]}),audit=emk.audit();
ck('EMK_source:confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name:'EMK_source_reflection',spec:emk.spec(),audit});systems.set(emk.hash(),emk);
const reflection=emk.word(['k']).scale(7).plus(emk.word(['r','k']).scale(24));
prove('words:source_reflection_square',reflection.times(reflection),emk.word([]).scale(625));
const chosen=proofs[0],bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_symbolic_witness',rejects(()=>systems.get(chosen.certificate.presentation).replay(bad)));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'HCM1_SOURCE_HEAT_OBSERVER_MEMORY_V1',status:'PASS_SOURCE_HEAT_OBSERVER_MEMORY',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results:{source_tensor:C,heat_block:Lambda,observer:u,retained_observer:O,retained_generator:descended,closure,hidden_initial_reading_derivative:'1/54',defect_leading_coefficient:'1/2916',memory_kernel:{amplitude:'1/2916',decay:'1/18',tail_prefactor:'1/162'},phase_law:{alpha,beta,reflection:S,slow_projector:Ps,fast_projector:Pf},phase_records,bounds,scope:{actual_source:'published YM54 free thermo fixture',observer:'declared mixed adjoint coefficient channel',full_interacting_row_defect_evaluated:false,physical_mass_gap:false,RH_proved:false,QG_complete:false,full_UGD_numeral_homomorphism:false,UGD_seam_identified_as_winding:false,source_UGD_implementation_rerun:false}},summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
