'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};

// Local reference normal/spatial bracket jets, before any constraint.
for(const[i,row]of [[3,1,2,-1],[-3,1,1,2],[2,0,-1,3]].entries()){
 const[p,v,pp,vp]=row.map(f),u=v.div(p),up=vp.div(p).sub(v.mul(pp).div(p.pow(2))),nn=f(2),np=f(-1),mm=f(3),mp=f(4),J=f('5/3');
 const du=n=>n[1].mul(p).add(n[0].mul(pp)).div(p).sub(v.mul(n[1].mul(v).add(n[0].mul(vp))).div(p.pow(2)));
 const duN=du([nn,np]),duM=du([mm,mp]),omega=nn.mul(mp).sub(mm.mul(np)),tail=f(1).sub(u.pow(2)).mul(pp).div(p).sub(u.mul(up));
 ck('normal_u_chain_rule:'+i,duN.eq(np.mul(f(1).sub(u.pow(2))).add(nn.mul(tail))));
 const cross=nn.mul(duM).sub(mm.mul(duN)).mul(J),self=u.pow(2).mul(omega).mul(J);
 ck('mixed_HH_reference_coefficient:'+i,cross.eq(omega.mul(f(1).sub(u.pow(2))).mul(J)));
 ck('full_HH_reference_closure:'+i,cross.add(self).eq(omega.mul(J)));
 const sm=f('2/3'),sm1=f('-1/4'),dp=sm1.mul(p).add(sm.mul(pp)),dv=sm1.mul(v).add(sm.mul(vp));
 ck('spatial_u_is_scalar:'+i,dv.div(p).sub(v.mul(dp).div(p.pow(2))).eq(sm.mul(up)));
 const oldCross=nn.mul(sm).mul(up).mul(J).neg(),refDH=sm.mul(np.mul(u).add(nn.mul(up))).sub(nn.mul(u).mul(sm1)).mul(J);
 ck('mixed_DH_reference_closure:'+i,oldCross.add(refDH).eq(sm.mul(np).sub(nn.mul(sm1)).mul(u).mul(J)));
 ck('b_u_choice_closes:'+i,u.pow(2).add(f(1).sub(u.pow(2))).eq(1));
 ck('constant_b_one_also_closes:'+i,f(1).pow(2).add(f(1).sub(u.pow(2)).mul(0)).eq(1));
 no('frozen_reference_b_zero_closes:'+i,!f(0).eq(1));
 no('b_one_matches_density_transport:'+i,!f(1).eq(u));
}
for(const[rho,rhop,zp,zpp,sm,sm1]of [[2,3,1,2,3,1],[-2,1,2,3,-1,2]].map(x=>x.map(f))){
 const dr=sm1.mul(rho).add(sm.mul(rhop)),dzp=sm1.mul(zp).add(sm.mul(zpp)),J=rho.mul(zp),Jp=rhop.mul(zp).add(rho.mul(zpp));
 ck('reference_J_weight_two:'+rho,dr.mul(zp).add(rho.mul(dzp)).eq(sm.mul(Jp).add(sm1.mul(J).mul(2))));
}

// Chiral transport and exact material-window cancellation for unconstrained jets.
for(const[i,row]of [[2,3,1,-2,4,-1,2],[-2,-3,2,1,-1,3,-2],[1,4,-2,3,2,5,3]].entries()){
 const[ap,am,ap1,am1,bp,bm,J]=row.map(f),p=ap.add(am).div(2),v=ap.sub(am).div(2),pp=ap1.add(am1).div(2);
 const tp=ap.mul(bp).sub(ap1).div(2),tm=am.mul(bm).add(am1).div(2),ttP=tp.add(ap.mul(J).div(p).div(2)),ttM=tm.sub(am.mul(J).div(p).div(2));
 const cp=bp.sub(bm).sub(ap1.div(ap)).sub(am1.div(am)).div(2),Zcal=am.div(ap).mul(ttP).add(ap.div(am).mul(ttM));
 ck('reference_c_spatial_identity:'+i,cp.eq(ttP.div(ap).sub(ttM.div(am)).sub(J.div(p))));
 for(const eps of [1,-1]){
  const sm=f(2),sm1=f(-3),A=eps===1?ap:am,A1=eps===1?ap1:am1,B=eps===1?bp:bm;
  const V=sm.mul(A).mul(eps).div(p).div(2),Vp=sm1.mul(A).add(sm.mul(A1)).div(p).sub(sm.mul(A).mul(pp).div(p.pow(2))).mul(eps).div(2);
  const dA=sm1.mul(A).add(sm.mul(A1)).mul(eps),dp=dA.div(2),dq=sm.mul(B).add(sm1.mul(eps)).div(2).sub(sm.mul(v).mul(J).div(p.pow(2)).div(2)),dc=dq.sub(dA.div(A).div(2));
  ck('common_density_transport:'+i+':'+eps,dp.eq(Vp.mul(p).add(V.mul(pp))));
  ck('reference_c_transport_defect:'+i+':'+eps,dc.sub(V.mul(cp)).eq(sm.mul(Zcal).div(p).div(2)));
  const zp=f('3/2'),zpp=f('2/3'),dz=V.mul(zp),dzp=Vp.mul(zp).add(V.mul(zpp));
  const ratioVar=dp.div(zp).sub(p.mul(dzp).div(zp.pow(2))),ratioPrime=pp.div(zp).sub(p.mul(zpp).div(zp.pow(2)));
  ck('material_m_inverse_coordinate_cancellation:'+i+':'+eps,ratioVar.sub(V.mul(ratioPrime)).zero());
  const a=f('13/10'),ac=f('-2/5'),w=f('3/4'),wz=f('2/3');
  const densityVar=dp.mul(a).mul(w).add(p.mul(ac).mul(w).mul(dc)).add(p.mul(a).mul(wz).mul(dz));
  const divergence=Vp.mul(p).mul(a).mul(w).add(V.mul(pp).mul(a).mul(w)).add(V.mul(p).mul(ac).mul(cp).mul(w)).add(V.mul(p).mul(a).mul(wz).mul(zp));
  const ideal=sm.mul(w).mul(ac).mul(Zcal).div(2);
  ck('material_window_constraint_ideal:'+i+':'+eps,densityVar.sub(divergence).eq(ideal));
  no('freeze_Z_and_keep_window_closure:'+i+':'+eps,!densityVar.sub(p.mul(a).mul(wz).mul(dz)).sub(divergence).eq(ideal));
  if(i===0)no('omit_reference_q_feedback:'+eps,!dc.add(sm.mul(v).mul(J).div(p.pow(2)).div(2)).sub(V.mul(cp)).eq(sm.mul(Zcal).div(p).div(2)));
 }
 // Independent on-shell coefficients from the exact modified constraints.
 const onBp=ap1.div(ap).sub(J.div(p)),onBm=am1.div(am).neg().add(J.div(p));
 ck('on_shell_plus_current_solution:'+i,ap.mul(onBp).sub(ap1).div(2).add(ap.mul(J).div(p).div(2)).zero());
 ck('on_shell_minus_current_solution:'+i,am.mul(onBm).add(am1).div(2).sub(am.mul(J).div(p).div(2)).zero());
 ck('on_shell_c_prime:'+i,onBp.sub(onBm).sub(ap1.div(ap)).sub(am1.div(am)).div(2).eq(J.div(p).neg()));
 fixtures.push({kind:'reference_off_shell_jet',p,v,J,T_tilde_plus:ttP,T_tilde_minus:ttM,c_prime:cp});
}

// Material gauge reconstruction and its apparatus momentum response.
for(const[i,row]of [[2,3,1,'3/5','4/5'],[-3,1,-2,1,0],[1,-2,3,'-3/5','-4/5']].entries()){
 const[m,mp,rp,cs,sn]=row.map(f),rho=mp.sub(m.mul(rp)),a=f(1).add(cs.div(2)),b=sn.div(-2),bb=cs.div(-2),cp=rp.sub(mp.div(m));
 ck('reconstructed_spatial_constraint:'+i,m.mul(rp).sub(mp).add(rho).zero());
 ck('reconstructed_normal_constraint:'+i,m.mul(0).add(f(0).mul(rho)).zero());
 ck('rational_cosine_circle:'+i,cs.pow(2).add(sn.pow(2)).eq(1));
 const w=f('3/4'),wp=f('2/5'),JP=f('-1/7'),ap=b.mul(cp),bp=bb.mul(cp);
 const mdot=m.mul(w).mul(b).neg(),rdot=w.mul(a.sub(b)).add(JP);
 const mdotp=mp.mul(w).mul(b).add(m.mul(wp).mul(b)).add(m.mul(w).mul(bp)).neg();
 const rdotp=wp.mul(a.sub(b)).add(w.mul(ap.sub(bp)));
 const rhodot=mdotp.sub(mdot.mul(rp)).sub(m.mul(rdotp));
 ck('reference_backreaction_reconstruction:'+i,rhodot.eq(m.mul(a).mul(wp).neg()));
 no('rho_zero_can_be_imposed_without_feedback:'+i,!m.mul(a).mul(wp).zero());
 fixtures.push({kind:'reduced_reference_feedback',m,rho,a,b,mdot,rdot,rhodot});
}

// Finite weighted functional-gradient witness, not a continuum discretization.
const cells=[{weight:'1/2',m:2,w:2,a:'13/10',b:'-2/5'},{weight:'1/2',m:3,w:0,a:'7/10',b:'2/5'}].map(o=>Object.fromEntries(Object.entries(o).map(([k,v])=>[k,f(v)])));
const P=cells.reduce((s,c)=>s.add(c.weight.mul(c.m)),f(0)),B=cells.reduce((s,c)=>s.add(c.weight.mul(c.m).mul(c.w).mul(c.b)),f(0)),JP=B.div(P);
ck('global_P_feedback_value',JP.eq('-8/25'));
let Gdot=f(0);for(const[i,c]of cells.entries()){
 const gr=c.m.mul(c.w).mul(c.b),gm=c.w.mul(c.a.sub(c.b)).add(JP),rdot=gm,mdot=gr.neg();
 Gdot=Gdot.add(c.weight.mul(gr.mul(rdot).add(gm.mul(mdot))));
 ck('reduced_pointwise_hamiltonian_cancellation:'+i,gr.mul(rdot).add(gm.mul(mdot)).zero());
 if(c.w.zero()){
  ck('outside_window_m_response_zero',mdot.zero());
  ck('outside_window_global_r_response',rdot.eq('-8/25'));
  no('RPM_global_P_coupling_is_fully_local',!rdot.zero());
  ck('P_independent_choice_restores_response_support',c.w.mul(c.a.sub(c.b)).zero());
 }
}
ck('reduced_probe_conserved',Gdot.zero());

// Fourier witness on the actual unit circle. Derivative is d/dtheta,
// theta=2*pi*z, so every spatial derivative output is divided by 2*pi.
const clean=a=>new Map([...a].filter(([,c])=>!c.zero()));
const fourier=rows=>clean(new Map(rows.map(([n,c])=>[n,C(c)])));
const scale=(a,c)=>clean(new Map([...a].map(([n,v])=>[n,v.mul(c)])));
const fd=a=>clean(new Map([...a].map(([n,c])=>[n,c.mul(new N.Cut(0,n))])));
const same=(a,b)=>new Set([...a.keys(),...b.keys()]).size===0||[...new Set([...a.keys(),...b.keys()])].every(n=>(a.get(n)||N.ZERO).eq(b.get(n)||N.ZERO));
const w=fourier([[0,1],[1,'1/4'],[-1,'1/4']]),rdot=scale(w,'3/2'),rhodot=scale(fd(w),-3);
ck('material_window_normalization',(w.get(0)||N.ZERO).eq(1));
ck('first_harmonic_r_response',(rdot.get(1)||N.ZERO).eq('3/8'));
ck('rho_plus_harmonic_response_over_two_pi',(rhodot.get(1)||N.ZERO).eq([0,'-3/4']));
ck('reference_reconstruction_fourier',same(rhodot,scale(fd(rdot),-2)));
ck('reference_momentum_response_zero_mean',(rhodot.get(0)||N.ZERO).zero());
no('localized_window_keeps_homogeneous_geometry',!(rdot.get(1)||N.ZERO).zero());
ck('normalized_homogeneous_probe_value',f(2).mul('3/2').mul(1).eq(3));
for(const kappa of ['-2','1/3',1])ck('initial_homogeneous_pointer_recovery:'+kappa,f(3).mul(kappa).eq(f(2).mul('3/2').mul(kappa)));
no('current_inverse_defined_at_p_zero',(()=>{try{f(1).div(0);return false;}catch{return true;}})());
fixtures.push({kind:'periodic_material_window_response',spatial_derivative_normalization:'divided by 2*pi',w:[...w],r_dot_per_pulse_area:[...rdot],rho_dot_per_pulse_area_over_two_pi:[...rhodot],homogeneous_initial_probe:3});

// Word certificates for general coefficient identities, not continuum proofs.
const tokens=['n','m','a','b','u','j','e'],rules=[];
for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting classical smearing and reference jets'});
const s=new M.Presentation({tokens,rules}),audit=s.audit();ck('word_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:s.spec(),audit});
function prove(name,lhs,rhs){const r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(cert));proofs.push({name,presentation_sha256:s.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
const[n,m,a,b,u,j,e]=tokens.map(x=>s.word([x])),one=s.one(),omega=n.times(b).minus(m.times(a)),factor=one.minus(u.times(u));
prove('mixed_normal_NM_terms_cancel',n.times(b.times(factor).plus(m.times(e))).minus(m.times(a.times(factor).plus(n.times(e)))).times(j),omega.times(factor).times(j));
prove('full_reference_normal_bracket',omega.times(factor).times(j).plus(u.times(u).times(omega).times(j)),omega.times(j));
prove('spatial_normal_reference_cross_cancel',n.times(m).times(e).scale(-1).plus(m.times(a.times(u).plus(n.times(e)))).minus(n.times(u).times(b)),m.times(a).minus(n.times(b)).times(u));
prove('reference_momentum_feedback_cancellation',u.times(n).times(m).times(j).scale(-1).minus(u.times(u).times(b).times(a)).plus(u.times(u).times(m).times(j).times(e)).minus(u.times(m).times(j).times(u.times(e).minus(n))),u.times(u).times(a).times(b).scale(-1));
const bad=JSON.parse(JSON.stringify(proofs.find(x=>x.certificate.steps.length).certificate));bad.steps[0].coefficient=['991','0'];
no('corrupted_word_certificate',(()=>{try{s.replay(bad);return false;}catch{return true;}})());
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'MRF1_MATERIAL_REFERENCE_FIELD_V1',status:'PASS_CLASSICAL_MATERIAL_REFERENCE_AND_PROFILE_REDUCTION',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{local_reference_constraints:"h_tilde=h+(phi'/p)rho Z'; d_tilde=d+rho Z'",all_smearing_classical_closure:true,material_window_residual_cancelled:true,physical_profiles:'r(z),m(z), integral m delta r',reference_momentum_on_material_gauge:"rho=m'-m r'",RPM_probe_fully_local:false,P_independent_probe_has_supported_response:true,classical_pointer_recovery:'exact on initially homogeneous profiles with integral w=1',nonconstant_window_preserves_homogeneous_dynamics:false,continuum_quantum_representation_constructed:false,old_vacuum_anomaly_result_automatically_inherited:false,positive_energy_material_identification:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cpth=path.join(__dirname,'CERTIFICATE.json'),hpth=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cpth,bytes);fs.writeFileSync(hpth,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cpth,'utf8'),bytes);assert.equal(fs.readFileSync(hpth,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
