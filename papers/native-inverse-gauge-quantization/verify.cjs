'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const sum=a=>a.reduce((s,v)=>s.add(v),f(0));

// Actual source jets, keeping both orientations and independent gauge momenta.
for(const[i,row]of [[3,2,1,'1/2',2,-1],[-3,-2,2,1,-1,3],[2,3,'2/3',-1,4,2]].entries()){
 const[P,m,zp,v,K,R]=row.map(f),p=m.mul(zp),ap=p.add(v),am=p.sub(v),Fp=ap.div(P),Fm=am.div(P);
 ck('regular_embedding_derivatives_positive:'+i,f(0).le(Fp)&&f(0).le(Fm));
 const pip=P.mul(K).div(2).add(P.mul(R).div(m).div(2)),pim=P.mul(K).div(-2).add(P.mul(R).div(m).div(2));
 const H=p.mul(K).add(v.div(p).mul(zp).mul(R)),D=v.mul(K).add(zp.mul(R));
 ck('source_plus_cotangent_moment_map:'+i,H.add(D).div(2).eq(Fp.mul(pip)));
 ck('source_minus_cotangent_moment_map:'+i,H.sub(D).div(2).eq(Fm.mul(pim).neg()));
 ck('embedding_K_inverse:'+i,pip.sub(pim).div(P).eq(K));
 ck('embedding_R_inverse:'+i,pip.add(pim).mul(m).div(P).eq(R));
 const dz=f('2/5'),dphi=f('-3/4'),dFp=m.div(P).mul(dz).add(dphi.div(P)),dFm=m.div(P).mul(dz).sub(dphi.div(P));
 ck('fixed_m_cotangent_one_form:'+i,K.mul(dphi).add(R.mul(dz)).eq(pip.mul(dFp).add(pim.mul(dFm))));
 for(const eps of [1,-1]){
  const sm=f('3/2'),A=eps===1?ap:am,V=sm.mul(A).mul(eps).div(p).div(2),dZ=V.mul(zp),dph=sm.mul(A).div(2);
  const plus=m.div(P).mul(dZ).add(dph.div(P)),minus=m.div(P).mul(dZ).sub(dph.div(P));
  ck('embedding_chiral_active:'+i+':'+eps,(eps===1?plus:minus).eq(sm.mul(A).mul(eps).div(P)));
  ck('embedding_chiral_spectator:'+i+':'+eps,(eps===1?minus:plus).zero());
 }
 fixtures.push({kind:'source_embedding_jet',P,m,p,v,F_plus_prime:Fp,F_minus_prime:Fm,Pi_plus:pip,Pi_minus:pim});
}

// Mean-zero primitive is skew-adjoint. Finite Fourier witnesses use
// theta=2*pi*z and R_theta, so the common physical 1/(2*pi) is suppressed.
const clean=p=>new Map([...p].filter(([,v])=>!v.zero()));
const F=rows=>clean(new Map(rows.map(([n,v])=>[n,C(v)])));
const fm=(a,b)=>{const z=new Map();for(const[i,x]of a)for(const[j,y]of b)z.set(i+j,(z.get(i+j)||N.ZERO).add(x.mul(y)));return clean(z);};
const fscl=(a,c)=>clean(new Map([...a].map(([n,v])=>[n,v.mul(c)])));
const fpow=(a,n)=>{let r=F([[0,1]]);for(let k=0;k<n;k++)r=fm(r,a);return r;};
const Rtheta=a=>clean(new Map([...a].filter(([n])=>n!==0).map(([n,v])=>[n,v.div(C([0,n]))])));
const mean=a=>a.get(0)||N.ZERO;
const ff=F([[0,2],[1,'1/3'],[-1,'1/3']]),hh=F([[0,1],[1,[0,'1/2']],[-1,[0,'-1/2']]]);
ck('periodic_primitive_skew_adjoint',mean(fm(ff,Rtheta(hh))).eq(mean(fm(Rtheta(ff),hh)).neg()));
ck('periodic_primitive_removes_constant',mean(Rtheta(ff)).zero());
no('primitive_is_self_adjoint',!mean(fm(ff,Rtheta(hh))).eq(mean(fm(Rtheta(ff),hh))));
for(const[P,rm,phi,pis,pid,dRm,dP]of [[2,3,1,4,-2,1,3],[-3,1,-2,2,3,-1,2]].map(a=>a.map(f))){
 const A=pis.mul(dRm).div(P).sub(pis.mul(rm).add(pid.mul(phi)).mul(dP).div(P.pow(2)));
 const pplus=pis.add(pid).div(2),pminus=pis.sub(pid).div(2);
 const dFp=dRm.div(P).sub(rm.add(phi).mul(dP).div(P.pow(2))),dFm=dRm.div(P).sub(rm.sub(phi).mul(dP).div(P.pow(2)));
 ck('off_shell_physical_conjugate_correction:'+P,A.eq(pplus.mul(dFp).add(pminus.mul(dFm))));
 no('omit_off_shell_conjugate_correction:'+P,!A.zero());
}

// Brownian comparison: exact coefficient identities and summable bounds.
for(let n=1;n<=10;n++){
 const nn=f(n),base=nn.pow(2).add(1),r=nn.pow(4).div(base.pow(2)),rem=f(1).sub(r);
 ck('Brownian_coefficient_ratio:'+n,rem.eq(nn.pow(2).mul(2).add(1).div(base.pow(2))));
 ck('differentiable_remainder_bound:'+n,rem.le(f(3).div(nn.pow(2))));
 ck('coefficient_ratio_in_unit_interval:'+n,f(0).le(r)&&r.le(1));
 const bound=nn.add(1).pow(3).div(base.pow(2));
 ck('inverse_shift_Hilbert_Schmidt_basis_bound:'+n,bound.le(f(8).div(nn)));
 no('finite_mode_ratio_already_exact_Brownian:'+n,!r.eq(1));
}

// Spatial variation contrast: exact Fourier integration, not sampled paths.
for(const amp of [f(0),f('1/4'),f('1/2'),f('-1/2')]){
 const gp=F([[0,1],[1,amp.div(2)],[-1,amp.div(2)]]),qv=mean(fpow(gp,7));
 ck('circle_diffeomorphism_mean_derivative:'+amp,mean(gp).eq(1));
 if(amp.zero())ck('rotation_variation_control',qv.eq(1));
 else no('nonrigid_forward_measure_preserves_variation:'+amp,!qv.eq(1));
 if(amp.abs().eq('1/2'))ck('exact_nonrigid_variation_contrast:'+amp,qv.eq('4587/1024'));
 fixtures.push({kind:'deterministic_variation_multiplier',g_prime_fourier:[...gp],total_variation_multiplier:qv,baseline:1,scope:'exact integral of g_prime^7; stochastic singularity is a written theorem'});
}
for(const gprime of ['1/2',1,'3/2','4/3']){
 const a=f(gprime),amplitude=a.pow(3),qv=amplitude.pow(2).mul(a);
 ck('third_derivative_variation_chain:'+gprime,qv.eq(a.pow(7)));
 if(!a.eq(1))no('omit_parameter_change_factor:'+gprime,!amplitude.pow(2).eq(qv));
}

// Gaussian affinity fixture: a scalar covariance-ratio model, not a
// numerical estimate for the full diffeomorphism marginals in T2.
for(const ratio of [f(1),f(4),f('1/4')]){
 const affinity4=f(4).mul(ratio).div(ratio.add(1).pow(2));
 for(let n=1;n<=5;n++){
  ck('finite_Gaussian_affinity_positive:'+ratio+':'+n,f(0).le(affinity4.pow(n)));
  if(ratio.eq(1))ck('equal_covariance_affinity:'+n,affinity4.pow(n).eq(1));
  else ck('scalar_covariance_affinity_decay:'+ratio+':'+n,affinity4.pow(n).eq(f('16/25').pow(n)));
 }
}

// Inverse-embedding variation has the normalization term required by (17).
const weights=['1/4','1/4','1/2'].map(f),density=[1,2,'1/2'].map(f),dir=[2,-1,3].map(f);
ck('normalized_embedding_density',sum(weights.map((w,i)=>w.mul(density[i]))).eq(1));
const mh=sum(weights.map((w,i)=>w.mul(density[i]).mul(dir[i]))),dk=density.map((v,i)=>v.mul(dir[i].sub(mh)));
ck('inverse_embedding_variation_has_zero_period',sum(weights.map((w,i)=>w.mul(dk[i]))).zero());
no('drop_normalization_from_inverse_variation',!sum(weights.map((w,i)=>w.mul(density[i]).mul(dir[i]))).zero());
fixtures.push({kind:'weighted_normalization_derivative',scope:'finite derivative fixture only',mean_direction:mh,density_variation:dk});

// Exact Radon--Nikodym square-root representation in a weighted finite space.
const roots=[1,2,3,4].map(f),wt=roots.map(x=>x.pow(2).div(30));
const gram=N.matrix(wt.map((w,i)=>wt.map((_,j)=>i===j?w:0)));
const ident=[0,1,2,3],g=[1,2,3,0],h=[1,0,2,3],compose=(p,q)=>q.map(x=>p[x]);
const inverse=p=>p.map((_,i)=>p.indexOf(i));
const repr=perm=>N.matrix(roots.map((ri,i)=>roots.map((rj,j)=>perm[j]===i?rj.div(ri):f(0))));
for(const[name,perm]of [['e',ident],['g',g],['h',h],['gh',compose(g,h)]]){
 const U=repr(perm);
 ck('weighted_gauge_unitary:'+name,N.equal(N.mul(N.mul(N.dagger(U),gram),U),gram));
 ck('weighted_gauge_inverse:'+name,N.equal(N.mul(U,repr(inverse(perm))),N.identity(4)));
}
ck('exact_gauge_group_law',N.equal(N.mul(repr(g),repr(h)),repr(compose(g,h))));
no('reverse_noncommuting_group_order',!N.equal(N.mul(repr(g),repr(h)),repr(compose(h,g))));
const gp=[1,0,3,2],gm=[2,3,0,1];
ck('independent_gauge_actions_commute',N.equal(N.mul(repr(gp),repr(gm)),N.mul(repr(gm),repr(gp))));
const plain=N.matrix(roots.map((_,i)=>roots.map((_,j)=>g[j]===i?1:0)));
no('omit_square_root_density_unitarity',!N.equal(N.mul(N.mul(N.dagger(plain),gram),plain),gram));
const wrong=N.matrix(roots.map((ri,i)=>roots.map((rj,j)=>g[j]===i?rj.div(ri).pow(2):f(0))));
no('use_density_instead_of_square_root',!N.equal(N.mul(N.mul(N.dagger(wrong),gram),wrong),gram));

// Relative-rotation chart and Gaussian half-density cancellation.
for(const[beta,delta,t]of [[1,2,3],['1/3','-2/5','3/4'],[-2,'1/2',1]].map(row=>row.map(f))){
 ck('relative_rotation_invariant_coordinate:'+beta,beta.add(t).add(delta.sub(t)).eq(beta.add(delta)));
 const exponent=x=>x.pow(2).div(-2),logFactor=exponent(delta.add(t)).sub(exponent(delta)).div(2);
 ck('relative_rotation_density_exponent:'+beta,logFactor.eq(delta.mul(t).div(-2).sub(t.pow(2).div(4))));
 ck('flat_gauge_wave_half_density_cancellation:'+beta,logFactor.sub(exponent(delta.add(t)).div(2)).eq(exponent(delta).div(-2)));
}

// General word identities. They certify the displayed finite algebra only.
const tokens=['p','m','z','v','k','r','f','g','a','b','B'],rules=[];
for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting source, smearing and variation jets'});
const alg=new M.Presentation({tokens,rules}),audit=alg.audit();ck('word_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:alg.spec(),audit});
function prove(name,lhs,rhs){const out=alg.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',alg.replay(cert));proofs.push({name,presentation_sha256:alg.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
const[p,m,z,v,k,r,ff0,gg,a,b,BB]=tokens.map(t=>alg.word([t]));
prove('source_plus_moment_map_cleared',m.times(m.times(z).plus(v)).times(k).plus(m.times(z).plus(v).times(r)),m.times(z).plus(v).times(m.times(k).plus(r)));
prove('source_minus_moment_map_cleared',m.times(m.times(z).minus(v)).times(k).minus(m.times(z).minus(v).times(r)),m.times(z).minus(v).times(m.times(k).minus(r)));
prove('cotangent_one_form_fixed_m',k.times(a).plus(r.times(b)),k.plus(r).times(a.plus(b)).plus(r.minus(k).times(b.minus(a))).scale('1/2'));
prove('physical_conjugate_off_shell_sign',m.times(a).scale(-1).plus(p.times(b)).minus(BB.times(a)),m.plus(BB).times(a).scale(-1).plus(p.times(b)));
prove('left_vector_field_Lie_coefficient',ff0.times(b.times(p).plus(gg.times(r))).minus(gg.times(a.times(p).plus(ff0.times(r)))),ff0.times(b).minus(gg.times(a)).times(p));
const bad=JSON.parse(JSON.stringify(proofs.find(x=>x.certificate.steps.length).certificate));bad.steps[0].coefficient=['991','0'];
no('corrupted_rewrite_certificate',(()=>{try{alg.replay(bad);return false;}catch{return true;}})());
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'IGQ1_INVERSE_GAUGE_QUANTIZATION_V1',status:'PASS_INVERSE_GAUGE_UNITARIES_WITH_RECORDED_RIGGING_BOUNDARY',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{source_two_embedding_canonical_chart:true,off_shell_physical_conjugate_correction_retained:true,forward_Gaussian_gauge_lift:'mutually singular for every non-rigid diffeomorphism',variation_witness:'4587/1024 versus 1',inverse_embedding_measure_quasi_invariant:true,finite_gauge_unitary_representation:true,common_integer_identification_retained:true,individual_flow_generators_self_adjoint:true,full_common_infinitesimal_domain_constructed:false,kinematic_invariant_subspace:'zero; not generalized nonexistence',positive_rotation_subgroup_rigging:true,full_gauge_rigging_constructed:false,MPQ_physical_equivalence:'conditional on one-channel positive full gauge rigging',new_probe_representative:'differs from original by the classical constraint ideal; quantum ordering choice explicit',source_action_changed:false,physical_three_plus_one_gravity_closed:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cpth=path.join(__dirname,'CERTIFICATE.json'),hpth=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cpth,bytes);fs.writeFileSync(hpth,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cpth,'utf8'),bytes);assert.equal(fs.readFileSync(hpth,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
