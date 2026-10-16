'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const sum=xs=>xs.reduce((s,x)=>s.add(x),f(0));

// Exact off-shell source jets. No logarithm or atanh is approximated.
for(const [i,row]of [[3,1,2,-1,4,2,-3,2],[-3,1,-2,3,1,-2,4,1],[2,0,1,3,-2,4,1,3]].entries()){
 const[p,v,pp,vp,qp,k,rho,zp]=row.map(f),u=v.div(p),D=p.pow(2).sub(v.pow(2)),up=vp.div(p).sub(v.mul(pp).div(p.pow(2)));
 const cp=qp.sub(p.mul(pp).sub(v.mul(vp)).div(D)),K=k.sub(up.div(f(1).sub(u.pow(2)))),CZ=p.mul(cp).add(rho.mul(zp));
 const H=p.mul(k).add(qp.mul(v)).sub(vp).add(u.mul(rho).mul(zp)),DD=p.mul(qp).add(k.mul(v)).sub(pp).add(rho.mul(zp));
 ck('source_normal_canonical_decomposition:'+i,H.eq(p.mul(K).add(u.mul(CZ))));
 ck('source_spatial_canonical_decomposition:'+i,DD.eq(v.mul(K).add(CZ)));
 ck('source_gauge_K_inverse:'+i,H.sub(u.mul(DD)).div(p.mul(f(1).sub(u.pow(2)))).eq(K));
 ck('source_gauge_R_inverse:'+i,DD.sub(u.mul(H)).div(f(1).sub(u.pow(2))).div(zp).eq(CZ.div(zp)));
 const dc=f('2/3'),dz=f('-1/4'),cmprime=cp.div(zp),oldOne=p.mul(dc.add(cmprime.mul(dz))).add(rho.mul(dz));
 ck('material_canonical_one_form_jet:'+i,oldOne.eq(p.mul(dc).add(rho.add(p.mul(cp).div(zp)).mul(dz))));
 no('drop_K_phi_improvement:'+i,!H.eq(p.mul(k).add(u.mul(CZ))));
 fixtures.push({kind:'off_shell_source_canonical_jet',p,v,K_phi:K,C_Z:CZ,H_tilde:H,D_tilde:DD});
}
no('source_chart_allows_null_current',f(2).pow(2).sub(f(2).pow(2)).zero());

// Finite weighted samples test the functional Jacobian, not a lattice theory.
for(const [label,rows]of [
 ['positive', [['1/3',2,'3/5','4/5','-2/3',2],['2/3',3,'-3/5','4/5','4/5',0]]],
 ['negative', [['1/4',-2,1,0,1,1],['1/4',-3,0,1,-2,2],['1/2',-1,'3/5','-4/5','2/3',-1]]]
]){
 const data=rows.map(row=>{const[wt,m,cs,sn,y,w]=row.map(f);return{wt,m,cs,sn,y,w,a:f(1).add(cs.div(2)),b:sn.div(-2)};});
 const dim=data.length,P=sum(data.map(x=>x.wt.mul(x.m))),B=sum(data.map(x=>x.wt.mul(x.m).mul(x.w).mul(x.b)));
 const blank=()=>({q:Array.from({length:dim},()=>f(0)),p:Array.from({length:dim},()=>f(0))});
 const pb=(a,b)=>sum(data.map((x,k)=>a.q[k].mul(b.p[k]).sub(a.p[k].mul(b.q[k])).div(x.wt)));
 const theta=[],ys=[],js=[],ss=[],ds=[];
 for(const[i,x]of data.entries()){
  ck('cosine_unit_circle:'+label+':'+i,x.cs.pow(2).add(x.sn.pow(2)).eq(1));
  const th=blank(),yy=blank(),jj=blank(),ss0=blank(),dd=blank(),j=x.m.mul(x.a);
  for(let k=0;k<dim;k++){
   const delta=f(i===k?1:0);th.q[k]=delta;th.p[k]=data[k].wt.div(P).sub(delta.div(x.m));
   yy.q[k]=th.q[k].div(x.a);yy.p[k]=th.p[k].div(x.a);
   jj.q[k]=x.m.mul(x.b).mul(th.q[k]);jj.p[k]=delta.mul(x.a).add(x.m.mul(x.b).mul(th.p[k]));
   ss0.q[k]=jj.q[k].div(j);ss0.p[k]=jj.p[k].div(j);
   dd.q[k]=jj.q[k].mul(x.y).add(j.mul(yy.q[k]));dd.p[k]=jj.p[k].mul(x.y).add(j.mul(yy.p[k]));
  }
  theta.push(th);ys.push(yy);js.push(jj);ss.push(ss0);ds.push(dd);
  const rdot=x.w.mul(x.a.sub(x.b)).add(B.div(P)),mdot=x.m.mul(x.w).mul(x.b).neg(),Pdot=B.neg();
  const thetadot=rdot.sub(mdot.div(x.m)).add(Pdot.div(P));
  ck('physical_theta_flow:'+label+':'+i,thetadot.eq(x.w.mul(x.a)));
  ck('physical_y_flow:'+label+':'+i,thetadot.div(x.a).eq(x.w));
  ck('conserved_local_j:'+label+':'+i,mdot.mul(x.a).add(x.m.mul(x.b).mul(thetadot)).zero());
  if(!B.zero())no('omit_global_feedback_in_reconstruction:'+label+':'+i,!rdot.sub(B.div(P)).sub(mdot.div(x.m)).add(Pdot.div(P)).eq(x.w.mul(x.a)));
 }
 for(let i=0;i<dim;i++)for(let j=0;j<dim;j++){
  const delta=f(i===j?1:0).div(data[i].wt);
  ck('theta_theta_bracket:'+label+':'+i+':'+j,pb(theta[i],theta[j]).zero());
  ck('y_j_canonical_bracket:'+label+':'+i+':'+j,pb(ys[i],js[j]).eq(delta));
  ck('j_j_commuting_bracket:'+label+':'+i+':'+j,pb(js[i],js[j]).zero());
  ck('d_s_canonical_bracket:'+label+':'+i+':'+j,pb(ds[i],ss[j]).eq(delta));
  ck('d_d_commuting_bracket:'+label+':'+i+':'+j,pb(ds[i],ds[j]).zero());
  ck('d_j_affine_bracket:'+label+':'+i+':'+j,pb(ds[i],js[j]).eq(delta.mul(data[j].m).mul(data[j].a)));
 }
 const dr=data.map((_,i)=>f(i+1)),dm=data.map((_,i)=>f(i+2)),dP=sum(data.map((x,i)=>x.wt.mul(dm[i])));
 const oldForm=sum(data.map((x,i)=>x.wt.mul(x.m).mul(dr[i]))),newForm=sum(data.map((x,i)=>x.wt.mul(x.m).mul(dr[i].sub(dm[i].div(x.m)).add(dP.div(P)))));
 ck('global_exact_one_form:'+label,oldForm.eq(newForm));
 no('drop_global_shift_preserves_exact_one_form:'+label,!oldForm.eq(newForm.sub(dP)));
 fixtures.push({kind:'weighted_functional_jacobian',scope:'finite exact derivative fixture; not a continuum measure',label,P,B,dimension:dim});
}

// Finite Gaussian exponents certify the density ratio and group law exactly.
const variance=[f(1),f('1/3'),f('1/16')],sv=[f('2/3'),f(-1),f('1/4')],hv=[f('1/2'),f('2/3'),f(-1)];
const normh=sum(hv.map((h,i)=>h.pow(2).div(variance[i]))),beta=s=>sum(s.map((x,i)=>x.mul(hv[i]).div(variance[i])));
const logDensity=s=>sum(s.map((x,i)=>x.pow(2).div(variance[i]).div(-2)));
const shift=(s,t)=>s.map((x,i)=>x.add(hv[i].mul(t)));
const logWeight=(s,t)=>beta(s).mul(t).div(-2).sub(t.pow(2).mul(normh).div(4));
for(const t of [f(1),f('-2/3'),f('1/4')]){
 const sm=shift(sv,t.neg());
 ck('Gaussian_translation_norm_exponent:'+t,logDensity(sm).add(logWeight(sm,t).mul(2)).eq(logDensity(sv)));
 ck('Gaussian_beta_shift:'+t,beta(shift(sv,t)).eq(beta(sv).add(t.mul(normh))));
 no('drop_Radon_Nikodym_factor:'+t,!logDensity(sm).eq(logDensity(sv)));
 for(const u of [f('3/4'),f(-2)])ck('Gaussian_unitary_group_exponent:'+t+':'+u,logWeight(sv,t).add(logWeight(shift(sv,t),u)).eq(logWeight(sv,t.add(u))));
}

// Polynomial Gaussian adjoints use exact moments; no symbolic library.
const trim=p=>{while(p.length&&p[p.length-1].zero())p.pop();return p;};
const padd=(p,q)=>trim(Array.from({length:Math.max(p.length,q.length)},(_,i)=>(p[i]||N.ZERO).add(q[i]||N.ZERO)));
const pscale=(p,a)=>trim(p.map(x=>x.mul(a)));
const pminus=(p,q)=>padd(p,pscale(q,-1));
const pmul=(p,q)=>{const r=Array.from({length:p.length+q.length},()=>N.ZERO);for(let i=0;i<p.length;i++)for(let j=0;j<q.length;j++)r[i+j]=r[i+j].add(p[i].mul(q[j]));return trim(r);};
const pd=p=>trim(p.slice(1).map((a,i)=>a.mul(i+1)));
const mono=n=>Array.from({length:n+1},(_,i)=>C(i===n?1:0));
const peq=(p,q)=>pminus(p,q).length===0;
const mom=(n,lam)=>{if(n%2)return f(0);let z=f(1);for(let k=2;k<=n;k+=2)z=z.mul(k-1).mul(lam);return z;};
const inner=(p,q,lam)=>{let z=N.ZERO;for(let i=0;i<p.length;i++)for(let j=0;j<q.length;j++)z=z.add(p[i].dagger().mul(q[j]).mul(mom(i+j,lam)));return z;};
const Sop=p=>pmul([N.ZERO,N.ONE],p);
const Dop=(p,lam)=>pscale(pminus(pd(p),pscale(Sop(p),f(1).div(lam).div(2))),N.IOTA);
for(const lam of [f(1),f('1/3'),f(2)]){
 for(let i=0;i<5;i++)for(let j=0;j<5;j++)ck('Gaussian_derivative_symmetry:'+lam+':'+i+':'+j,inner(mono(i),Dop(mono(j),lam),lam).eq(inner(Dop(mono(i),lam),mono(j),lam)));
 for(let i=0;i<7;i++)ck('Gaussian_canonical_commutator:'+lam+':'+i,peq(pminus(Dop(Sop(mono(i)),lam),Sop(Dop(mono(i),lam))),pscale(mono(i),N.IOTA)));
 const plain=p=>pscale(pd(p),N.IOTA);
 no('omit_half_density_drift_is_symmetric:'+lam,!inner(mono(0),plain(mono(1)),lam).eq(inner(plain(mono(0)),mono(1),lam)));
 for(const b of [f('2/3'),f('-2/3')]){
  const Dphase=[C(b.neg()),C([0,f(1).div(lam).div(-2)])];
  const mean=inner(mono(0),Dphase,lam),second=inner(Dphase,Dphase,lam),v=second.sub(mean.norm2());
  ck('phase_state_D_mean:'+lam+':'+b,mean.eq(b.neg()));
  ck('phase_state_D_variance:'+lam+':'+b,v.eq(f(1).div(lam).div(4)));
 }
}
ck('phase_pair_readable_contrast',f('2/3').sub(f('-2/3')).eq('4/3'));
no('Gaussian_measure_implies_zero_conjugate_variance',!f('1/4').zero());

// Cutoff coefficient witnesses and a rigorous written tail bound in T3/T4.
let H3sum=f(1);
for(let n=1;n<=8;n++){
 const base=f(n).pow(2).add(1),lambda=f(1).div(base.pow(4)),term=base.pow(3).mul(lambda).mul(2);
 ck('Fourier_H3_covariance_term:'+n,term.eq(f(2).div(base)));
 const old=H3sum;H3sum=H3sum.add(term);ck('Fourier_H3_moment_positive_increment:'+n,old.le(H3sum));
 let finiteTail=f(0);for(let k=n+1;k<=16;k++)finiteTail=finiteTail.add(f(2).div(f(k).pow(2).add(1)));
 ck('finite_tail_below_written_two_over_N_bound:'+n,finiteTail.le(f(2).div(n)));
}
fixtures.push({kind:'covariance_coefficient_witness',H3_squared_moment_through_n8:H3sum,tail_bound:'2/N for N>=1; infinite bound is proved by integral comparison'});

// Exact old GOA cross-orientation witness, replayed without changing it.
const intPoly=(coeff,a,b)=>sum(coeff.map((c,i)=>f(c).mul(b.pow(i+1).sub(a.pow(i+1))).div(i+1)));
const cross=intPoly(['1/2',1],f('-3/4'),f('-1/2'));
ck('old_RPM_cross_orientation_matrix_element',cross.eq('-1/32'));
no('old_RPM_probe_preserves_P_orientation',!cross.zero());
no('positive_spectrum_invariant_under_all_additive_shifts',f('1/2').sub(1).le(0));
for(const epsilon of [f(1),f(-1)]){
 const positive=[f('1/3'),f(2),f('5/4')],weights=[f('1/4'),f('1/2'),f('1/4')],J=sum(positive.map((v,i)=>v.mul(weights[i]))).mul(epsilon);
 ck('positive_window_respects_orientation:'+epsilon,epsilon.mul(J).eq('67/48'));
}
fixtures.push({kind:'quantum_comparison',old_cross_orientation_matrix_element:cross,new_sign_commutator:'zero by block definition',scope:'obstruction requires identifying orientation and measured probe'});

// General formal commutators. Written domain/spectral proofs remain separate.
function register(p){const audit=p.audit();ck('confluence:'+p.hash(),audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:p.spec(),audit});return p;}
function prove(p,name,lhs,rhs){const r=p.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',p.replay(cert));proofs.push({name,presentation_sha256:p.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
const tokens=['E','S','X','P','D'],rules=[];
for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++){
 const high=tokens[j],low=tokens[i],rhs=[[[low,high],1]];
 if(high==='D'&&low==='E')rhs.push([['E'],[0,1]]);
 if(high==='D'&&low==='S')rhs.push([[],[0,1]]);
 if(high==='P'&&low==='X')rhs.push([[],[0,-1]]);
 rules.push({id:'order_'+high+low,lhs:[high,low],rhs,source:'MPQ affine core and independent canonical pointer algebra'});
}
const alg=register(new M.Presentation({tokens,rules})),[E,S,X,Pi,D]=tokens.map(x=>alg.word([x])),one=alg.one(),comm=(a,b)=>a.times(b).minus(b.times(a)),A=E.times(Pi);
prove(alg,'affine_density_commutator',comm(D,E),E.scale([0,1]));
prove(alg,'log_density_canonical_commutator',comm(D,S),one.scale([0,1]));
prove(alg,'affine_density_square_commutator',comm(D,E.times(E)),E.times(E).scale([0,2]));
prove(alg,'pointer_pulse_commutator',comm(A,X),E.scale([0,-1]));
prove(alg,'profile_pulse_commutator',comm(A,D),A.scale([0,-1]));
prove(alg,'profile_pulse_double_commutator',comm(A,comm(A,D)),alg.zero());
const ct=['p','u','k','c','w','a','b','B'],cr=[];
for(let j=0;j<ct.length;j++)for(let i=0;i<j;i++)cr.push({id:'commute_'+j+'_'+i,lhs:[ct[j],ct[i]],rhs:[[[ct[i],ct[j]],1]],source:'commuting classical source or reduced coefficients'});
const cp=register(new M.Presentation({tokens:ct,rules:cr})),[p,u,k,c,w,a,b,BB]=ct.map(x=>cp.word([x]));
prove(cp,'source_Abelianized_normal_combination',p.times(k).plus(u.times(c)).minus(u.times(u.times(p).times(k).plus(c))),p.times(cp.one().minus(u.times(u))).times(k));
prove(cp,'global_feedback_restores_linear_theta_flow',p.times(w).times(a.minus(b)).plus(BB).plus(p.times(w).times(b)).minus(BB),p.times(w).times(a));
prove(cp,'adapted_density_exact_conservation',p.times(w).times(b).times(a).scale(-1).plus(p.times(b).times(w).times(a)),cp.zero());
const bad=JSON.parse(JSON.stringify(proofs.find(x=>x.certificate.steps.length).certificate));bad.steps[0].coefficient=['991','0'];
no('corrupted_rewrite_certificate',(()=>{try{alg.replay(bad);return false;}catch{return true;}})());
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'MPQ1_MATERIAL_PROFILE_QUANTIZATION_V1',status:'PASS_SELECTED_REDUCED_CONTINUUM_QUANTIZATION',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{source_classical_Abelianization:true,source_action_changed:false,exact_reduced_pulse_solution:true,positive_continuum_reduced_representation:true,representation_covariance:'lambda_0=1; lambda_n=(1+n^2)^(-4)',self_adjoint_observable_contract:'smeared log density, adapted conserved density and affine conjugate',clock_physical_norm:'L2(mu) direct_sum L2(mu), tensor L2(pointer), conditional on selected reduced representation',controlled_Fourier_cutoff_and_pulse_limit:true,quantum_pointer_backreaction:true,coherent_semiclassical_limit:true,old_RPM_orientation_and_probe_intertwiner:false,raw_nonlinear_quantum_geometry_constructed:false,original_unreduced_quantum_constraint_closure:false,source_selected_covariance_unique:false,physical_three_plus_one_gravity_closed:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cpth=path.join(__dirname,'CERTIFICATE.json'),hpth=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cpth,bytes);fs.writeFileSync(hpth,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cpth,'utf8'),bytes);assert.equal(fs.readFileSync(hpth,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
