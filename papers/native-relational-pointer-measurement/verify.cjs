'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};

// Exact piecewise-polynomial action of the inherited *full* shifted operator.
// Interval j means [1/4+j,1/2+j]; all coefficients use the native Cut field.
const poly=xs=>xs.map(C);
const add=(a,b)=>Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]||N.ZERO).add(b[i]||N.ZERO));
const scale=(a,c)=>a.map(x=>x.mul(c));
function mul(a,b){const r=Array(a.length+b.length-1).fill(N.ZERO);a.forEach((x,i)=>b.forEach((y,j)=>r[i+j]=r[i+j].add(x.mul(y))));return r;}
function shift(a,s){let r=[],p=[N.ONE];for(const c of a){r=add(r,scale(p,c));p=mul(p,poly([s,1]));}return r;}
function integrate(a,lo,hi){lo=f(lo);hi=f(hi);return a.reduce((r,c,i)=>r.add(c.mul(hi.pow(i+1).sub(lo.pow(i+1)).div(i+1))),N.ZERO);}
function action(a){const r=new Map(),put=(j,b)=>r.set(j,add(r.get(j)||[],b));for(const[j,c]of a){put(j,mul(c,poly([0,1])));put(j+1,scale(mul(shift(c,-1),poly([-1,2])),'1/8'));put(j-1,scale(mul(shift(c,1),poly([1,2])),'1/8'));}return r;}
function inner(a,b){let r=N.ZERO;for(const[j,c]of a){const d=b.get(j);if(d)r=r.add(integrate(mul(c.map(x=>x.dagger()),d),f('1/4').add(j),f('1/2').add(j)));}return r;}
const pp=new Map([[0,poly([2])]]),pm=new Map([[-1,poly([2])]]),fp=action(pp),fm=action(pm);
ck('P_packet_normalized',inner(pp,pp).eq(1));ck('M_packet_normalized',inner(pm,pm).eq(1));ck('packets_orthogonal',inner(pp,pm).zero());
ck('GOA_plus_mean',inner(pp,fp).eq('3/8'));ck('GOA_minus_mean',inner(pm,fm).eq('-5/8'));
ck('GOA_cross_F',inner(pm,fp).eq('-1/32'));ck('F_cross_adjoint',inner(pm,fp).eq(inner(pp,fm).dagger()));
ck('full_F_plus_second_moment',inner(fp,fp).eq('25/128'));ck('full_F_minus_second_moment',inner(fm,fm).eq('61/128'));
ck('full_F_cross_second_moment',inner(fp,fm).eq('1/96'));
ck('central_cross_integral_P',integrate(poly([0,'-1/2',1]),'1/4','1/2').eq('-1/96'));
ck('central_cross_integral_M',integrate(poly([0,'1/2',1]),'-3/4','-1/2').eq('1/48'));
const moments=[];
for(const[s,mean,m2,varF]of [[1,'-59/200','1233/3200','23863/80000'],[-1,'-47/200','1169/3200','24807/80000']]){
 const psi=new Map([[0,poly(['6/5'])],[-1,poly([C('8/5').mul(s)])]]),Fpsi=action(psi),mu=inner(psi,Fpsi),second=inner(Fpsi,Fpsi),variance=second.sub(mu.mul(mu));
 ck('phase_packet_norm:'+s,inner(psi,psi).eq(1));ck('phase_mean:'+s,mu.eq(mean));ck('phase_second_moment:'+s,second.eq(m2));ck('phase_variance:'+s,variance.eq(varF));
 ck('second_moment_via_F_squared:'+s,inner(psi,action(Fpsi)).eq(second));
 // All these compact-support packets happen also to be in Dom F^2. The
 // continuum moment theorem only requires Dom F, as stated in the paper.
 no('omit_second_moment_interference:'+s,!second.eq(C('25/128').mul('9/25').add(C('61/128').mul('16/25'))));
 moments.push({sign:s,mean:mu.rad,second_moment:second.rad,variance:variance.rad,action:[...Fpsi]});
}
const contrast=moments[0].mean.sub(moments[1].mean);
ck('inherited_phase_contrast',contrast.eq('-3/50'));
const truncPlus=inner(pp,fp).norm2().add(inner(pm,fp).norm2());
no('two_packet_compression_preserves_F_second_moment',!truncPlus.eq(inner(fp,fp).rad));
no('phase_mixture_equals_either_pure_readout',!moments[0].mean.eq(moments[1].mean));

// Rational coefficients of Gaussian translation moments and Chebyshev bounds.
// Gaussian normalization/integration and the continuum theorem are written,
// not replaced with finite quadrature or a rational mock exponential.
for(const kappa of ['-2','1/2',1,3])for(const sigma2 of ['1/4',1,4]){
 const k=f(kappa),s2=f(sigma2),delta=k.mul(contrast),vmax=moments[1].variance;
 ck('nonzero_pointer_contrast:'+kappa+':'+sigma2,!delta.zero());
 const distance2=delta.pow(2).div(4),bound=s2.add(k.pow(2).mul(vmax)).div(distance2);
 ck('Chebyshev_coefficient:'+kappa+':'+sigma2,bound.eq(f(10000).div(9).mul(s2.div(k.pow(2)).add(vmax))));
}
const error=f(104807).div(72).div(30000);
ck('unit_pointer_worst_variance',f(1).add(moments[1].variance).eq('104807/80000'));
ck('unit_pointer_error_bound',error.eq('104807/2160000'));
ck('thirty_thousand_run_bound_below_one_twentieth',error.le('1/20'));
ck('zero_coupling_mean_contrast',contrast.mul(0).zero());
no('zero_coupling_distinguishes_phase',contrast.mul(0).zero());
no('opposite_pointer_shift_sign_is_same_prediction',!contrast.eq(contrast.neg()));

// Complete the square behind the *analytic* Gaussian overlap formula.
for(const l of [-2,0,3])for(const m of [-1,2])for(const k of ['1/2',2]){
 const a=f(k).mul(l),b=f(k).mul(m),mid=a.add(b).div(2);
 const lhs=add(mul(poly([a.neg(),1]),poly([a.neg(),1])),mul(poly([b.neg(),1]),poly([b.neg(),1])));
 const rhs=add(scale(mul(poly([mid.neg(),1]),poly([mid.neg(),1])),2),poly([a.sub(b).pow(2).div(2)]));
 ck('Gaussian_overlap_completed_square:'+l+':'+m+':'+k,lhs.every((x,i)=>x.eq(rhs[i])));
}
no('wrong_Gaussian_dephasing_denominator',!f(1).div(8).eq(f(1).div(4)));

// Finite instrument analogue, explicitly not a discretization of F's spectrum.
const K0=N.matrix([['3/5',0],[0,'4/5']]),K1=N.matrix([['4/5',0],[0,'3/5']]),I=N.identity(2);
const V=N.matrix([...K0,...K1]),Vd=N.dagger(V),R=N.mul(V,Vd);
ck('finite_Kraus_normalization',N.equal(N.add(N.mul(N.dagger(K0),K0),N.mul(N.dagger(K1),K1)),I));
ck('finite_preparation_isometry',N.equal(N.mul(Vd,V),I));ck('finite_image_projection',N.equal(N.mul(R,R),R));
const J=O=>N.mul(N.mul(V,O),Vd),e01=N.matrix([[0,1],[0,0]]),e10=N.dagger(e01);
ck('finite_image_product',N.equal(J(N.mul(e01,e10)),N.mul(J(e01),J(e10))));
ck('finite_image_adjoint',N.equal(J(e10),N.dagger(J(e01))));
ck('finite_inverse_on_coherent_image',N.equal(N.mul(N.mul(Vd,J(e01)),V),e01));
ck('finite_image_identity_is_R',N.equal(J(I),R));
no('image_identity_is_whole_record_identity',!N.equal(R,N.identity(4)));
const channel=O=>N.add(N.mul(N.mul(K0,O),N.dagger(K0)),N.mul(N.mul(K1,O),N.dagger(K1)));
ck('finite_channel_off_diagonal_attenuation',N.equal(channel(e01),N.scale(e01,'24/25')));
ck('finite_channel_diagonal_preservation',N.equal(channel(N.mul(e01,e10)),N.mul(e01,e10)));
no('coarse_effect_identity_forces_identity_instrument',!N.equal(channel(e01),e01));
no('individual_Kraus_operator_is_projector',!N.equal(N.mul(K0,K0),K0));

// Native Weyl-word witnesses. Analytic self-adjointness is *not* established
// by these formal identities; the spectral/conjugation proofs supply it.
const weyl=new M.Presentation({tokens:['F','X','p'],rules:[
 {id:'XF',lhs:['X','F'],rhs:[[['F','X'],1]],source:'separate system and pointer tensor factors'},
 {id:'pF',lhs:['p','F'],rhs:[[['F','p'],1]],source:'separate system and pointer tensor factors'},
 {id:'pX',lhs:['p','X'],rhs:[[['X','p'],1],[[],[0,-1]]],source:'Pi=-i partial_X on Schwartz pointer core'}
]});
function register(p){const audit=p.audit();ck('confluence:'+presentations.length,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:p.spec(),audit});}
function prove(p,name,lhs,rhs){const r=p.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',p.replay(cert));proofs.push({name,presentation_sha256:p.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
register(weyl);
const [F,X,p]=['F','X','p'].map(x=>weyl.word([x])),A=F.times(p),comm=(a,b)=>a.times(b).minus(b.times(a));
prove(weyl,'pointer_Heisenberg_derivative',comm(A,X).scale(N.IOTA),F);
prove(weyl,'second_pointer_commutator_zero',comm(A,comm(A,X)),weyl.zero());
prove(weyl,'measured_observable_commutator_zero',comm(A,F),weyl.zero());
const wrong=weyl.reduce(comm(A,X).scale(N.IOTA).plus(F));
no('wrong_pointer_commutator_sign',wrong.normal.terms.size!==0);
const iso=new M.Presentation({tokens:['a','b','W','Wd'],rules:[{id:'isometry',lhs:['Wd','W'],rhs:[[[],1]],source:'RPM-T3/T5 physical record isometry'}]});
register(iso);const[a,b,W,Wd]=['a','b','W','Wd'].map(x=>iso.word([x])),enc=o=>W.times(o).times(Wd);
prove(iso,'coherent_corner_multiplication',enc(a).times(enc(b)),enc(a.times(b)));
prove(iso,'coherent_corner_recovery',Wd.times(enc(a)).times(W),a);

fixtures.push({kind:'GOA_full_shifted_packet_moments',moments,contrast,finite_run_error_bound:error,finite_run_count:30000});
fixtures.push({kind:'finite_instrument_analogue_not_spectral_truncation',K0,K1,V,off_diagonal_multiplier:'24/25'});
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'RPM1_RELATIONAL_POINTER_MEASUREMENT_V1',status:'PASS_IDEAL_RELATIONAL_POINTER_EXTENSION',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{clock_constraint:'W Pi_T W^* with domain W H^1',rigging_form:'<integral W^*f, integral W^*h>',physical_completion:'inherited gravitational H tensor L^2(pointer)',pointer_instrument:'K_x=phi_sigma(x-kappa F)',phase_pointer_contrast:'-3 kappa/50',observable_variances:['23863/80000','24807/80000'],added_ideal_clock_and_global_coupling:true,local_LNC_apparatus_derivation:false,source_selected_gravity_dynamics:false,physical_calibration:false,inter_branch_coupling:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cpth=path.join(__dirname,'CERTIFICATE.json'),hpth=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cpth,bytes);fs.writeFileSync(hpth,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cpth,'utf8'),bytes);assert.equal(fs.readFileSync(hpth,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
