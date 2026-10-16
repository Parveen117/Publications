'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const clean=x=>new Map([...x].filter(([,c])=>!c.zero()));
function add(...xs){const r=new Map();for(const x of xs)for(const[k,c]of x)r.set(k,(r.get(k)||N.ZERO).add(c));return clean(r);}
const scale=(x,c)=>clean(new Map([...x].map(([k,v])=>[k,v.mul(c)])));
const sub=(x,y)=>add(x,scale(y,-1)),same=(x,y)=>sub(x,y).size===0;

// Full homogeneous polynomial jets in (Q,P,Phi,K), computed over exact Cuts.
const poly=rows=>clean(new Map(rows.map(([e,c])=>[e.join(','),C(c)])));
const unit=poly([[[0,0,0,0],1]]),zero=new Map();
const [Q,P,Phi,K]=[0,1,2,3].map(i=>poly([[Array.from({length:4},(_,j)=>i===j?1:0),1]]));
function mul(x,y){const r=new Map();for(const[k,c]of x)for(const[l,d]of y){const a=k.split(',').map(Number),b=l.split(',').map(Number),v=a.map((n,i)=>n+b[i]).join(',');r.set(v,(r.get(v)||N.ZERO).add(c.mul(d)));}return clean(r);}
function diff(x,i){const r=new Map();for(const[k,c]of x){const a=k.split(',').map(Number),n=a[i];if(!n)continue;a[i]--;r.set(a.join(','),c.mul(n));}return clean(r);}
const pb=(x,y)=>add(sub(mul(diff(x,0),diff(y,1)),mul(diff(x,1),diff(y,0))),sub(mul(diff(x,2),diff(y,3)),mul(diff(x,3),diff(y,2))));
function atZero(x,indices){return clean(new Map([...x].filter(([k])=>{const a=k.split(',').map(Number);return indices.every(i=>a[i]===0);})));}
const constraint=mul(P,K),aList=[unit,Q,add(unit,mul(Q,Q)),add(mul(Q,P),scale(Q,2)),add(mul(P,P),mul(Q,mul(Q,P)))];
const bList=[unit,Phi,add(mul(Phi,Phi),K),add(mul(Phi,K),scale(Phi,-2))];
for(const[i,a]of aList.entries())for(const[j,b]of bList.entries()){
 const F=add(scale(unit,3),mul(P,a),mul(K,b));
 ck('weak_Dirac_identity:'+i+':'+j,same(pb(F,constraint),mul(constraint,add(diff(a,0),diff(b,2)))));
 ck('critical_value_constant:'+i+':'+j,same(atZero(F,[1,3]),scale(unit,3)));
 ck('P_branch_scalar_independence:'+i+':'+j,diff(atZero(F,[3]),2).size===0);
 ck('K_branch_geometry_independence:'+i+':'+j,diff(atZero(F,[1]),0).size===0);
 ck('cross_lift_bracket_zero:'+i+':'+j,pb(mul(P,a),mul(K,b)).size===0);
}
for(const[i,a]of aList.entries())for(const[j,u]of aList.entries()){
 const formula=mul(P,add(sub(mul(diff(a,0),u),mul(a,diff(u,0))),mul(P,sub(mul(diff(a,0),diff(u,1)),mul(diff(a,1),diff(u,0))))));
 ck('normal_form_Poisson_closure:'+i+':'+j,same(pb(mul(P,a),mul(P,u)),formula));
}
no('bare_Q_is_global_Dirac',atZero(pb(Q,constraint),[1]).size!==0);
no('bare_Phi_is_global_Dirac',atZero(pb(Phi,constraint),[3]).size!==0);
no('arbitrary_origin_label_is_global_invariant',diff(Q,0).size!==0);

// Periodic mean-zero primitive in finite Fourier coordinates. The common
// spatial factor 2 pi is handled in writing; theta derivative is i*n.
const fourier=rows=>clean(new Map(rows.map(([n,c])=>[n,C(c)])));
const fd=x=>clean(new Map([...x].map(([n,c])=>[n,c.mul(new N.Cut(0,n))])));
const mean=x=>x.get(0)||N.ZERO;
function primitive(x){return clean(new Map([...x].filter(([n])=>n!==0).map(([n,c])=>[n,c.div(new N.Cut(0,n))])));}
for(const meanK of [-3,0,2])for(const n of [1,2,4]){
 const k=fourier([[0,meanK],[n,[2,1]],[-n,[2,-1]]]),r=primitive(k);
 ck('periodic_inverse_derivative:'+meanK+':'+n,same(fd(r),sub(k,fourier([[0,meanK]]))));
 ck('primitive_zero_mean:'+meanK+':'+n,mean(r).zero());
 ck('primitive_uniform_k_shift_invariance:'+meanK+':'+n,same(primitive(add(k,fourier([[0,7]]))),r));
}

// Self-adjoint transport's exact derivative factors (no numerical trig).
for(const row of [['1/2',2,3,5],['3/2',-1,2,4],[1,0,-2,3]]){
 const[a,ap,fp,f0]=row.map(f);
 const transported=a.mul(fp).add(ap.mul(f0).div(2));
 ck('half_density_transport_jet:'+row.join(','),transported.eq(a.mul(fp).add(ap.div(2).mul(f0))));
 if(!ap.zero())no('omit_half_density_correction:'+row.join(','),!transported.eq(a.mul(fp)));
}
for(const eps of ['1/3','1/2','3/4']){
 ck('strictly_positive_flow_speed:'+eps,f(0).le(f(1).sub(eps))&&!f(1).sub(eps).zero());
 ck('bounded_speed_ratio:'+eps,f(1).add(eps).div(f(1).sub(eps)).le(7));
}
no('epsilon_one_has_same_positive_lower_bound',f(1).sub(1).zero());

// One-variable polynomial integration for exact normalized physical packets.
// Step packets are legitimate Dom F vectors; smooth-test regularity is not
// asserted. Interval endpoints and every coefficient are exact rationals.
const cp=rows=>rows.map(C);
function cadd(x,y){return Array.from({length:Math.max(x.length,y.length)},(_,i)=>(x[i]||N.ZERO).add(y[i]||N.ZERO));}
const cscale=(x,s)=>x.map(c=>c.mul(s));
function cmul(x,y){const r=Array(x.length+y.length-1).fill(N.ZERO);for(let i=0;i<x.length;i++)for(let j=0;j<y.length;j++)r[i+j]=r[i+j].add(x[i].mul(y[j]));return r;}
function integrate(x,lo,hi){lo=f(lo);hi=f(hi);return x.reduce((s,c,i)=>s.add(c.mul(hi.pow(i+1).sub(lo.pow(i+1)).div(i+1))),N.ZERO);}
const abs2=x=>cmul(x.map(c=>c.dagger()),x);
const loP='1/4',hiP='1/2',loM='-3/4',hiM='-1/2',eps=f('1/2'),packet=cp([2]);
const negativeAction=cp(['1/4','1/2']),positiveShiftAction=cp(['-1/4','1/2']);
const normP=integrate(abs2(packet),loP,hiP),normM=integrate(abs2(packet),loM,hiM);
ck('positive_packet_unit_norm',normP.eq(1));ck('negative_packet_unit_norm',normM.eq(1));
const transition=integrate(abs2(negativeAction),loM,hiM);
ck('orientation_transition_coefficient',transition.eq('1/768'));
const cross=integrate(cmul(packet,negativeAction),loM,hiM);
ck('off_diagonal_observable_element',cross.eq('-1/32'));
const diagP=integrate(cmul(abs2(packet),cp([0,1])),loP,hiP),diagM=integrate(cmul(abs2(packet),cp([0,1])),loM,hiM);
ck('positive_packet_diagonal',diagP.eq('3/8'));ck('negative_packet_diagonal',diagM.eq('-5/8'));
const mainNorm=integrate(abs2(cp([0,2])),loP,hiP),otherNorm=integrate(abs2(positiveShiftAction),'5/4','3/2');
const F2=mainNorm.add(otherNorm).add(transition);
ck('complete_F_packet_norm_squared',F2.eq('25/128'));
ck('positive_packet_F_variance',F2.sub(diagP.mul(diagP)).eq('7/128'));
const alpha=f('3/5'),beta=f('4/5'),base=diagP.mul(alpha.pow(2)).add(diagM.mul(beta.pow(2))),interference=cross.mul(alpha.mul(beta).mul(2));
ck('coherent_packets_unit_norm',alpha.pow(2).add(beta.pow(2)).eq(1));
ck('phase_plus_expectation',base.add(interference).eq('-59/200'));
ck('phase_minus_expectation',base.sub(interference).eq('-47/200'));
ck('phase_contrast',interference.mul(2).eq('-3/50'));
ck('same_positive_probability',alpha.pow(2).eq('9/25'));ck('same_negative_probability',beta.pow(2).eq('16/25'));
no('affine_algebra_reads_relative_sign',!cross.zero());
no('orientation_projection_reduces_F',!transition.zero());
no('phase_blind_mixture_has_same_F_reading',!base.eq(base.add(interference)));
// Reverse matrix element: positive support sees psi_-(P-1).
const reverseAction=cp(['-1/4','1/2']);
ck('symmetric_shift_coefficients',integrate(cmul(packet,reverseAction),loP,hiP).eq(cross));
const wrongLeft=cp(['1/2','1/2']); // epsilon*(P+1) on the negative output, unsymmetrized a(Q)P
const wrongCross=integrate(cmul(packet,wrongLeft),loM,hiM);
no('unsymmetric_ordering_preserves_cross_element',!wrongCross.eq(cross));
fixtures.push({kind:'physical_packets',epsilon:eps,normP,normM,transition,cross,diagP,diagM,F2,variance:F2.sub(diagP.mul(diagP)),phase_plus:base.add(interference),phase_minus:base.sub(interference),contrast:interference.mul(2)});

// Finite block witnesses, separately scoped from the infinite commutant proof.
const E=(i,j)=>{const m=N.zeros(5);m[i][j]=N.ONE;return m;};
const projections=Array.from({length:5},(_,i)=>E(i,i)),XP=N.add(E(0,1),E(1,0)),XK=N.add(E(2,3),E(3,2));
const minimal=N.operatorAlgebra(projections),extended=N.operatorAlgebra([...projections,XP,XK]);
ck('finite_affine_center_witness_dimension',minimal.dimension===5);
ck('finite_extended_block_witness_dimension',extended.dimension===9);
ck('finite_center_equalities_rank',N.rank([[1,-1,0,0,0],[0,0,1,-1,0]])===2);
const branchP=N.add(projections[0],projections[1]),branchK=N.add(projections[2],projections[3]);
for(const[g,x]of [['XP',XP],['XK',XK]]){
 ck('P_branch_projection_commutes:'+g,N.isZero(N.commutator(branchP,x)));
 ck('K_branch_projection_commutes:'+g,N.isZero(N.commutator(branchK,x)));
 ck('origin_projection_commutes:'+g,N.isZero(N.commutator(projections[4],x)));
}
no('positive_P_projection_remains_central',!N.isZero(N.commutator(projections[0],XP)));
no('positive_K_projection_remains_central',!N.isZero(N.commutator(projections[2],XK)));
no('arbitrary_cross_branch_matrix_is_in_constructed_algebra',!N.isZero(N.commutator(branchP,N.add(E(0,2),E(2,0)))));
no('arbitrary_origin_matrix_is_in_constructed_algebra',!N.isZero(N.commutator(projections[4],N.add(E(0,4),E(4,0)))));

// Native word replays for normal-form multiplication and the Poisson factors.
const tokens=['a','b','c','d','p','k'],rules=[];
for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting classical observable coefficients'});
const scalar=new M.Presentation({tokens,rules}),audit=scalar.audit();ck('word_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:scalar.spec(),audit});
function prove(name,lhs,rhs){const r=scalar.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',scalar.replay(cert));proofs.push({name,presentation_sha256:scalar.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
const[a,b,c,d,p,k]=tokens.map(x=>scalar.word([x]));
prove('normal_form_product_mod_constraint',c.plus(p.times(a)).plus(k.times(b)).times(d.plus(p.times(b)).plus(k.times(a))),c.times(d).plus(p.times(c.times(b).plus(d.times(a)).plus(p.times(a).times(b)))).plus(k.times(c.times(a).plus(d.times(b)).plus(k.times(a).times(b)))).plus(p.times(k).times(a.times(a).plus(b.times(b)))));
prove('source_weak_constraint_factor',k.times(p.times(a)).plus(p.times(k.times(b))),p.times(k).times(a.plus(b)));
prove('branch_bracket_factor',p.times(a).times(b.plus(p.times(c))).minus(a.plus(p.times(d)).times(p.times(b))),p.times(p).times(a.times(c).minus(d.times(b))));

const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'GOA1_GLOBAL_OBSERVABLE_ALGEBRA_V1',status:'PASS_SMOOTH_OBSERVABLE_ALGEBRA_AND_ORIENTATION_MIXING',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{smooth_classical_value_algebra:'c+P a(Q,P) on K=0; c+K b(Phi,K) on P=0; actual periodic-primitive local lift',minimal_affine_quantum_center:'C^5 for four momentum orientations and the origin',chosen_nonlinear_observable:'symmetric P[1+(1/2)cos Q], self-adjoint on Dom P after physical weight transport',orientation_transition_coefficient:'1/768 at small observable-unitary parameter',phase_expectations:['-59/200','-47/200'],phase_contrast:'-3/50',extended_quantum_algebra:'B(H_P) direct_sum B(H_K) direct_sum C_origin',extended_center:'C^3',P_to_K_or_origin_transition_constructed:false,general_quantization_of_all_smooth_observables:false,physical_time_dynamics_constructed:false,calibrated_gravitational_prediction:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cpth=path.join(__dirname,'CERTIFICATE.json'),hpth=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cpth,bytes);fs.writeFileSync(hpth,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cpth,'utf8'),bytes);assert.equal(fs.readFileSync(hpth,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
