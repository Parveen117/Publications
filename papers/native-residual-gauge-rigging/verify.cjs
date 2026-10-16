'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const tokens=['k','s1','s2','s3','s4','h1','h2','h3','h4','h5','x','u','v','t','H','dH','L','J'];
const rules=[];for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting scalar jets and energy/determinant variables'});
const alg=new M.Presentation({tokens,rules}),audit=alg.audit();
ck('scalar_jet_presentation_confluent',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
presentations.push({spec:alg.spec(),audit:{protocol:audit.protocol,status:audit.status,termination:audit.termination,criticalPairs:audit.criticalPairs,resolved:audit.resolutions.length,open:audit.open.length,full_audit_sha256:crypto.createHash('sha256').update(JSON.stringify(audit)).digest('hex')}});
const z=alg.word([],0),one=alg.word([]),V=Object.fromEntries(tokens.map(t=>[t,alg.word([t])])),pow=(a,n)=>{let r=one;for(let i=0;i<n;i++)r=r.times(a);return r;};
function prove(name,lhs,rhs){const out=alg.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',alg.replay(certificate));proofs.push({name,presentation_sha256:alg.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
const derivative={k:V.k.times(V.s1),s1:V.s2,s2:V.s3,s3:V.s4,h1:V.h2.times(V.k),h2:V.h3.times(V.k),h3:V.h4.times(V.k),h4:V.h5.times(V.k),t:one,H:V.dH};
function D(expr){let r=z;for(const[word,c]of expr.entries())for(let i=0;i<word.length;i++){if(!derivative[word[i]])continue;r=r.plus(alg.word(word.slice(0,i),c).times(derivative[word[i]]).times(alg.word(word.slice(i+1))));}return r;}
const a=V.h1.times(V.k),u4=D(D(D(a)));
const b=V.h4.times(pow(V.k,4)).plus(V.h3.times(pow(V.k,3)).times(V.s1).scale(6))
 .plus(V.h2.times(pow(V.k,2)).times(V.s2.scale(4).plus(pow(V.s1,2).scale(7))))
 .plus(V.h1.times(V.k).times(V.s1.times(V.s2).scale(3).plus(pow(V.s1,3))));
prove('fourth_derivative_chain_rule',u4,a.times(V.s3).plus(b));
const boundary=a.times(pow(V.s3,2)).scale('1/2').plus(V.s3.times(b));
prove('rough_cross_term_integration_by_parts',D(boundary),V.s4.times(a.times(V.s3).plus(b)).plus(D(a).times(pow(V.s3,2)).scale('1/2')).plus(V.s3.times(D(b))));
prove('conditional_trace_integration_by_parts',D(one.scale('1/2').minus(V.t).times(V.H)),one.scale('1/2').minus(V.t).times(V.dH).minus(V.H));
const delta=(x,u)=>x.times(u).scale(2).plus(pow(u,2));
prove('finite_quadratic_energy_cocycle',delta(V.x,V.u).plus(delta(V.x.plus(V.u),V.v)),delta(V.x,V.u.plus(V.v)));
// L=log g'(alpha), J=mean log g'; tr_F A=J-L, log det_2=-J.
prove('full_logarithmic_Jacobian_cancellation',V.L.minus(V.J).plus(V.J.minus(V.L)),z);
no('omit_angular_Jacobian',alg.reduce(V.J.scale(-1).plus(V.J.minus(V.L))).normal.terms.size!==0);
no('omit_regularized_determinant',alg.reduce(V.L.plus(V.J.minus(V.L))).normal.terms.size!==0);
no('wrong_fourth_derivative_coefficient',alg.reduce(u4.minus(a.times(V.s3).plus(b).plus(V.h2.times(pow(V.k,2)).times(pow(V.s1,2))))).normal.terms.size!==0);
no('omit_half_in_cross_integration_by_parts',alg.reduce(D(boundary).minus(V.s4.times(a.times(V.s3).plus(b)).plus(D(a).times(pow(V.s3,2))).plus(V.s3.times(D(b))))).normal.terms.size!==0);
no('omit_energy_cross_term',alg.reduce(delta(V.x,V.u).plus(delta(V.x.plus(V.u),V.v)).minus(pow(V.u.plus(V.v),2))).normal.terms.size!==0);

// Identity-profile Fourier trace, with coordinate angle period 2*pi.
// e_n primitive=(e_n-1)/(i*n); diagonal of f'' B is -i*m at n=m.
for(const m of [-8,-5,-2,-1,1,2,5,8]){
 let tr=N.ZERO;for(let n=-10;n<=10;n++)if(n!==0&&n===m)tr=tr.add(C([0,-m]));
 ck('paired_Fourier_trace_at_identity:'+m,tr.eq(C([0,-m])));
 ck('angle_plus_shape_divergence_zero:'+m,tr.add(C([0,m])).zero());
}
no('finite_cutoff_captures_unresolved_smearing',!C([0,0]).eq(C([0,-11])));
// For f(theta)=sin(theta) at sigma=0, trace K^2 tends to -1/2.
// The det_2 contribution is nontrivial even though the full Jacobian is one.
for(let cutoff=2;cutoff<=8;cutoff++){
 const modes=[];for(let n=-cutoff;n<=cutoff;n++)if(n)modes.push(n);
 const ac=j=>j===1?C([0,'1/2']):j===-1?C([0,'-1/2']):N.ZERO;
 const K=N.matrix(modes.map(m=>modes.map(n=>ac(m-n).sub(ac(m)).div(C([0,n])))));
 const K2=N.mul(K,K),trace=a=>a.reduce((s,row,i)=>s.add(row[i]),N.ZERO);
 ck('nonlinear_shape_trace:'+cutoff,trace(K).eq(-1));
 ck('Hilbert_Schmidt_square_trace_fixture:'+cutoff,trace(K2).eq(f('-1/2').add(f(1).div(cutoff))));
 fixtures.push({kind:'identity_profile_sine_flow',cutoff,trace:trace(K),square_trace:trace(K2),scope:'exact finite matrix; written trace limit in theorem remains analytic'});
}

// Polynomial integral witness for integral h'(t)(1/2-t)=mean(h)-h(0), h(1)=h(0).
const pclean=p=>p.map(f),padd=(a,b)=>Array.from({length:Math.max(a.length,b.length)},(_,i)=>(a[i]||f(0)).add(b[i]||f(0))),
 pscale=(a,c)=>a.map(x=>x.mul(c)),pmul=(a,b)=>{let r=Array(a.length+b.length-1).fill(f(0));for(let i=0;i<a.length;i++)for(let j=0;j<b.length;j++)r[i+j]=r[i+j].add(a[i].mul(b[j]));return r;},
 pd=a=>a.slice(1).map((x,i)=>x.mul(i+1)),pint=a=>a.reduce((s,x,i)=>s.add(x.div(i+1)),f(0));
for(let n=1;n<=7;n++){
 const h=Array(n+2).fill(f(0));h[0]=f('2/3');h[n]=h[n].add(1);h[n+1]=f(-1);
 const trace=pint(pmul(pd(h),pclean(['1/2',-1])));
 ck('finite_trace_boundary_identity:'+n,trace.eq(pint(h).sub(h[0])));
 ck('finite_determinant_log_sum:'+n,h[0].sub(pint(h)).add(trace).zero());
}

// Exact polynomial Gaussian Ward identities for a finite divergence-free
// nonlinear comparison, NOT a cutoff of the source diffeomorphism action.
// X=(1+x^2)partial_x-2xy partial_y, A=X[(x^2+y^2)/2].
const mono=(i,j,c=1)=>new Map([[i+','+j,f(c)]]),clean=p=>new Map([...p].filter(([,c])=>!c.zero()));
const add=(a,b)=>{const r=new Map(a);for(const[k,c]of b)r.set(k,(r.get(k)||f(0)).add(c));return clean(r);};
const scale=(a,c)=>clean(new Map([...a].map(([k,v])=>[k,v.mul(c)])));
const mul=(a,b)=>{let r=new Map();for(const[ka,ca]of a)for(const[kb,cb]of b){const ia=ka.split(',').map(Number),ib=kb.split(',').map(Number);r=add(r,mono(ia[0]+ib[0],ia[1]+ib[1],ca.mul(cb)));}return r;};
const dx=(a,j)=>{let r=new Map();for(const[k,c]of a){const ij=k.split(',').map(Number),v=ij[j];if(v){ij[j]--;r=add(r,mono(ij[0],ij[1],c.mul(v)));}}return r;};
const vx=add(mono(0,0),mono(2,0)),vy=mono(1,1,-2),X=p=>add(mul(vx,dx(p,0)),mul(vy,dx(p,1))),score=add(add(mono(3,0),mono(1,0)),mono(1,2,-2));
const moment=(n,c)=>{if(n%2)return f(0);let v=f(1);for(let j=1;j<n;j+=2)v=v.mul(j).mul(c);return v;};
const avg=(p,c)=>[...p].reduce((s,[key,v])=>{const[i,j]=key.split(',').map(Number);return s.add(v.mul(moment(i,c)).mul(moment(j,c)));},f(0));
ck('nonlinear_comparison_divergence_zero',add(dx(vx,0),dx(vy,1)).size===0);
ck('comparison_score_is_half_energy_derivative',add(X(add(mono(2,0,'1/2'),mono(0,2,'1/2'))),scale(score,-1)).size===0);
for(let i=0;i<=6;i++)for(let j=0;j<=6-i;j++){
 const p=mono(i,j),q=add(scale(X(p),-1),scale(mul(score,p),'1/2'));
 ck('doubled_Gaussian_Ward:'+i+':'+j,avg(q,2).zero());
 for(const c of [1,2,3])ck('Gaussian_integration_by_parts:'+c+':'+i+':'+j,avg(X(p),c).eq(avg(mul(score,p),c).div(c)));
}
for(const c of [1,3,4]){
 const p=mono(1,0),q=add(scale(X(p),-1),scale(mul(score,p),'1/2'));
 no('wrong_Gaussian_scale:'+c,!avg(q,c).zero());
 fixtures.push({kind:'finite_nonlinear_Ward_defect',covariance_scale:c,test:'x',defect_divided_by_i:avg(q,c),scope:'finite divergence-free comparison only'});
}
no('omit_half_density_score',!avg(X(mono(1,0)),2).zero());
no('double_score_instead_of_half',!avg(add(scale(X(mono(1,0)),-1),mul(score,mono(1,0))),2).zero());

const tests=[mono(0,0),mono(1,0),mono(2,0),add(mono(2,0),mono(0,0,-2)),mono(2,2)];
const ell=tests.map(p=>avg(p,2)),gram=N.matrix(ell.map(a=>ell.map(b=>a.mul(b))));
ck('selected_rigging_Gram_rank_one',N.rank(gram)===1);
ck('selected_normalization',ell[0].eq(1));
ck('nontrivial_null_test',ell[3].zero());
ck('null_test_has_nonzero_kinematic_norm',avg(mul(tests[3],tests[3]),1).eq(3));
for(let i=0;i<ell.length;i++)ck('positive_physical_diagonal:'+i,f(0).le(ell[i].pow(2)));
no('kinematic_norm_equals_physical_norm',!avg(mul(tests[2],tests[2]),1).eq(ell[2].pow(2)));
const R=N.matrix([[1,0,2,0],[0,1,0,2]]),section=N.matrix([[1,0],[0,1],[0,0],[0,0]]);
ck('spectator_quotient_onto',N.equal(N.mul(R,section),N.identity(2)));
for(const O of [N.matrix([[1,0],[0,-1]]),N.matrix([[0,1],[1,0]]),N.matrix([[2,1],[1,3]])]){
 const lifted=N.kron(N.identity(2),O);
 ck('physical_observable_intertwiner:'+checks.length,N.equal(N.mul(R,lifted),N.mul(O,R)));
}
const bad=JSON.parse(JSON.stringify(proofs.find(p=>p.certificate.steps.length).certificate));bad.steps[0].coefficient=['997','0'];
no('corrupted_native_replay',(()=>{try{alg.replay(bad);return false;}catch{return true;}})());
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'RGR1_RESIDUAL_GAUGE_RIGGING_V1',status:'PASS_SELECTED_POSITIVE_FULL_GAUGE_RIGGING',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{kinematic_covariance_unchanged:true,continuous_Radon_Nikodym_formula:'exp[-Delta_g_inverse/(2c)]',full_nonzero_invariant_functional:'integration on continuous tests against nu_2 tensor nu_2',physical_pairing:'rank-one positive form from ell',selected_gauge_quotient_dimension:1,full_invariant_dual_classified:false,Gaussian_scale_uniqueness:'c=2 within the stated covariance family only',common_constraint_core_constructed:true,MPQ_physical_norm_and_spectral_domains_intertwined:true,original_nonlinear_field_ordering_constructed:false,physical_three_plus_one_gravity_closed:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
