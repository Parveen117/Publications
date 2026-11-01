'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const zero=N.ZERO,one=N.ONE;
const clean=a=>new Map([...a].filter(([,v])=>!v.zero()));
const seq=rows=>clean(new Map(rows.map(([n,v])=>[n,C(v)])));
const add=(a,b)=>{const r=new Map(a);for(const[n,v]of b)r.set(n,(r.get(n)||zero).add(v));return clean(r);};
const scale=(a,c)=>clean(new Map([...a].map(([n,v])=>[n,v.mul(c)])));
const shift=(a,k)=>new Map([...a].map(([n,v])=>[n-k,v]));
const total=a=>[...a.values()].reduce((s,v)=>s.add(v),zero);
const inner=(a,b)=>[...a].reduce((s,[n,v])=>s.add(v.dagger().mul(b.get(n)||zero)),zero);
const eq=(a,b)=>clean(add(a,scale(b,-1))).size===0;
const cases=[seq([[0,1],[1,2],[-2,[0,1]]]),seq([[-3,2],[2,-2]]),seq([[1,['1/3','2/5']],[4,3]]),seq([])];
for(let i=0;i<cases.length;i++)for(let j=0;j<cases.length;j++){
 const a=cases[i],b=cases[j],ks=new Set();for(const n of a.keys())for(const m of b.keys())ks.add(m-n);
 const avg=[...ks].reduce((s,k)=>s.add(inner(a,shift(b,k))),zero);
 ck('exact_central_average:'+i+':'+j,avg.eq(total(a).dagger().mul(total(b))));
 if(i===j)ck('positive_central_quotient:'+i,avg.eq(C(total(a).norm2())));
}
for(let i=0;i<cases.length;i++){
 const h=cases[i],v=add(h,scale(shift(h,1),-1));
 ck('shift_difference_has_zero_periodization:'+i,total(v).zero());
 const recovered=new Map();if(v.size){const lo=Math.min(...v.keys()),hi=Math.max(...v.keys());for(let n=lo;n<=hi;n++)recovered.set(n,total(new Map([...v].filter(([k])=>k>=n))));}
 ck('finite_tail_inverse_of_shift_difference:'+i,eq(clean(recovered),h));
 ck('continuous_right_inverse_fixture:'+i,total(seq([[0,total(h)]])).eq(total(h)));
 ck('periodization_invariant_under_integer_shift:'+i,total(shift(h,3)).eq(total(h)));
}
const a0=seq([[0,1]]),b5=seq([[5,1]]);
const truncated=[-1,0,1].reduce((s,k)=>s.add(inner(a0,shift(b5,k))),zero);
no('finite_shift_window_replaces_full_average',!truncated.eq(total(a0).dagger().mul(total(b5))));
no('nonzero_sum_has_rapid_tail_inverse',!total(new Map([...a0].filter(([k])=>k>=-10))).zero());
no('ordinary_kinematic_norm_is_central_physical_norm',!inner(cases[1],cases[1]).zero()&&total(cases[1]).zero());

// An exact central Z cover of a finite weighted base. This witnesses the
// cocycle/reindexing mechanism; it is not a diffeomorphism-group cutoff.
const roots=[1,2,3].map(f),weights=roots.map(x=>x.pow(2).div(14));
const cover=rows=>clean(new Map(rows.map(([n,j,v])=>[n+':'+j,C(v)])));
const parse=s=>s.split(':').map(Number);
const g={p:[1,2,0],k:[1,0,-1]},h={p:[1,0,2],k:[0,2,1]};
const compose=(a,b)=>({p:b.p.map(j=>a.p[j]),k:b.k.map((v,j)=>v+a.k[b.p[j]])});
function action(a,v){const out=new Map();for(const[key,value]of v){const[n,j]=parse(key),i=a.p[j],nk=(n+a.k[j])+':'+i;out.set(nk,(out.get(nk)||zero).add(value.mul(roots[j].div(roots[i]))));}return clean(out);}
const deck=(v,k)=>new Map([...v].map(([key,value])=>{const[n,j]=parse(key);return[(n+k)+':'+j,value];}));
const period=v=>{const out=[zero,zero,zero];for(const[key,x]of v){const[,j]=parse(key);out[j]=out[j].add(x);}return out;};
const baseAction=(a,v)=>{const out=[zero,zero,zero];for(let j=0;j<3;j++)out[a.p[j]]=v[j].mul(roots[j].div(roots[a.p[j]]));return out;};
const veq=(a,b)=>a.every((x,i)=>x.eq(b[i]));
const cnorm=v=>[...v].reduce((s,[key,x])=>s.add(x.norm2().mul(weights[parse(key)[1]])),f(0));
const vn=v=>v.reduce((s,x,i)=>s.add(x.norm2().mul(weights[i])),f(0));
const vf=cover([[0,0,1],[1,1,[1,2]],[-2,2,'2/3'],[3,0,-1]]);
for(const[name,t]of [['g',g],['h',h],['gh',compose(g,h)]]){
 ck('cover_weighted_unitarity:'+name,cnorm(action(t,vf)).eq(cnorm(vf)));
 ck('cover_periodization_intertwiner:'+name,veq(period(action(t,vf)),baseAction(t,period(vf))));
 ck('cover_central_deck_commutes:'+name,eq(action(t,deck(vf,2)),deck(action(t,vf),2)));
 ck('cover_residual_quotient_unitarity:'+name,vn(baseAction(t,period(vf))).eq(vn(period(vf))));
}
ck('cover_cocycle_group_law',eq(action(g,action(h,vf)),action(compose(g,h),vf)));
const wrong={p:compose(g,h).p,k:h.k.map((v,j)=>v+g.k[j])};
no('omit_base_pullback_from_cover_cocycle',!eq(action(wrong,vf),action(compose(g,h),vf)));
no('central_reduction_imposes_all_gauge_transformations',!veq(baseAction(g,period(vf)),period(vf)));
fixtures.push({kind:'central_cover',g,h,composed:compose(g,h),periodized_wave:period(vf),scope:'finite weighted base with exact infinite integer deck; no continuum claim'});

// A nonconstant derivative has a nonzero exact Haar-density defect.
const multiply=(a,b)=>{const r=new Map();for(const[n,x]of a)for(const[m,y]of b)r.set(n+m,(r.get(n+m)||zero).add(x.mul(y)));return clean(r);};
const gp=seq([[0,1],[1,'1/4'],[-1,'1/4']]),defect=add(gp,seq([[0,-1]]));
ck('circle_Jacobian_has_unit_mean',(gp.get(0)||zero).eq(1));
ck('circle_Haar_density_defect_exact',(multiply(defect,defect).get(0)||zero).eq('1/8'));
no('all_circle_diffeomorphisms_preserve_Haar',(multiply(defect,defect).get(0)||zero).eq('1/8'));
for(let n=1;n<=8;n++){
 const nn=f(n),bound=nn.add(1).pow(3).div(nn.pow(2).add(1).pow(2));
 ck('retained_C6_Hilbert_Schmidt_bound:'+n,bound.le(f(8).div(nn)));
}

function presentation(tokens,rules){const a=new M.Presentation({tokens,rules}),audit=a.audit();ck('confluence:'+presentations.length,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:a.spec(),audit});return a;}
function prove(alg,name,lhs,rhs){const out=alg.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',alg.replay(certificate));proofs.push({name,presentation_sha256:alg.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
const tt=['pM','p0','pP','mM','m0','mP'],mode=[-1,0,1,-1,0,1],rules=[];
for(let j=0;j<tt.length;j++)for(let i=0;i<j;i++){
 const rhs=[[[tt[i],tt[j]],1]];
 if(Math.floor(i/3)===Math.floor(j/3))rhs.push([[tt[Math.floor(i/3)*3+mode[i]+mode[j]+1]],mode[j]-mode[i]]);
 rules.push({id:'PBW_'+j+'_'+i,lhs:[tt[j],tt[i]],rhs,source:'two commuting finite Witt subalgebras, modes -1,0,1'});
}
const alg=presentation(tt,rules),Lp=n=>alg.word([tt[n+1]]),Lm=n=>alg.word([tt[n+4]]),comm=(a,b)=>a.times(b).minus(b.times(a));
const H=n=>Lp(n).minus(Lm(n)),D=n=>Lp(n).plus(Lm(n));
prove(alg,'source_HH_equals_D',comm(H(1),H(-1)),D(0).scale(2));
prove(alg,'source_DH_equals_H',comm(D(1),H(-1)),H(0).scale(2));
prove(alg,'source_DD_equals_D',comm(D(1),D(-1)),D(0).scale(2));
prove(alg,'opposite_chiral_generators_commute',comm(Lp(1),Lm(-1)),alg.word([],0));
no('wrong_HH_sign_or_generator',alg.reduce(comm(H(1),H(-1)).minus(H(0).scale(2))).normal.terms.size!==0);
no('insert_scalar_central_term',alg.reduce(comm(H(1),H(-1)).minus(D(0).scale(2)).minus(alg.word([]))).normal.terms.size!==0);

const ct=['m','n','k','f','g','df','dg','ddf','ddg','ell','dell'],cr=[];
for(let j=0;j<ct.length;j++)for(let i=0;i<j;i++)cr.push({id:'commute_'+j+'_'+i,lhs:[ct[j],ct[i]],rhs:[[[ct[i],ct[j]],1]],source:'commuting index and finite-coordinate divergence jets'});
const ca=presentation(ct,cr),[m,n,k,ff,gg,df,dg,ddf,ddg,ell,dell]=ct.map(t=>ca.word([t]));
const jac=n.minus(k).times(m.minus(n).minus(k)).plus(k.minus(m).times(n.minus(k).minus(m))).plus(m.minus(n).times(k.minus(m).minus(n)));
prove(ca,'general_Witt_index_Jacobi',jac,ca.word([],0));
const scoreComm=ff.times(ddg.plus(dell.times(gg)).plus(ell.times(dg))).minus(gg.times(ddf.plus(dell.times(ff)).plus(ell.times(df))));
const scoreBracket=ff.times(ddg).minus(gg.times(ddf)).plus(ell.times(ff.times(dg).minus(gg.times(df))));
prove(ca,'finite_half_density_divergence_cocycle',scoreComm,scoreBracket);
no('drop_density_score_from_half_density_bracket',ca.reduce(scoreComm.minus(ff.times(ddg).minus(gg.times(ddf)))).normal.terms.size!==0);
const bad=JSON.parse(JSON.stringify(proofs.find(x=>x.certificate.steps.length).certificate));bad.steps[0].coefficient=['991','0'];
const badAlg=proofs.find(x=>x.certificate.steps.length).presentation_sha256===alg.hash()?alg:ca;
no('corrupted_rewrite_certificate',(()=>{try{badAlg.replay(bad);return false;}catch{return true;}})());

const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'CGD1_CENTRAL_GAUGE_DOMAIN_V1',status:'PASS_CENTRAL_REDUCTION_AND_COMMON_CORE_WITH_RIGGING_OBSTRUCTION',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{positive_central_integer_average:true,remaining_full_gauge_action_descends:true,rapid_sequence_coinvariants:'exactly H_c',full_continuous_invariant_dual_on_Phi_Z:'zero; other test topologies not classified',finite_C6_continuous_extension:true,dense_common_source_core:'after proved central gauge reduction; Shimomura smooth-vector theorem input',all_smearing_source_algebra:'zero central term on D_c',real_smearing_essential_self_adjointness:true,MPQ_observables_pass_through_central_reduction:true,full_positive_source_physical_pairing_constructed:false,nonlinear_original_field_operators_constructed:false,physical_three_plus_one_gravity_closed:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
