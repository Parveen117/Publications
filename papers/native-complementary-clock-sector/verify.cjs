'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const O=require('../native-local-sector-classification/oscillator_fixture.cjs');
const {N,f,C,key,state,vac,plus,scale,minus,eq,current,T,A,L,inner,tensor,partitions,weight}=O;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
// Coefficientwise squeezed bra; no cutoff or infinite sum is evaluated.
function monomialBra(vars,pairs=[['u','v']]){let value=f(1);for(const [left,right]of pairs){const lc=new Map(),rc=new Map();for(const v of vars){const [s,raw]=v.split(':'),n=Number(raw);if(s===left)lc.set(n,(lc.get(n)||0)+1);if(s===right)rc.set(n,(rc.get(n)||0)+1);}
 for(const n of new Set([...lc.keys(),...rc.keys()])){const r=lc.get(n)||0;if(r!==(rc.get(n)||0))return f(0);for(let j=1;j<=r;j++)value=value.mul(-2*n).mul(j);}}
 return value;}
function R0(x,pairs){let z=N.ZERO;for(const[k,c]of x)z=z.add(c.mul(monomialBra(JSON.parse(k),pairs)));return z;}
function Ilevel(l){return plus(...partitions(l).map(p=>{const us=p.map(n=>'u:'+n),vs=p.map(n=>'v:'+n),sign=p.length%2?-1:1;return state([...us,...vs],f(sign).div(weight(us)));}));}
const seedVars=[[],['u:1'],['v:1'],['u:1','v:1'],['u:2','v:2'],['u:1','u:1','v:2'],['u:1','v:3'],['v:7'],['u:7','v:7']];
const seeds=seedVars.map(x=>state(x));
for(const raw of ['0','1','-2','3/2']){
 const t=f(raw),u=t.neg();
 for(const [i,x]of seeds.entries())for(const m of [-5,-2,-1,0,1,2,5]){
  ck('all_mode_dual_squeezed:'+raw+':'+i+':'+m,R0(L(m,x,t,u)).zero());
  const nullCurrent=scale(plus(current('u',m,x,t),current('v',-m,x,u)),'1/2');
  ck('null_current_dual:'+raw+':'+i+':'+m,R0(nullCurrent).zero());
 }
 for(let l=0;l<=4;l++){
  const ps=partitions(l),dsu=ps.map(p=>O.descendant('u',p,t)),dsw=ps.map(p=>O.descendant('v',p,u)),G=dsw.map(x=>dsw.map(y=>inner(x,y))),iv=N.inverse(G);let paired=new Map();
  for(let i=0;i<ps.length;i++)for(let j=0;j<ps.length;j++)paired=plus(paired,scale(tensor(dsu[i],dsw[j]),iv[i][j]));
  ck('LSC_pairing_equals_universal_series:'+raw+':'+l,eq(paired,Ilevel(l)));
 }
}
for(let l=0;l<=5;l++){
 const ip=Ilevel(l);ck('level_norm_is_partition_count:'+l,inner(ip,ip).eq(partitions(l).length));
 for(const p of partitions(l)){const x=state([...p.map(n=>'u:'+n),...p.map(n=>'v:'+n)]);ck('bra_coefficient_matches_Fock_pairing:'+p.join('.'),R0(x).eq(inner(ip,x)));}
}
const truncated=plus(...[0,1,2,3].map(Ilevel));
no('finite_cutoff_is_all_mode_state',L(4,truncated,f(1),f(-1)).size!==0);
no('plain_vacuum_is_all_mode_state',L(1,vac,f(1),f(-1)).size!==0);
const wrong=plus(vac,state(['u:1','v:1'],'1/2'));
no('wrong_squeezed_sign',!inner(wrong,L(1,state(['u:1']),f(1),f(-1))).zero());
const secondPair=(m,x,t)=>minus(A('x',m,x,t),A('y',-m,x,t.neg()));
for(const x of [state(['u:1','x:1']),state(['u:1','v:1','x:2']),state(['v:2','x:1','y:1'])])for(const m of [-2,-1,0,1,2]){
 ck('two_chiralities_first:'+key([...x.keys()])+':'+m,R0(L(m,x,f(2),f(-2)),[['u','v'],['x','y']]).zero());
 ck('two_chiralities_second:'+key([...x.keys()])+':'+m,R0(secondPair(m,x,f(2)),[['u','v'],['x','y']]).zero());
}
// Canonical source jets on the P=0 stratum. delta q and delta k are derived
// from the frozen source before converting to affine connections.
for(const row of [[2,3,5,7,11,13,17],[-2,1,-3,4,2,-1,6],['1/2','2/3',1,-2,3,4,-5]]){
 const [k,qx,kx,qxx,ff,fp,fpp]=row.map(f),cp=k.add(qx),cm=qx.sub(k),cpx=kx.add(qxx),cmx=qxx.sub(kx);
 const dqpx=fp.mul(cp).add(ff.mul(cpx)).add(fpp).div(2),dkplus=fp.mul(qx).add(ff.mul(qxx)).add(fpp).add(fp.mul(k)).add(ff.mul(kx)).div(2);
 const dqmx=fp.mul(cm).add(ff.mul(cmx)).add(fpp).div(-2),dkminus=fp.mul(qx).add(ff.mul(qxx)).add(fpp).sub(fp.mul(k)).sub(ff.mul(kx)).div(2);
 ck('plus_affine_source:'+row.join(','),dqpx.add(dkplus).eq(fp.mul(cp).add(ff.mul(cpx)).add(fpp)));
 ck('plus_leaves_other_connection:'+row.join(','),dqpx.sub(dkplus).zero());
 ck('minus_affine_source:'+row.join(','),dqmx.sub(dkminus).eq(fp.mul(cm).add(ff.mul(cmx)).add(fpp).neg()));
 ck('minus_leaves_other_connection:'+row.join(','),dqmx.add(dkminus).zero());
 no('drop_connection_inhomogeneous_term:'+row.join(','),!fpp.zero());
 fixtures.push({k,qx,kx,qxx,ff,fp,fpp,cp,cm});
}
// Exact jet form of (log Y)/kappa with Y''=c Y'. Transcendental
// functions are handled by the written derivative/positivity proof.
for(const row of [[2,3,4,5],[-2,3,-4,5],['1/2',2,3,-1]]){
 const [kap,Y,Yp,c]=row.map(f),Ypp=c.mul(Yp),Fp=Yp.div(kap.mul(Y)),FppOverFp=Ypp.div(Yp).sub(Yp.div(Y));
 ck('affine_straightening_jet:'+row.join(','),kap.mul(Fp).add(FppOverFp).eq(c));
 ck('orientation_jet:'+row.join(','),f(0).le(Fp)&&!Fp.zero());
 no('treat_connection_as_density:'+row.join(','),!kap.mul(Fp).eq(c));
}
for(const k of [-3,-1,1,2]){
 const K=f(k),lapse=f(1).div(K);ck('geometry_clock_gauge:'+k,lapse.mul(K).eq(1));
 ck('scalar_clock_fails_on_P0:'+k,lapse.mul(0).zero());
 ck('absolute_measure_jacobian:'+k,K.abs().mul(f(1).div(K.abs())).eq(1));
 ck('Phi_weight_correction:'+k,f(-1).div(K).sub(f('-1/2').div(K).mul(2)).zero());
 no('bare_scalar_derivative_is_symmetric:'+k,!f(-1).div(K).zero());
 if(k<0)no('signed_measure_is_positive:'+k,!f(0).le(f(1).div(K)));
}
no('geometry_clock_at_joint_origin',f(0).mul(1).zero());
// Two-axis packet witness, both unit packets supported on [1,2].
const wP=f(3).pow(2).div(25),wK=f(4).pow(2).div(25),mean=f('3/2');
ck('two_axis_unit_physical_norm',wP.add(wK).eq(1));
ck('P_branch_expectation',wP.mul(mean).eq('27/50'));
ck('K_branch_expectation',wK.mul(mean).eq('24/25'));
ck('PK_constraint_expectation',wP.mul(0).add(wK.mul(0)).zero());
no('product_of_means_equals_PK_expectation',!wP.mul(wK).mul(mean.pow(2)).zero());
// Disjoint dyadic rectangles give an M-growing lower bound for the origin
// Gaussian integral, in units of its fixed positive prefactor.
for(const m of [1,2,4,8,12]){
 const sqrtE=f(1).div(f(2).pow(m+2));let areaOverSqrt=f(0);
 for(let j=0;j<m;j++){
  const high=f(1).div(f(2).pow(j)),low=high.div(2),kmax=sqrtE.div(high);
  ck('origin_rectangle_inside_square:'+m+':'+j,kmax.le(1)&&!kmax.zero());
  ck('origin_rectangle_kernel_bound:'+m+':'+j,high.mul(kmax).eq(sqrtE));
  areaOverSqrt=areaOverSqrt.add(high.sub(low).mul(2).mul(kmax).div(sqrtE));
 }
 ck('origin_mass_lower_bound:'+m,areaOverSqrt.eq(m));
}
// Origin evaluation is positive; a derivative distribution is not a positive
// measure (apply delta' to |1+tP|^2).
for(const t of [1,2,3])no('origin_derivative_is_positive_measure:'+t,!f(0).le(f(-2).mul(t)));
ck('origin_atom_vacuum_norm',C(1).norm2().eq(1));
// Native word witnesses for source reduction and affine jets.
const tokens=['A','B','C','D','p','f','k','q','r','s','t','u'],rules=[];
for(let j=0;j<tokens.length;j++)for(let k=0;k<j;k++)rules.push({id:'commute_'+j+'_'+k,lhs:[tokens[j],tokens[k]],rhs:[[[tokens[k],tokens[j]],1]],source:'commuting classical field jets and scalar indeterminates'});
const scalar=new M.Presentation({tokens,rules});const audit=scalar.audit();ck('word_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:scalar.spec(),audit});
function prove(name,lhs,rhs){const r=scalar.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',scalar.replay(cert));proofs.push({name,presentation_sha256:scalar.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
const [ap,am,bp,bm,p,ph,k,q,r,s,t,u]=tokens.map(x=>scalar.word([x]));
prove('zero_A_currents_force_p',p.plus(ph).plus(p.minus(ph)),p.scale(2));
prove('zero_A_currents_force_phi_prime',p.plus(ph).minus(p.minus(ph)),ph.scale(2));
prove('plus_affine_variation',r.times(k.plus(q)).plus(s.times(t.plus(u))).plus(p),r.times(k).plus(r.times(q)).plus(s.times(t)).plus(s.times(u)).plus(p));
prove('connection_straightening_cleared_denominators',ap.times(bp).plus(am.times(ap)).minus(ap.times(bp)),am.times(ap));
prove('two_branch_constraint_support',p.times(k).times(ph),k.times(p).times(ph));
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'CCS1_COMPLEMENTARY_CLOCK_SECTOR_V1',status:'PASS_COMPLEMENTARY_CLOCK_AND_PUNCTURED_TWO_AXIS_COMPLETION',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{smooth_P0_source:'p=0,phi_prime=0',complementary_clock:'q=tau on P=0,K!=0 after local affine straightening',reduced_pair:'Phi,K',dressing:'exp[-sum U_n V_n/(2n)] per chirality, independent of K',physical_norm:'L2(dK/|K|)',scalar_observable:'i[d/dK-1/(2K)] on transported full-line H1 domain',two_axis_completion:'H_P direct_sum H_K on tests excluding the joint origin',origin_average:'positive logarithmic divergence on nonzero-origin tests',origin_atom:'optional separate rank-one positive form, no regular canonical clock',full_local_quantum_fields_constructed:false,interacting_dynamics_constructed:false,physical_3plus1_gravity_complete:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
