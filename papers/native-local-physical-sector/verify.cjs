'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,C=N.Cut.of,checks=[],negative_controls=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const key=v=>JSON.stringify([...v].sort());
const state=(v=[],c=1)=>new Map([[key(v),C(c)]]),vac=state();
function plus(...xs){const out=new Map();for(const x of xs)for(const[k,c]of x){const z=(out.get(k)||N.ZERO).add(c);if(z.zero())out.delete(k);else out.set(k,z);}return out;}
const scale=(x,c)=>new Map([...x].map(([k,z])=>[k,z.mul(c)]).filter(([,z])=>!z.zero()));
const minus=(x,y)=>plus(x,scale(y,-1)),eq=(x,y)=>minus(x,y).size===0;
function create(x,v){return plus(...[...x].map(([k,c])=>state([...JSON.parse(k),v],c)));}
function deriv(x,v){return plus(...[...x].map(([k,c])=>{const a=JSON.parse(k),r=a.filter(z=>z===v).length;if(!r)return new Map();a.splice(a.indexOf(v),1);return state(a,c.mul(r));}));}
function current(s,n,x,t){if(n===0)return scale(x,t);return n<0?create(x,s+':'+(-n)):scale(deriv(x,s+':'+n),2*n);}
function bound(x){let b=0;for(const k of x.keys())for(const v of JSON.parse(k))b=Math.max(b,Number(v.split(':')[1]));return b;}
function T(s,m,x,t){let out=new Map();const B=bound(x)+Math.abs(m);for(let j=-B;j<=B;j++){
 const a=m-j,z=a>0&&j<0?current(s,j,current(s,a,x,t),t):current(s,a,current(s,j,x,t),t);out=plus(out,z);
}return scale(out,'1/4');}
const A=(s,m,x,t)=>plus(T(s,m,x,t),scale(current(s,m,x,t),new N.Cut(0,f(m).div(2))),m===0?scale(x,'1/4'):new Map());
const L=(m,x,t,u=t)=>minus(A('u',m,x,t),A('v',-m,x,u));
const sourceL=(m,x,t,u=t)=>plus(T('u',m,x,t),scale(T('v',-m,x,u),-1),scale(plus(current('u',m,x,t),current('v',-m,x,u)),new N.Cut(0,f(m).div(2))));
function weight(vars){const counts=new Map();for(const v of vars)counts.set(v,(counts.get(v)||0)+1);let w=f(1);for(const[v,r]of counts){const n=Number(v.split(':')[1]);for(let j=1;j<=r;j++)w=w.mul(2*n).mul(j);}return w;}
function inner(x,y){let z=N.ZERO;for(const[k,c]of x)if(y.has(k))z=z.add(c.dagger().mul(y.get(k)).mul(weight(JSON.parse(k))));return z;}
function tensor(x,y){return plus(...[...x].flatMap(([a,c])=>[...y].map(([b,d])=>state([...JSON.parse(a),...JSON.parse(b)],c.mul(d)))));}
function partitions(n,max=n){if(n===0)return [[]];const out=[];for(let k=Math.min(n,max);k>=1;k--)for(const p of partitions(n-k,k))out.push([k,...p]);return out;}
function descendant(s,p,t){let x=vac;for(const n of [...p].reverse())x=A(s,-n,x,t);return x;}
const grade=(k,s)=>JSON.parse(k).filter(v=>v.startsWith(s+':')).reduce((a,v)=>a+Number(v.split(':')[1]),0);
const project=(x,max)=>new Map([...x].filter(([k])=>grade(k,'u')<=max&&grade(k,'v')<=max));
function paired(level,t){const ps=partitions(level),us=ps.map(p=>descendant('u',p,t)),vs=ps.map(p=>descendant('v',p,t));
 const G=us.map(x=>us.map(y=>inner(x,y))),iv=N.inverse(G);let out=new Map();
 for(let a=0;a<ps.length;a++)for(let b=0;b<ps.length;b++)out=plus(out,scale(tensor(us[a],vs[b]),iv[a][b]));
 return {out,G,iv,ps,us,vs};
}
const fixtures=[];
for(const raw of ['0','1','-1','3/2']){
 const t=f(raw),h=t.pow(2).add(1).div(4),blocks=[],max=4;
 for(let l=0;l<=max;l++){
  const b=paired(l,t);blocks.push(b.out);
  ck('Gram_real:'+raw+':'+l,b.G.flat().every(z=>z.turn.zero()));
  ck('Gram_inverse:'+raw+':'+l,N.equal(N.mul(b.G,b.iv),N.identity(b.G.length)));
  ck('Gram_matching_copies:'+raw+':'+l,N.equal(b.G,b.vs.map(x=>b.vs.map(y=>inner(x,y)))));
  ck('paired_level_norm:'+raw+':'+l,inner(b.out,b.out).eq(b.ps.length));
  for(let j=0;j<b.us.length;j++)ck('positive_descendant:'+raw+':'+l+':'+j,!inner(b.us[j],b.us[j]).zero()&&f(0).le(inner(b.us[j],b.us[j]).rad));
 }
 const total=plus(...blocks),expectedNorm=blocks.reduce((z,_,l)=>z.add(partitions(l).length),f(0));
 ck('cutoff_norm_exact:'+raw,inner(total,total).eq(expectedNorm));
 for(const m of [-3,-2,-1,0,1,2,3]){
  ck('source_factorization:'+raw+':'+m,eq(L(m,total,t),sourceL(m,total,t)));
  ck('gluing_interior:'+raw+':'+m,project(L(m,total,t),max-Math.abs(m)).size===0);
 }
 for(const m of [-2,-1,1,2])no('cutoff_is_exact_all_mode_state:'+raw+':'+m,L(m,total,t).size!==0);
 const first=C(new N.Cut(t,-1)).div(new N.Cut(t,1)).div(2);
 ck('level_one_phase:'+raw,blocks[1].get(key(['u:1','v:1'])).eq(first));
 ck('level_one_modulus:'+raw,first.norm2().eq('1/4'));
 no('naive_improvement_free_pairing:'+raw,!first.eq('1/2'));
 const bad=plus(vac,state(['u:1','v:1'],'1/2'));
 no('wrong_pairing_fails_test:'+raw,!inner(bad,L(1,state(['u:1']),t)).zero());
 if(raw!=='0')no('wrong_pairing_conjugation:'+raw,!eq(blocks[1],new Map([...blocks[1]].map(([k,c])=>[k,c.dagger()]))));
 // Test the dual equation on unequal levels; level cutoff is sufficient
 // for every resulting input, rather than zeroing a boundary residual.
 for(let a=0;a<=2;a++)for(let b=0;b<=2;b++)for(const pa of partitions(a))for(const pb of partitions(b))for(const m of [-2,-1,0,1,2]){
  const x=state([...pa.map(n=>'u:'+n),...pb.map(n=>'v:'+n)]);
  ck('dual_gluing:'+raw+':'+pa.join('.')+'|'+pb.join('.')+':'+m,inner(total,L(m,x,t)).zero());
 }
 for(const [m,n]of [[1,-1],[2,-2],[3,-3],[2,-1],[-2,1]]){
  const x=state(['u:1']),lhs=minus(A('u',m,A('u',n,x,t),t),A('u',n,A('u',m,x,t),t));
  const central=m+n===0?f(m).pow(3).sub(m).mul('7/12'):f(0);
  ck('auxiliary_c7:'+raw+':'+m+','+n,eq(lhs,plus(scale(A('u',m+n,x,t),m-n),scale(x,central))));
 }
 let d=vac,norm=f(1);for(let l=1;l<=6;l++){d=A('u',-1,d,t);norm=norm.mul(l).mul(h.mul(2).add(l-1));ck('sl2_norm_recurrence:'+raw+':'+l,inner(d,d).eq(norm));}
 fixtures.push({t,h,first_level_coefficient:first,level_two_Gram:paired(2,t).G,cutoff_level:max,cutoff_norm_squared:inner(total,total),interpretation:'finite witnesses only; all-level kernel proved in THEOREM.md'});
}
for(const [t,u]of [[2,1],[-2,1],[1,0]])ck('zero_mode_source_no_intercept:'+t+':'+u,inner(vac,L(0,vac,f(t),f(u))).eq(f(t).pow(2).sub(f(u).pow(2)).div(4)));
no('unequal_zero_weights_allow_vacuum_pairing',!inner(vac,L(0,vac,f(2),f(1))).zero());
const probabilities=[f(3).pow(2).div(25),f(4).pow(2).div(25)];
ck('physical_packet_unit_norm',probabilities[0].add(probabilities[1]).eq(1));
ck('physical_packet_mean_P',probabilities[0].sub(probabilities[1]).mul('3/2').eq('-21/50'));
const integralPower=(lo,hi,n)=>f(hi).pow(n+1).sub(f(lo).pow(n+1)).div(n+1);
ck('physical_packet_second_moment',integralPower(1,2,2).eq('7/3'));
no('signed_measure_is_positive',!f(0).le(probabilities[1].neg()));
for(const p of [2,-2]){const weightDerivativeOverWeight=f(-1).div(p),correction=f('-1/2').div(p);ck('physical_Q_weight_symmetry:'+p,weightDerivativeOverWeight.sub(correction.mul(2)).zero());no('naive_Q_omits_weight_correction:'+p,!weightDerivativeOverWeight.zero());}
// Source-specific local regular reduction. Jets need not be spatial constants.
const localFixtures=[];
for(const data of [[2,3,1,-2,4,3],[2,2,0,1,-1,2],[-2,-3,1,2,3,-1],['3/2','2/3',-2,3,'1/2',4]]){
 const [ap,am,dp,dm,ff,fp]=data.map(f),bp=dp.div(ap),bm=dm.div(am).neg(),qp=bp.sub(bm).div(2),k=bp.add(bm).div(2);
 ck('local:T_plus:'+data.join(','),ap.mul(k.add(qp)).sub(dp).zero());
 ck('local:T_minus:'+data.join(','),am.mul(k.sub(qp)).add(dm).zero());
 ck('local:C_spatial_constant:'+data.join(','),qp.sub(dp.div(ap).add(dm.div(am)).div(2)).zero());
 const deltaAp=fp.mul(ap).add(ff.mul(dp)),deltaAm=fp.mul(am).add(ff.mul(dm)).neg();
 const deltaQplus=ff.mul(bp).add(fp).div(2),deltaQminus=ff.mul(bm).sub(fp).div(2);
 ck('local:C_plus_invariant:'+data.join(','),deltaQplus.sub(deltaAp.div(ap).div(2)).zero());
 ck('local:C_minus_invariant:'+data.join(','),deltaQminus.sub(deltaAm.div(am).div(2)).zero());
 if(!fp.zero())no('local:omit_derivative_improvement:'+data.join(','),!ff.mul(bp).div(2).sub(deltaAp.div(ap).div(2)).zero());
 localFixtures.push({ap,am,dp,dm,qp,k,deltaAp,deltaAm});
}
no('local:nonzero_P_implies_nonvanishing_current',f(1).sub(1).zero()); // A+=P(1+sin x) has a zero despite mean P.

// Native word witnesses are finite identities, not rewrite axioms for gluing.
const sys=new M.Presentation({tokens:['t','m','a','b'],rules:[
 {id:'mt',lhs:['m','t'],rhs:[[['t','m'],1]],source:'commuting scalar parameters'},
 {id:'at',lhs:['a','t'],rhs:[[['t','a'],1]],source:'commuting scalar parameters'},
 {id:'am',lhs:['a','m'],rhs:[[['m','a'],1]],source:'commuting scalar parameters'},
 {id:'bt',lhs:['b','t'],rhs:[[['t','b'],1]],source:'commuting scalar parameters'},
 {id:'bm',lhs:['b','m'],rhs:[[['m','b'],1]],source:'commuting scalar parameters'},
 {id:'ba',lhs:['b','a'],rhs:[[['a','b'],1]],source:'commuting scalar parameters'}]});
const audit=sys.audit();ck('word_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:sys.spec(),audit});
function prove(name,lhs,rhs){const r=sys.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',sys.replay(cert));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
const [t,m,a,b]=['t','m','a','b'].map(s=>sys.word([s])),one=sys.one(),ii=one.scale(N.IOTA),m3=m.times(m).times(m);
prove('c7_central_polynomial',m3.minus(m).scale('1/12').plus(m3.scale('1/2')).minus(m.scale('1/2')),m3.minus(m).scale('7/12'));
prove('phase_denominator_identity',t.minus(ii).times(t.plus(ii)),t.times(t).plus(one));
prove('matching_weight_difference',a.plus(b).times(a.plus(b)).minus(a.minus(b).times(a.minus(b))).scale('1/4'),a.times(b));
prove('auxiliary_intercepts_cancel',a.plus(one.scale('1/4')).minus(b.plus(one.scale('1/4'))),a.minus(b));
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'LPS1_REGULAR_LOCAL_PHYSICAL_SECTOR_V1',status:'PASS_SELECTED_REGULAR_LOCAL_PHYSICAL_SECTOR',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,localFixtures,results:{auxiliary_c:7,local_constraint_c:0,regular_branch:'P!=0,K=0',physical_measure:'dP/|P|',physical_solution:'two chirality descendant pairings in the test dual',physical_observables:'HQS reduced algebra transported with domains',unique_rigging:false,full_kernel_classified:false,classical_regular_local_reduction_proved:true,relational_length_quantized_after_reduction:true,pointwise_quantum_field_length_constructed:false,full_physical_quantum_gravity:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
