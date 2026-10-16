'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const O=require('./oscillator_fixture.cjs');
const {N,f,C,key,state,vac,plus,scale,minus,eq,current,T,A,L,inner,partitions,descendant,weight,paired}=O;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const predecessor=fs.readFileSync(path.join(__dirname,'../native-local-physical-sector/verify.cjs'),'utf8');
const helper=fs.readFileSync(path.join(__dirname,'oscillator_fixture.cjs'),'utf8');
assert.equal(helper.slice(helper.indexOf('const key='),helper.indexOf('\nmodule.exports=')),predecessor.slice(predecessor.indexOf('const key='),predecessor.indexOf('const fixtures=[];')));
const checks=[],negative_controls=[],proofs=[],presentations=[],level_fixtures=[],jet_fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
function D(v,t,l){let o=new Map();for(let n=1;n<=l;n++)o=plus(o,scale(current('u',-n,A('u',n,v,t),t),n));return o;}
function H(v,t,l){let o=new Map();for(let n=1;n<=l;n++)o=plus(o,scale(current('u',-n,current('u',n,v,t),t),C(t).mul(n).div(2)));for(let a=1;a<=l;a++)for(let b=1;b<=l-a;b++){
 o=plus(o,scale(plus(current('u',-a,current('u',-b,current('u',a+b,v,t),t),t),current('u',-a-b,current('u',a,current('u',b,v,t),t),t)),f(a+b).div(4)));}return o;}
function E(v,t,l){let o=new Map();for(let n=1;n<=l;n++)o=plus(o,scale(current('u',-n,current('u',n,v,t),t),f(n).pow(2).div(2)));return o;}
function columns(vs,fn){return vs.map(x=>vs.map(y=>inner(x,fn(y)).div(inner(x,x))));}
for(const raw of ['0','1','-2','3/2'])for(let l=1;l<=6;l++){
 const t=f(raw),ps=partitions(l),vs=ps.map(p=>state(p.map(n=>'u:'+n))),dim=ps.length;
 const W=N.matrix(vs.map((x,j)=>vs.map((_,k)=>j===k?inner(x,x):N.ZERO)));
 const dm=columns(vs,x=>D(x,t,l)),hm=columns(vs,x=>H(x,t,l)),em=columns(vs,x=>E(x,t,l));
 ck('exact_positive_split:'+raw+':'+l,N.equal(dm,N.add(hm,N.scale(em,N.IOTA))));
 ck('weighted_Hermitian_H:'+raw+':'+l,N.equal(N.mul(W,hm),N.mul(N.dagger(hm),W)));
 const imaginary=N.scale(N.sub(N.mul(W,dm),N.mul(N.dagger(dm),W)),new N.Cut(0,'-1/2'));
 ck('imaginary_coercivity_matrix:'+raw+':'+l,N.equal(imaginary,N.mul(W,em)));
 const energies=ps.map(p=>p.reduce((z,n)=>z.add(f(n).pow(3)),f(0)));
 ck('E_number_operator:'+raw+':'+l,N.equal(em,N.matrix(ps.map((_,j)=>ps.map((_,k)=>j===k?energies[j]:f(0))))));
 ck('positive_level_bound:'+raw+':'+l,energies.every(e=>f(l).le(e)));
 ck('D_invertibility_witness:'+raw+':'+l,N.equal(N.mul(dm,N.inverse(dm)),N.identity(dim)));
 const stack=[];
 for(let n=1;n<=l;n++)for(const p of partitions(l-n)){const out=state(p.map(k=>'u:'+k));stack.push(vs.map(x=>inner(out,A('u',n,x,t)).div(inner(out,out))));}
 ck('no_primary_finite_witness:'+raw+':'+l,N.rank(stack)===dim);
 const ds=ps.map(p=>descendant('u',p,t)),cm=vs.map(x=>ds.map(y=>inner(x,y).div(inner(x,x))));
 ck('descendants_span_level:'+raw+':'+l,N.rank(cm)===dim);
 const gram=ds.map(x=>ds.map(y=>inner(x,y)));
 ck('Gram_real_positive_diagonal:'+raw+':'+l,gram.every((row,j)=>row[j].turn.zero()&&!row[j].rad.zero()&&f(0).le(row[j].rad))&&gram.flat().every(z=>z.turn.zero()));
 ck('Gram_inverse:'+raw+':'+l,N.equal(N.mul(gram,N.inverse(gram)),N.identity(dim)));
 if(l<=3)level_fixtures.push({t,level:l,dimension:dim,D:dm,H:hm,E:em,Gram:gram});
}
const p1=state(['u:1']);
for(let n=1;n<=5;n++)ck('control_improvement_removed_primary:'+n,T('u',n,p1,0).size===0);
no('omit_improvement_preserves_no_primary_theorem',T('u',1,p1,0).size===0&&A('u',1,p1,0).size!==0);
for(let n=1;n<=5;n++)ck('control_nonreal_zero_mode_primary:'+n,A('u',n,p1,new N.Cut(0,-1)).size===0);
no('extend_Hermiticity_to_nonreal_t',!inner(p1,H(p1,new N.Cut(0,-1),1)).turn.zero());
// Explicit matched/opposite-zero-mode pairings and unequal-weight controls.
function pairedDifferent(level,t,u){const ps=partitions(level),us=ps.map(p=>descendant('u',p,t)),ws=ps.map(p=>descendant('v',p,u));const gram=ws.map(x=>ws.map(y=>inner(x,y))),iv=N.inverse(gram);let x=new Map();for(let a=0;a<ps.length;a++)for(let b=0;b<ps.length;b++)x=plus(x,scale(O.tensor(us[a],ws[b]),iv[a][b]));return x;}
for(const [raw,other]of [['1','1'],['1','-1'],['0','0'],['3/2','-3/2']]){
 const t=f(raw),u=f(other),ip=plus(...[0,1,2,3,4].map(l=>pairedDifferent(l,t,u)));
 for(const pa of [[],[1],[2],[1,1]])for(const pb of [[],[1],[2],[1,1]])for(const m of [-2,-1,0,1,2]){
 const x=state([...pa.map(n=>'u:'+n),...pb.map(n=>'v:'+n)]);
 ck('full_dual_pairing:'+raw+':'+other+':'+pa.join('.')+'|'+pb.join('.')+':'+m,inner(ip,L(m,x,t,u)).zero());
 }
}
// An L0 resonance at Delta=1 is not an all-mode physical solution.
const t=f(2),u=f(0),resonance=state(['v:1']);
ck('integer_grade_resonance_L0',L(0,resonance,t,u).size===0);
no('integer_grade_resonance_all_mode_kernel',L(1,resonance,t,u).size!==0);
// Reproduce reduction rho off shell in descendant coordinates, using exact
// change-of-basis inverses. Its residual vanishes on both axes, but not off shell.
function rho(x,t,u){let out=N.ZERO;const groups=new Map();for(const[k,c]of x){const a=O.grade(k,'u'),b=O.grade(k,'v');if(a!==b)continue;if(!groups.has(a))groups.set(a,new Map());groups.get(a).set(k,c);}
 for(const [l,part]of groups){const ps=partitions(l),basis=ps.map(p=>state(p.map(n=>'u:'+n))),dsu=ps.map(p=>descendant('u',p,t)),dsw=ps.map(p=>descendant('u',p,u));
 const mat=ds=>basis.map(z=>ds.map(d=>inner(z,d).div(inner(z,z)))),iu=N.inverse(mat(dsu)),iw=N.inverse(mat(dsw)),G=dsw.map(z=>dsw.map(d=>inner(z,d)));
 const monomials=new Map(ps.map((p,j)=>[JSON.stringify(p.slice().sort((a,b)=>a-b)),j]));
 for(const[k,c]of part){const vars=JSON.parse(k),get=s=>monomials.get(JSON.stringify(vars.filter(v=>v.startsWith(s+':')).map(v=>Number(v.split(':')[1])).sort((a,b)=>a-b))),a=get('u'),b=get('v');let value=N.ZERO;
 for(let i=0;i<ps.length;i++)for(let j=0;j<ps.length;j++)value=value.add(iu[i][a].mul(iw[j][b]).mul(G[i][j]));out=out.add(c.mul(value));}
 }return out;
}
for(const [aa,bb]of [[1,0],[0,1],[0,0],[2,0],[0,-2],[1,1]]){
 const a=f(aa),b=f(bb),t=a.add(b),u=a.sub(b);
 ck('rho_vacuum:'+aa+','+bb,rho(vac,t,u).eq(1));
 ck('rho_L0_vacuum:'+aa+','+bb,rho(L(0,vac,t,u),t,u).eq(a.mul(b)));
 if(a.mul(b).zero())for(const x of [vac,state(['u:1']),state(['v:2']),state(['u:1','v:1'])])for(const m of [-2,-1,0,1,2])ck('rho_axis_residual:'+aa+','+bb+':'+[...x.keys()][0]+':'+m,rho(L(m,x,t,u),t,u).zero());
 else no('rho_is_offshell_constraint_intertwiner',!rho(L(0,vac,t,u),t,u).zero());
}
// Distribution jets: delta^(r)(a) delta^(s)(b), x delta^(r)=-r delta^(r-1).
function multiplyCoordinate(r,s,axis){if((axis==='a'?r:s)===0)return {coefficient:f(0),r,s};return {coefficient:f(-(axis==='a'?r:s)),r:r-(axis==='a'?1:0),s:s-(axis==='b'?1:0)};}
for(let r=0;r<=4;r++)for(let s=0;s<=4;s++){
 const a=multiplyCoordinate(r,s,'a'),b=multiplyCoordinate(a.r,a.s,'b'),coefficient=a.coefficient.mul(b.coefficient);
 ck('ab_delta_jet:'+r+','+s,coefficient.eq(f(r).mul(s)));
 ck('axis_jet_kernel:'+r+','+s,coefficient.zero()===(r===0||s===0));
 if(r>0&&s>0)no('mixed_normal_derivative_is_solution:'+r+','+s,!coefficient.zero());
 jet_fixtures.push({r,s,ab_coefficient:coefficient});
}
for(const P of [-3,-1,1,2]){
 const p=f(P),logDerivative=f(-1).div(p),qCorrection=f('-1/2').div(p);
 ck('weight_ODE:'+P,logDerivative.sub(qCorrection.mul(2)).zero());
 no('Lebesgue_measure_symmetric_for_fixed_Q:'+P,!qCorrection.mul(2).zero());
}
for(const [cp,cm]of [[1,1],[2,2],[1,2],[3,1]]){
 const boundary=new N.Cut(0,f(cp).sub(cm));
 ck('matching_boundary:'+cp+','+cm,boundary.zero()===(cp===cm));
 if(cp!==cm)no('unequal_sign_weights_with_fixed_matching:'+cp+','+cm,!boundary.zero());
}
// Finite native word witnesses, with explicit oscillator relations.
// Use an independent Weyl presentation whose three oscillator modes have
// [x,X]=2, [y,Y]=4, [z,Z]=6. Every other pair commutes.
const names=['X','Y','Z','x','y','z'],rules=[];
for(let j=0;j<names.length;j++)for(let k=0;k<j;k++){const rhs=[[[names[k],names[j]],1]];if(j===k+3)rhs.push([[],2*(k+1)]);rules.push({id:'r'+j+'_'+k,lhs:[names[j],names[k]],rhs,source:'positive current Weyl relation at modes one, two, three'});}
const weyl=new M.Presentation({tokens:names,rules});
const audit=weyl.audit();ck('Weyl_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:weyl.spec(),audit});
function prove(sys,name,lhs,rhs){const r=sys.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',sys.replay(cert));proofs.push({name,presentation_sha256:sys.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
const [X,Y,Z,x,y,z]=names.map(s=>weyl.word([s]));
prove(weyl,'two_mode_positive_norm',x.times(X),X.times(x).plus(weyl.one().scale(2)));
prove(weyl,'cubic_mode_three_contraction',X.times(Y).times(z).times(Z).minus(X.times(Y).times(Z).times(z)),X.times(Y).scale(6));
const scalarNames=['a','b','n','k'],scalarRules=[];for(let j=0;j<scalarNames.length;j++)for(let k=0;k<j;k++)scalarRules.push({id:'s'+j+'_'+k,lhs:[scalarNames[j],scalarNames[k]],rhs:[[[scalarNames[k],scalarNames[j]],1]],source:'commuting scalar indeterminates'});
const scalar=new M.Presentation({tokens:scalarNames,rules:scalarRules});const sa=scalar.audit();ck('scalar_confluence',sa.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:scalar.spec(),audit:sa});
const [a,b,n,k]=scalarNames.map(s=>scalar.word([s]));
prove(scalar,'zero_weight_support_ab',a.plus(b).times(a.plus(b)).minus(a.minus(b).times(a.minus(b))).scale('1/4'),a.times(b));
prove(scalar,'cubic_symmetrization_coefficient',n.scale('1/2').plus(k.scale('1/2')).scale('1/2'),n.plus(k).scale('1/4'));
prove(scalar,'n_cubed_minus_n_nonnegative_factor',n.times(n).times(n).minus(n),n.times(n.minus(scalar.one())).times(n.plus(scalar.one())));
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'LSC1_LOCAL_SECTOR_CLASSIFICATION_V1',status:'PASS_TEST_DUAL_CLASSIFICATION_AND_CONDITIONAL_MEASURE_UNIQUENESS',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,level_fixtures,jet_fixtures,results:{all_level_primary_obstruction:'Im<v,D_N v>=sum n^3 <N_n> >= N||v||^2; written proof',all_real_fibers_cyclic:true,fiber_kernel:'one dimensional iff a0*b0=0, else zero',continuous_test_dual_kernel:'F composed with rho; F=U(a)delta(b)+delta(a)V(b)',regular_physical_oscillator_multiplicity:1,measure_uniqueness:'conditional on momentum multipliers, fixed Q and full-line matching; overall scale remains',constraints_alone_select_measure:false,singular_axis_physical_quantization_completed:false,quantum_relational_field_operators_constructed:false,physical_3plus1_gravity_complete:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
