'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,C=N.Cut.of;
const checks=[],negative_controls=[],presentations=[],proofs=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
// Finite polynomial vectors in the declared infinite oscillator module.
// All coefficients use the unchanged native Cut/F arithmetic. Mode indices
// and array counts are combinatorial integers, not approximate field values.
const key=vars=>JSON.stringify([...vars].sort());
const state=(vars=[],coefficient=1)=>new Map([[key(vars),C(coefficient)]]);
const vacuum=state();
function plus(...xs){const out=new Map();for(const x of xs)for(const [k,c]of x){const v=(out.get(k)||N.ZERO).add(c);if(v.zero())out.delete(k);else out.set(k,v);}return out;}
const scale=(x,c)=>new Map([...x].map(([k,v])=>[k,v.mul(c)]).filter(([,v])=>!v.zero()));
const minus=(x,y)=>plus(x,scale(y,-1));
const equal=(x,y)=>minus(x,y).size===0;
const serialize=x=>[...x].sort(([a],[b])=>a.localeCompare(b)).map(([k,c])=>({creators:JSON.parse(k),coefficient:c}));
function create(x,v){return plus(...[...x].map(([k,c])=>state([...JSON.parse(k),v],c)));}
function derivative(x,v){return plus(...[...x].map(([k,c])=>{const vars=JSON.parse(k),count=vars.filter(t=>t===v).length;if(!count)return new Map();vars.splice(vars.indexOf(v),1);return state(vars,c.mul(f(count)));}));}
const settings=(o={})=>({a0:f(0),b0:f(0),hbar:f(1),r:f(1),s:f(0),gamma:f(0),chir:'+',...o});
function current(kind,n,x,o){
 if(n===0)return scale(x,kind==='a'?o.a0:o.b0);
 if(n<0)return create(x,o.chir+':'+kind+':'+(-n));
 return scale(derivative(x,o.chir+':'+(kind==='a'?'b':'a')+':'+n),f(n).mul(o.hbar));
}
function occupiedBound(x){let b=0;for(const k of x.keys())for(const v of JSON.parse(k))b=Math.max(b,Number(v.split(':')[2]));return b;}
function S(n,x,o,forcedBound=null){
 const bound=forcedBound===null?occupiedBound(x)+Math.abs(n):forcedBound;
 let out=new Map();
 for(let j=-bound;j<=bound;j++){
  const a=n-j;
  // :a_a b_j:; the rightmost operator acts first.
  const term=a>0&&j<0?current('b',j,current('a',a,x,o),o):current('a',a,current('b',j,x,o),o);
  out=plus(out,term);
 }
 return out;
}
function L(n,x,o){return plus(S(n,x,o),scale(plus(scale(current('a',n,x,o),o.r),scale(current('b',n,x,o),o.s)),new N.Cut(0,n)),n===0?scale(x,o.gamma):new Map());}
const F=m=>f(m).pow(3).sub(m).div(6);
const anomaly=(m,o)=>o.hbar.pow(2).mul(F(m)).add(o.hbar.mul(2).mul(o.r).mul(o.s).mul(f(m).pow(3))).sub(o.hbar.mul(2).mul(o.gamma).mul(m));
const commAction=(A,B,x)=>minus(A(B(x)),B(A(x)));
const baseline=settings(),repair=settings({s:f('-1/12'),gamma:f('-1/12')});
const scaledRepair=settings({hbar:f(2),s:f('-1/6'),gamma:f('-1/6'),a0:f(2),b0:f('3/5')});
const seeds=[vacuum,state(['+:a:1']),state(['+:a:1','+:b:2']),state(['+:b:7'])];
const pairs=[[2,-2],[3,-3],[5,-5],[2,-1],[-2,1],[1,2],[-1,-2],[0,2],[3,0]];
const modeExamples=[];
for(const [label,o]of [['original',baseline],['counterterm',repair],['restored_hbar',scaledRepair]]){
 for(const [j,x]of seeds.entries())for(const [m,n]of pairs){
  const lhs=commAction(y=>L(m,y,o),y=>L(n,y,o),x);
  const rhs=plus(scale(L(m+n,x,o),o.hbar.mul(f(m).sub(n))),m+n===0?scale(x,anomaly(m,o)):new Map());
  ck('mode_commutator:'+label+':seed'+j+':'+m+','+n,equal(lhs,rhs));
  if(j===0&&m+n===0)modeExamples.push({representation:label,m,n,result:serialize(lhs),scalar_defect:anomaly(m,o)});
 }
}
// Basic currents and primary transformation are computed by action, not by
// asserting the desired Virasoro relation as a rewrite rule.
for(const n of [1,2,5]){
 for(const kind of ['a','b']){
  const other=kind==='a'?'b':'a';
  ck('current_cross:'+kind+':'+n,equal(commAction(x=>current(kind,n,x,baseline),x=>current(other,-n,x,baseline),vacuum),scale(vacuum,n)));
  ck('current_null:'+kind+':'+n,commAction(x=>current(kind,n,x,baseline),x=>current(kind,-n,x,baseline),vacuum).size===0);
 }
}
for(const [m,n]of [[2,-3],[-2,1],[0,-2]])for(const kind of ['a','b']){
 const x=state(['+:a:1','+:b:2']);
 ck('Sugawara_current:'+kind+':'+m+','+n,equal(commAction(y=>S(m,y,baseline),y=>current(kind,n,y,baseline),x),scale(current(kind,m+n,x,baseline),f(n).neg())));
}
const contractions=[];
for(let m=1;m<=12;m++){
 let sum=f(0);for(let j=1;j<m;j++)sum=sum.add(f(j).mul(f(m).sub(j)));
 ck('finite_contraction_sum:'+m,sum.eq(F(m)));
 ck('vacuum_contraction:'+m,equal(commAction(x=>S(m,x,baseline),x=>S(-m,x,baseline),vacuum),scale(vacuum,sum)));
 contractions.push({m,double_contraction:sum});
}
ck('central_examples',F(1).eq(0)&&F(2).eq(1)&&F(3).eq(4));
no('zero_mode_intercept_repairs_original',!F(2).div(4).eq(F(3).div(6)));
no('drop_double_contractions',!equal(commAction(x=>S(2,x,baseline),x=>S(-2,x,baseline),vacuum),new Map()));
no('fixed_mode_zero_truncation_certifies_central_term',!equal(commAction(x=>S(2,x,baseline,0),x=>S(-2,x,baseline,0),vacuum),scale(vacuum,F(2))));
no('reusing_input_bound_for_composed_operators',!equal(S(2,S(-2,vacuum,baseline),baseline,0),S(2,S(-2,vacuum,baseline),baseline)));
// Both chiral copies: central cancellation in DD and HH is insufficient.
const plusOptions=settings(),minusOptions=settings({chir:'-'});
const H=(m,x)=>plus(L(m,x,plusOptions),L(-m,x,minusOptions));
const D=(m,x)=>minus(L(m,x,plusOptions),L(-m,x,minusOptions));
for(const [j,x]of [vacuum,state(['+:a:1','-:b:2'])].entries())for(const [m,n]of pairs){
 ck('chiral_commute:'+j+':'+m+','+n,commAction(y=>L(m,y,plusOptions),y=>L(n,y,minusOptions),x).size===0);
 ck('DD:'+j+':'+m+','+n,equal(commAction(y=>D(m,y),y=>D(n,y),x),scale(D(m+n,x),f(m).sub(n))));
 ck('HH:'+j+':'+m+','+n,equal(commAction(y=>H(m,y),y=>H(n,y),x),scale(D(m+n,x),f(m).sub(n))));
 ck('DH:'+j+':'+m+','+n,equal(commAction(y=>D(m,y),y=>H(n,y),x),plus(scale(H(m+n,x),f(m).sub(n)),m+n===0?scale(x,F(m).mul(2)):new Map())));
}
no('equal_chiral_central_charges_cancel_all_brackets',!equal(commAction(x=>D(2,x),x=>H(-2,x),vacuum),scale(H(0,vacuum),4)));
ck('mixed_defect_at_two',equal(minus(commAction(x=>D(2,x),x=>H(-2,x),vacuum),scale(H(0,vacuum),4)),scale(vacuum,2)));
// Exact metric/adjoint diagnostic; not interpreted as positive probabilities.
const gramExamples=[];
for(const n of [1,2,5]){
 const G=N.matrix([[0,n],[n,0]]),v=N.matrix([[1],[-1]]),w=N.matrix([[1],[1]]);
 const norm=N.mul(N.mul(N.dagger(v),G),v)[0][0],positive=N.mul(N.mul(N.dagger(w),G),w)[0][0];
 ck('indefinite_Gram_rank:'+n,N.rank(G)===2);
 ck('negative_norm_witness:'+n,norm.eq(f(n).mul(-2)));
 ck('positive_norm_witness:'+n,positive.eq(f(n).mul(2)));
 no('chosen_vacuum_is_positive:'+n,!f(0).le(norm.rad));
 gramExamples.push({mode:n,Gram:G,difference_norm:norm,sum_norm:positive});
}
for(const epsilon of [-1,1])for(const eta of [-1,1])ck('classical_current_cross_sign:'+epsilon+','+eta,f(epsilon).add(eta).eq(epsilon===eta?f(epsilon).mul(2):0));
// Solve the m=1,2 central cancellation equations for (r*s,gamma).
const coefficients=N.matrix([[2,-2],[16,-4]]),rhs=N.matrix([[0],[-1]]);
const solved=N.mul(N.inverse(coefficients),rhs);
ck('counterterm_unique_in_family',N.equal(solved,N.matrix([['-1/12'],['-1/12']])));
for(const r of [1,2,-3])for(const m of [1,2,3,7]){
 const o=settings({r:f(r),s:f('-1/12').div(r),gamma:f('-1/12')});
 ck('counterterm_coefficient:'+r+':'+m,anomaly(m,o).zero());
 if(m>1)no('null_improvement_cancels:'+r+':'+m,!anomaly(m,settings({r:f(r)})).zero());
}
no('cubic_cancellation_without_intercept',!anomaly(2,settings({s:f('-1/12')})).zero());
no('wrong_counterterm_sign',!anomaly(2,settings({s:f('1/12'),gamma:f('-1/12')})).zero());
no('unchanged_HQS_zero_constraint_after_correction',!repair.gamma.mul(4).zero());
ck('homogeneous_shift_in_units_of_pi',repair.gamma.mul(4).eq('-1/3'));
ck('restored_hbar_homogeneous_shift_in_units_of_pi',scaledRepair.gamma.mul(4).eq('-2/3'));
// Alternative positive Fock representation. This changes the vacuum/order,
// not the inherited classical currents or their null linear improvement.
function revCurrent(kind,n,x,o){
 if(n===0)return scale(x,kind==='u'?o.a0.add(o.b0):o.a0.sub(o.b0));
 const creator=kind==='u'?n<0:n>0,v=o.chir+':'+kind+':'+Math.abs(n);
 return creator?create(x,v):scale(derivative(x,v),f(Math.abs(n)).mul(2));
}
function revPart(kind,n,x,o){
 const bound=occupiedBound(x)+Math.abs(n);let out=new Map();
 for(let j=-bound;j<=bound;j++){
  const left=n-j,leftAnn=kind==='u'?left>0:left<0,rightCreator=kind==='u'?j<0:j>0;
  const t=leftAnn&&rightCreator?revCurrent(kind,j,revCurrent(kind,left,x,o),o):revCurrent(kind,left,revCurrent(kind,j,x,o),o);
  out=plus(out,t);
 }
 return scale(out,kind==='u'?'1/4':'-1/4');
}
const revS=(n,x,o)=>plus(revPart('u',n,x,o),revPart('v',n,x,o));
const revL=(n,x,o)=>plus(revS(n,x,o),scale(plus(revCurrent('u',n,x,o),revCurrent('v',n,x,o)),new N.Cut(0,f(n).div(2))));
function positiveWeight(vars){const counts=new Map();for(const v of vars)counts.set(v,(counts.get(v)||0)+1);let w=f(1);for(const [v,count]of counts){const mode=Number(v.split(':')[2]);for(let j=1;j<=count;j++)w=w.mul(f(mode).mul(2).mul(j));}return w;}
function inner(x,y){let out=N.ZERO;for(const [k,c]of x)if(y.has(k))out=out.add(c.dagger().mul(y.get(k)).mul(positiveWeight(JSON.parse(k))));return out;}
const revSeeds=[vacuum,state(['+:u:1']),state(['+:u:1','+:v:2']),state(['+:v:7'])];
for(const [label,o]of [['zero_momenta',baseline],['shared_real_momenta',settings({a0:f(2),b0:f('3/5')})]])for(const [j,x]of revSeeds.entries())for(const [m,n]of pairs){
 ck('reversed:Witt:'+label+':seed'+j+':'+m+','+n,equal(commAction(y=>revL(m,y,o),y=>revL(n,y,o),x),scale(revL(m+n,x,o),f(m).sub(n))));
}
const reverseContractions=[];
for(const m of [1,2,3,5]){
 const up=commAction(x=>revPart('u',m,x,baseline),x=>revPart('u',-m,x,baseline),vacuum);
 const vp=commAction(x=>revPart('v',m,x,baseline),x=>revPart('v',-m,x,baseline),vacuum);
 ck('reversed:u_central:'+m,equal(up,scale(vacuum,F(m).div(2))));
 ck('reversed:v_central:'+m,equal(vp,scale(vacuum,F(m).div(-2))));
 ck('reversed:central_cancel:'+m,plus(up,vp).size===0);
 const Ka=x=>scale(plus(revCurrent('u',m,x,baseline),revCurrent('v',m,x,baseline)),'1/2');
 const Kb=x=>scale(plus(revCurrent('u',-m,x,baseline),revCurrent('v',-m,x,baseline)),'1/2');
 ck('reversed:improvement_null:'+m,commAction(Ka,Kb,vacuum).size===0);
 reverseContractions.push({m,u_central:F(m).div(2),v_central:F(m).div(-2)});
}
for(const m of [-3,-2,-1,1,2,3]){
 const x=revL(m,vacuum,baseline),norm=inner(x,x);
 ck('reversed:adjoint_nonzero_probe:'+m,norm.eq(inner(revL(-m,x,baseline),vacuum)));
 ck('reversed:positive_constraint_image_norm:'+m,!norm.rad.zero()&&f(0).le(norm.rad)&&norm.turn.zero());
}
for(const [j,x]of revSeeds.entries())ck('reversed:positive_monomial:'+j,!inner(x,x).zero()&&f(0).le(inner(x,x).rad));
const revH=(m,x)=>plus(revL(m,x,plusOptions),revL(-m,x,minusOptions));
const revD=(m,x)=>minus(revL(m,x,plusOptions),revL(-m,x,minusOptions));
for(const [m,n]of pairs){const x=state(['+:u:1','-:v:2']);
 ck('reversed:DD:'+m+','+n,equal(commAction(y=>revD(m,y),y=>revD(n,y),x),scale(revD(m+n,x),f(m).sub(n))));
 ck('reversed:HH:'+m+','+n,equal(commAction(y=>revH(m,y),y=>revH(n,y),x),scale(revD(m+n,x),f(m).sub(n))));
 ck('reversed:DH:'+m+','+n,equal(commAction(y=>revD(m,y),y=>revH(n,y),x),scale(revH(m+n,x),f(m).sub(n))));
}
ck('reversed:zero_mode_no_shift',equal(revL(0,vacuum,settings({a0:f(2),b0:f('3/5')})),scale(vacuum,'6/5')));
ck('reversed:negative_mode_grade',equal(revS(0,state(['+:v:3']),baseline),state(['+:v:3'],-3)));
ck('reversed:positive_mode_grade',equal(revS(0,state(['+:u:3']),baseline),state(['+:u:3'],3)));
no('reversed:vacuum_is_physical_all_mode_state',revL(1,vacuum,baseline).size!==0);
no('reversed:source_improvement_can_be_omitted',!equal(revL(1,vacuum,baseline),revS(1,vacuum,baseline)));
no('reversed:S0_is_positive',!f(0).le(inner(state(['+:v:3']),revS(0,state(['+:v:3']),baseline)).rad));
const kernelSeed=state(['+:u:1','+:v:2']),kernelImage=revL(8,kernelSeed,baseline);
ck('reversed:finite_polynomial_kernel_obstruction',kernelImage.get(key(['+:u:1','+:v:2','+:v:3','+:v:5'])).eq('-1/2'));
no('reversed:finite_polynomial_fixture_solves_all_modes',kernelImage.size!==0);
// Serialized finite native word proofs. No Virasoro identity is a rewrite axiom.
function register(name,sys){const audit=sys.audit();ck('words:confluence:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:sys.spec(),audit});}
function prove(sys,name,lhs,rhs){const out=sys.reduce(lhs.minus(rhs),{witness:true});ck('words:'+name,out.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(out.certificate));ck('words:'+name+':replay',sys.replay(certificate));proofs.push({name,presentation_sha256:sys.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
const commutingTokens=['p','k','q','f','u','v','m'];
const commutingRules=[];for(let j=0;j<commutingTokens.length;j++)for(let k=0;k<j;k++)commutingRules.push({id:'commute_'+j+'_'+k,lhs:[commutingTokens[j],commutingTokens[k]],rhs:[[[commutingTokens[k],commutingTokens[j]],1]],source:'Classical field jets / scalar polynomial indeterminate commute'});
const scalar=new M.Presentation({tokens:commutingTokens,rules:commutingRules});register('classical_jets_and_cubic_polynomial',scalar);
const [p,k,q,ph,u,v,z]=commutingTokens.map(t=>scalar.word([t]));
const h=p.times(k).plus(q.times(ph)).minus(v),d=p.times(q).plus(k.times(ph)).minus(u);
for(const sign of [-1,1])prove(scalar,'source_chiral_factorization_'+sign,h.plus(d.scale(sign)),p.plus(ph.scale(sign)).times(k.plus(q.scale(sign))).minus(u.plus(v.scale(sign)).scale(sign)));
const one=scalar.one(),next=z.plus(one),cubic=x=>x.times(x).times(x).minus(x).scale('1/6');
prove(scalar,'central_sum_induction_increment',cubic(next).minus(cubic(z)),z.times(next).scale('1/2'));
const wt=['A','B','a','b'],wr=[];
for(let j=0;j<wt.length;j++)for(let k=0;k<j;k++){
 const lhs=[wt[j],wt[k]],rhs=[[[wt[k],wt[j]],1]];
 if((lhs[0]==='a'&&lhs[1]==='B')||(lhs[0]==='b'&&lhs[1]==='A'))rhs.push([[],1]);
 wr.push({id:'W_'+j+'_'+k,lhs,rhs,source:'Declared null oscillator pair at mode one: [a,B]=[b,A]=1, others zero'});
}
const weyl=new M.Presentation({tokens:wt,rules:wr});register('null_current_Weyl_pair',weyl);
const [A,B,a,b]=wt.map(t=>weyl.word([t])),comm=(x,y)=>x.times(y).minus(y.times(x));
prove(weyl,'cross_current_aB',comm(a,B),weyl.one());
prove(weyl,'cross_current_bA',comm(b,A),weyl.one());
prove(weyl,'null_current_aA',comm(a,A),weyl.zero());
prove(weyl,'double_contraction',comm(a.times(b),A.times(B)),A.times(b).plus(B.times(a)).plus(weyl.one()));
prove(weyl,'negative_metric_current',comm(a.minus(b),A.minus(B)),weyl.one().scale(-2));
prove(weyl,'positive_metric_current',comm(a.plus(b),A.plus(B)),weyl.one().scale(2));
const witness=proofs.find(x=>x.name==='cross_current_aB').certificate,bad=JSON.parse(JSON.stringify(witness));assert.ok(bad.steps.length);bad.steps[0].coefficient=['991','0'];
no('corrupted_native_word_witness',rejects(()=>weyl.replay(bad)));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'LQA1_LOCAL_QUANTUM_ANOMALY_V1',status:'PASS_COMPARATIVE_LOCAL_QUANTIZATION_AUDIT',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results:{original_central_charge_per_chirality:2,original_central_polynomial:'(m^3-m)/6',mixed_HD_central_polynomial:'(m^3-m)/3',contractions,mode_examples:modeExamples,Gram_examples:gramExamples,counterterm_solution:{r_times_s:solved[0][0],gamma:solved[1][0],source_A_coefficient:1,restored_hbar:'r*s=gamma=-hbar/12'},homogeneous_equation_original:'P*K=0',homogeneous_zero_mode_equation_corrected:'P*K-pi*hbar/3=0',reversed_vacuum:{central_charge_per_chirality:0,contractions:reverseContractions,kinematic_inner_product:'positive bosonic Fock',mode_grade:'unbounded in both directions; not identified with physical energy',homogeneous_zero_mode_equation:'P*K=0',strong_kernel_on_algebraic_test_domain:'zero; a physical extension and inner product remain required'},scope:{all_integer_algebra:'written proofs on declared common algebraic domains; finite native witnesses separately',first_vacuum_all_mode_strong_kernel:'zero in that representation',first_vacuum_kinematic_inner_product:'indefinite',counterterm_adopted_as_physical_theory:false,reversed_vacuum_adopted_as_next_working_candidate:true,positive_local_physical_space_constructed:false,homogeneous_HQS_changed:false,pointwise_coframe_exponential_constructed:false,quantum_gravity_complete:false,physical_FTL:false}},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,negative_controls:negative_controls.length,native_presentations:presentations.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),output=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(output,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(output,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha}));
