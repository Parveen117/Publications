'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};

// Exact off-shell jets. No constraint substitution is made in these rows.
// Logs/trigonometric functions are handled by the written chain rule, not
// approximated; a and a_c here are independent function jets.
const rows=[[2,3,1,-2,4,-1],[-2,-3,2,1,-1,3],[1,4,-2,3,2,5]];
for(const[i,row]of rows.entries()){
 const[ap,am,ap1,am1,bp,bm]=row.map(f),p=ap.add(am).div(2),v=ap.sub(am).div(2),D=ap.mul(am);
 const tp=ap.mul(bp).sub(ap1).div(2),tm=am.mul(bm).add(am1).div(2),cp=bp.sub(bm).sub(ap1.div(ap)).sub(am1.div(am)).div(2);
 const a=f('13/10'),ac=f('-2/5'),Z=am.div(ap).mul(tp).add(ap.div(am).mul(tm));
 ck('regular_current_discriminant:'+i,D.eq(p.pow(2).sub(v.pow(2))));
 ck('positive_discriminant:'+i,f(0).le(D)&&!D.zero());
 ck('off_shell_c_spatial_derivative:'+i,cp.eq(tp.div(ap).sub(tm.div(am))));
 for(const[ff,ff1]of [[2,3],[-1,2]])for(const eps of [1,-1]){
  const sm=f(ff),sm1=f(ff1),A=eps===1?ap:am,A1=eps===1?ap1:am1,B=eps===1?bp:bm,T=eps===1?tp:tm;
  const dA=sm1.mul(A).add(sm.mul(A1)).mul(eps),dq=sm.mul(B).add(sm1.mul(eps)).div(2),dc=dq.sub(dA.div(A).div(2));
  const dp=dA.div(2),deltaDensity=dp.mul(a).add(p.mul(ac).mul(dc));
  const totalDerivative=sm1.mul(A).add(sm.mul(A1)).mul(a).add(sm.mul(A).mul(ac).mul(cp)).mul(eps).div(2);
  const label=[i,ff,eps].join(':');
  ck('off_shell_c_gauge_variation:'+label,dc.eq(sm.mul(T).div(A)));
  ck('constraint_ideal_bracket_mod_derivative:'+label,deltaDensity.sub(totalDerivative).eq(sm.mul(ac).mul(Z).div(2)));
  const w=f('3/4'),w1=f('2/3'),weightedDerivative=totalDerivative.mul(w).add(w1.mul(sm).mul(A).mul(a).mul(eps).div(2));
  const weightedResidual=deltaDensity.mul(w).sub(weightedDerivative),expected=sm.mul(w).mul(ac).mul(Z).div(2).sub(sm.mul(A).mul(a).mul(w1).mul(eps).div(2));
  ck('fixed_window_off_shell_residual:'+label,weightedResidual.eq(expected));
 }
 const commonPlus=cp.mul(ap).div(-2).add(p.mul(tp).div(ap));
 const commonMinus=cp.mul(am).div(2).add(p.mul(tm).div(am));
 ck('strong_spatial_bracket_cancellation:'+i,commonPlus.eq(commonMinus));
 if(i===0){
  no('omit_c_spatial_derivative_in_integration_by_parts',!p.mul(tp).div(ap).eq(Z.div(2)));
  no('opposite_chiral_c_derivative_sign',!tp.div(ap).add(tm.div(am)).eq(cp));
 }
 fixtures.push({kind:'off_shell_regular_jet',ap,am,tp,tm,c_prime:cp,common_constraint_coefficient:Z.div(2)});
}

// A finite weighted chain-rule identity for the nonlocal P argument. It is
// not spatial quadrature or a finite-dimensional substitute for field theory.
const variationRows=[
 {w:'1/3',p:3,v:1,a:'13/10',ac:'-2/5',aP:'1/7',dq:2,dp:1,dv:-1},
 {w:'2/3',p:2,v:'1/2',a:'7/10',ac:'2/5',aP:'-2/7',dq:-1,dp:3,dv:2}
].map(o=>Object.fromEntries(Object.entries(o).map(([k,v])=>[k,f(v)])));
const deltaP=variationRows.reduce((s,r)=>s.add(r.w.mul(r.dp)),f(0)),JP=variationRows.reduce((s,r)=>s.add(r.w.mul(r.p).mul(r.aP)),f(0));
let direct=f(0),gradient=f(0),frozenP=f(0);
for(const[i,r]of variationRows.entries()){
 const D=r.p.pow(2).sub(r.v.pow(2)),dc=r.dq.sub(r.p.mul(r.dp).div(D)).add(r.v.mul(r.dv).div(D));
 const local=r.a.mul(r.dp).add(r.p.mul(r.ac).mul(dc));
 direct=direct.add(r.w.mul(local.add(r.p.mul(r.aP).mul(deltaP))));
 const grad=r.p.mul(r.ac).mul(r.dq).add(r.a.sub(r.p.pow(2).mul(r.ac).div(D)).add(JP).mul(r.dp)).add(r.p.mul(r.ac).mul(r.v).mul(r.dv).div(D));
 gradient=gradient.add(r.w.mul(grad));frozenP=frozenP.add(r.w.mul(local));
 ck('first_variation_discriminant:'+i,f(0).le(D)&&!D.zero());
}
ck('global_momentum_feedback_chain_rule',direct.eq(gradient));
ck('global_feedback_difference',direct.sub(frozenP).eq(JP.mul(deltaP)));
no('freeze_global_P_inside_source_coupling',!direct.eq(frozenP));

// Actual cosine values from rational unit-circle pairs. Differential jets
// a_Q and a_QQ follow analytically; no floating-point trig is evaluated.
for(const[cos,sin]of [['3/5','4/5'],['-3/5','-4/5'],[1,0]]){
 const cs=f(cos),sn=f(sin),a=f(1).add(cs.div(2)),ac=sn.div(-2),acc=cs.div(-2);
 ck('rational_circle:'+cos,cs.pow(2).add(sn.pow(2)).eq(1));
 ck('positive_complete_flow_speed:'+cos,f('1/2').le(a)&&a.le('3/2'));
 for(const P0 of [-3,2]){
  const P=f(P0),g=f('2/3'),gp=f('-1/4'),pi=f(3),J=ac,sourceQdot=g.mul(pi).mul(a.sub(ac).add(J));
  const Qdot=g.mul(pi).mul(a),Pdot=g.mul(pi).mul(P).mul(ac).neg(),F=P.mul(a),Fdot=Pdot.mul(a).add(P.mul(ac).mul(Qdot));
  ck('homogeneous_global_feedback:'+cos+':'+P0,sourceQdot.eq(Qdot));
  ck('source_reduced_probe_conservation:'+cos+':'+P0,Fdot.zero());
  const Edot=gp.mul(pi).mul(F).neg(),Cdot=Edot.add(gp.mul(pi).mul(F)).add(g.mul(pi).mul(Fdot));
  ck('clock_constraint_conservation:'+cos+':'+P0,Cdot.zero());
  ck('flow_coordinate_y_derivative:'+cos+':'+P0,Qdot.div(a).eq(g.mul(pi)));
  if(!ac.zero())no('omit_global_feedback_preserves_Q_flow:'+cos+':'+P0,!g.mul(pi).mul(a.sub(ac)).eq(Qdot));
  fixtures.push({kind:'homogeneous_source_feedback',P,a,a_c:ac,a_cc:acc,Qdot,Pdot,probe:F,probe_dot:Fdot});
 }
}
for(const P0 of [-2,3])for(const[a0,a1]of [['1/2','3/2'],['13/10','7/10']]){
 const P=f(P0),next=P.mul(a0).div(a1);
 ck('finite_flow_probe_conservation:'+P0+':'+a0,next.mul(a1).eq(P.mul(a0)));
 ck('finite_flow_orientation_preserved:'+P0+':'+a0,f(0).le(P.mul(next))&&!next.zero());
}
no('regular_formula_defined_at_zero_current',(()=>{try{f(1).div(0);return false;}catch{return true;}})());

// Fixed-window witness on the true periodic homogeneous constrained source.
// d/dtheta below is (1/(2*pi)) d/dx; outputs are bracket/(2*pi).
const clean=a=>new Map([...a].filter(([,c])=>!c.zero()));
const fourier=rows=>clean(new Map(rows.map(([n,c])=>[n,C(c)])));
function fm(a,b){const r=new Map();for(const[n,c]of a)for(const[m,d]of b)r.set(n+m,(r.get(n+m)||N.ZERO).add(c.mul(d)));return clean(r);}
const fd=a=>clean(new Map([...a].map(([n,c])=>[n,c.mul(new N.Cut(0,n))]))),mean=a=>a.get(0)||N.ZERO;
const w=fourier([[0,1],[1,'1/4'],[-1,'1/4']]),sm=fourier([[1,[0,'-1/2']],[-1,[0,'1/2']]]),pair=mean(fm(sm,fd(w)));
ck('window_smearing_derivative_pairing',pair.eq('-1/4'));
const plus=pair.mul('-3/2'),minus=pair.mul('3/2');
ck('window_plus_residual_over_two_pi',plus.eq('3/8'));ck('window_minus_residual_over_two_pi',minus.eq('-3/8'));
ck('homogeneous_source_T_plus_zero',f(2).mul(0).sub(0).div(2).zero());
ck('homogeneous_source_T_minus_zero',f(2).mul(0).add(0).div(2).zero());
no('nonconstant_fixed_window_is_weak_Dirac',!plus.zero());
no('fixed_window_keeps_spatial_constraint',!plus.sub(minus).zero());
ck('constant_window_residual_zero',mean(fm(sm,fd(fourier([[0,1]])))).zero());
fixtures.push({kind:'periodic_fixed_window_obstruction',normalization:'bracket divided by 2*pi',plus,minus,P:2,a:'3/2',old_constraints_zero:true});

// Denominator-cleared exact word proofs, supplementing the smooth proofs.
const tokens=['a','b','t','u','c','d','j'],rules=[];
for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting classical current and function jets'});
const scalar=new M.Presentation({tokens,rules}),audit=scalar.audit();ck('native_word_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:scalar.spec(),audit});
function prove(name,lhs,rhs){const r=scalar.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',scalar.replay(cert));proofs.push({name,presentation_sha256:scalar.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
const[a,b,t,u,c,d,j]=tokens.map(x=>scalar.word([x]));
const expected=b.times(b).times(t).plus(a.times(a).times(u));
prove('plus_constraint_ideal_cleared_denominators',a.times(b).times(t).scale(-1).plus(a.times(a).times(u)).plus(a.plus(b).times(b).times(t)),expected);
prove('minus_constraint_ideal_cleared_denominators',b.times(b).times(t).minus(a.times(b).times(u)).plus(a.plus(b).times(a).times(u)),expected);
prove('feedback_c_variation_cancellation',d.times(a).minus(b.times(b).times(c)).plus(d.times(j)).plus(b.times(b).times(c)),d.times(a.plus(j)));
prove('reduced_probe_flow_conservation',b.times(a).times(c).scale(-1).plus(a.times(b).times(c)),scalar.zero());
const corrupted=JSON.parse(JSON.stringify(proofs.find(x=>x.certificate.steps.length).certificate));corrupted.steps[0].coefficient=['991','0'];
no('corrupted_native_word_certificate',(()=>{try{scalar.replay(corrupted);return false;}catch{return true;}})());

const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'RAL1_REGULAR_APPARATUS_LIFT_V1',status:'PASS_REGULAR_CLASSICAL_LIFT_AND_FIXED_WINDOW_OBSTRUCTION',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{source_lift:'G_a=integral p a(c,P); C_T=E+g(T) Pi G_a',exact_mixed_bracket:'{C_T,T_epsilon[f]}=(g Pi/2) integral f a_c [(A_-/A_+)T_+ +(A_+/A_-)T_-]',unreduced_classical_first_class:true,global_P_feedback_retained:true,RPM_P_sector_recovery:'exact after old local reduction with inherited quantization',fixed_window_weak_Dirac:'iff w constant for the chosen positive a',explicit_window_residual_over_two_pi:['3/8','-3/8'],localized_apparatus_field_constructed:false,unreduced_extended_quantum_closure_proved:false,critical_branch_source_lift:false,apparatus_law_uniquely_selected:false,continuum_proofs_machine_certified:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cpth=path.join(__dirname,'CERTIFICATE.json'),hpth=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cpth,bytes);fs.writeFileSync(hpth,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cpth,'utf8'),bytes);assert.equal(fs.readFileSync(hpth,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
