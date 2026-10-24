'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs'),S=require('./full_source.cjs'),W=require('./window_adapter.cjs');
const {F}=N,{Interval:I}=O,argv=process.argv.slice(2);assert.ok(argv.length===1&&['--check','--write'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],presentations=[],proofs=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});},no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}},pos=x=>!F.of(x).le(0),width=x=>x.hi.sub(x.lo),overlap=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.hi)&&b.lo.le(a.hi);};
const inside=(x,lo,hi)=>pos(x.lo.sub(lo))&&pos(F.of(hi).sub(x.hi)),mm=N.mul,eq=N.equal;
function register(tokens,name){const rules=[];for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'declared central real source scalars'});const T=new P.Presentation({tokens,rules}),audit=T.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:T.spec(),audit});return T;}
function prove(T,name,lhs,rhs){const r=T.reduce(lhs.minus(rhs),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,T.replay(certificate));proofs.push({name,presentation_hash:T.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}

// Exact toy source, separately labelled from the actual completed xi source.
const atoms=[2,5,11,19].map(F.of),weights=[1,2,3,2].map(F.of);
const q=t=>atoms.reduce((sum,m,i)=>sum.add(weights[i].div(F.of(t).add(m))),new F(0)),dq=t=>atoms.reduce((sum,m,i)=>sum.sub(weights[i].div(F.of(t).add(m).pow(2))),new F(0));
const nodes=['25/4','81/4','169/4'].map(F.of),bank=W.exactPencil(nodes,nodes.map(q),nodes.map(dq));
function rawSigns(c,r,d){const ell=nodes.map((t,i)=>[bank.values[i].sub(r).div(c.x.sub(t))]),j=nodes.map((t,i)=>[c.x.mul(r).sub(t.mul(bank.values[i])).div(c.x.sub(t))]),quad=(M,z)=>mm(N.dagger(z),mm(N.inverse(M),z))[0][0].rad;return {alpha:d.neg().sub(quad(bank.L,ell)),beta:r.add(c.x.mul(d)).sub(quad(bank.J,j)),ell,j};}
for(const x of [1,10,'3249/4']){
 const c=W.exactCoefficients(bank,x),r=q(x),d=dq(x),delta=r.sub(c.r0),deltaPrime=d.sub(c.d0),alpha=deltaPrime.neg().add(c.b.mul(delta).mul(2)).sub(c.A.mul(delta.pow(2))),beta=delta.add(c.x.mul(deltaPrime)).sub(c.x.mul(c.b).mul(delta).mul(2)).sub(c.x.pow(2).mul(c.C).mul(delta.pow(2))),raw=rawSigns(c,r,d),gate=W.acceptance(bank,x,r,d),tag=String(x);
 ck('toy:residual_cross_columns:'+tag,eq(raw.ell,N.sub(mm(bank.L,c.v),N.scale(c.a,delta)))&&eq(raw.j,N.add(mm(bank.J,c.v),N.scale(c.a,c.x.mul(delta)))));
 ck('toy:both_residual_signs_equal_original:'+tag,alpha.eq(raw.alpha)&&beta.eq(raw.beta)&&pos(alpha)&&pos(beta));
 ck('toy:factored_budget:'+tag,beta.add(c.x.mul(alpha)).eq(delta.mul(new F(1).sub(c.kappa.mul(delta)))));
 ck('toy:determinant_log_derivative:'+tag,c.b.eq(c.a.reduce((s,row)=>s.add(row[0].rad),new F(0)).sub(c.traceInverseL)));
 ck('toy:outward_window_contains_exact_signs:'+tag,gate.accepted&&gate.alpha.contains(alpha)&&gate.beta.contains(beta)&&inside(gate.valueFraction,0,1)&&inside(gate.theta,0,1));
 const R=delta.div(c.g),Rp=deltaPrime.sub(c.b.mul(delta).mul(2)).div(c.g),lower=R.div(c.x).neg().add(c.x.mul(c.C).mul(c.g).mul(R.pow(2))),upper=c.A.mul(c.g).mul(R.pow(2)).neg();
 ck('toy:normalized_residual_slope_cone:'+tag,pos(R)&&pos(Rp.sub(lower))&&pos(upper.sub(Rp))&&Rp.le(0));
 fixtures.push({kind:'exact_toy_residual_window',x:c.x,r,d,r0:c.r0,d0:c.d0,A:c.A,C:c.C,b:c.b,kappa:c.kappa,h:c.h,g:c.g,delta,deltaPrime,alpha,beta,R,Rp,scope:'four declared rational atoms; not actual xi'});
}
{
 const c=W.exactCoefficients(bank,'3249/4'),unit=new F(1).div(c.kappa),mid=unit.div(2),low=c.b.mul(mid).mul(2).sub(mid.div(c.x)).add(c.x.mul(c.C).mul(mid.pow(2))),high=c.b.mul(mid).mul(2).sub(c.A.mul(mid.pow(2))),span=high.sub(low),r=c.r0.add(mid);
 ck('window:exact_width_at_interior_value',span.eq(mid.mul(new F(1).sub(c.kappa.mul(mid))).div(c.x))&&pos(span));
 for(const [name,d] of [['above',c.d0.add(high).add(span)],['below',c.d0.add(low).sub(span)]]){const signs=rawSigns(c,r,d);no('window:interior_value_with_bad_slope:'+name,name==='above'?!pos(signs.alpha):!pos(signs.beta));}
 for(const [name,delta] of [['negative',mid.neg()],['too_large',unit.add(mid)]])no('window:outside_value_can_be_repaired_by_slope:'+name,!pos(delta.mul(new F(1).sub(c.kappa.mul(delta)))));
 for(const [name,delta] of [['lower',new F(0)],['upper',unit]]){const dp=c.b.mul(delta).mul(2).sub(c.A.mul(delta.pow(2))),signs=rawSigns(c,c.r0.add(delta),c.d0.add(dp));no('window:strict_signs_at_endpoint:'+name,signs.alpha.zero()&&signs.beta.zero());}
 fixtures.push({kind:'sharp_value_and_slope_controls',x:c.x,lowerValue:c.r0,upperValue:c.valueUpper,interiorValue:r,lowerDerivative:c.d0.add(low),upperDerivative:c.d0.add(high),derivativeWidth:span});
}
{
 const c=W.exactCoefficients(bank,10),B=N.matrix([[1,2,0],[0,1,1],[0,0,2]]),Bt=N.dagger(B),L=mm(Bt,mm(bank.L,B)),J=mm(Bt,mm(bank.J,B)),u=mm(Bt,bank.values.map(z=>[z])),a=mm(Bt,c.a),Pencil=N.add(J,N.scale(L,10)),v=mm(N.inverse(Pencil),u),quad=(M,z)=>mm(N.dagger(z),mm(N.inverse(M),z))[0][0].rad;
 ck('covariance:all_retained_scalar_contractions',quad(Pencil,u).eq(c.r0)&&quad(L,a).eq(c.A)&&quad(J,a).eq(c.C)&&mm(N.dagger(a),v)[0][0].rad.eq(c.b)&&mm(N.dagger(v),mm(L,v))[0][0].rad.neg().eq(c.d0));
 ck('covariance:determinant_ratio',W.determinant(L).div(W.determinant(Pencil)).eq(W.determinant(bank.L).div(W.determinant(c.p))));
 no('covariance:uncarried_new_probe',!quad(L,c.a).eq(c.A));
 ck('normalization:old_double_zeros_have_nonzero_denominator',nodes.every(t=>pos(W.determinant(N.add(bank.J,N.scale(bank.L,t))))));
 fixtures.push({kind:'old_basis_congruence_with_carried_probe',B,L,J,u,a,scope:'sampling labels retained; the probe column moves with the old source basis'});
}
no('guard:duplicate_next_node',rejects(()=>W.coefficients(bank,nodes[0])));
no('guard:nonpositive_next_node',rejects(()=>W.coefficients(bank,0)));
no('guard:indefinite_old_form',rejects(()=>W.coefficients({...bank,L:bank.L.map(row=>row.map(x=>x.neg()))},10)));
no('guard:duplicate_old_nodes',rejects(()=>W.exactPencil([1,1],[1,1],[-1,-1])));

// Constants, complete special-function source and the new observation.
const constants=S.constants(),bernoulli=S.bernoulli(50),t2=new F(5,12),t4=t2.mul(2).div(new F(1).sub(t2.pow(2))),machinTangent=t4.sub('1/239').div(new F(1).add(t4.div(239)));
ck('constants:Machin_rational_tangent_and_range',machinTangent.eq(1)&&inside(constants.pi,'314/100','22/7'));
ck('constants:Bernoulli_reference_values',bernoulli[0].eq(1)&&bernoulli[1].eq('-1/2')&&bernoulli[2].eq('1/6')&&bernoulli[4].eq('-1/30')&&bernoulli[6].eq('1/42'));
ck('constants:logarithm_inverse_and_product_controls',S.logRational('1/2').add(S.logRational(2)).contains(0)&&S.logRational(3).add(S.logRational(5)).sub(S.logRational(15)).contains(0));
no('guard:nonpositive_logarithm',rejects(()=>S.logRational(0)));
no('guard:unsupported_even_source_argument',rejects(()=>S.sourceJet(28)));
no('guard:floating_source_argument',rejects(()=>S.sourceJet(29.5)));
{
 const p=4093,ps=new F(BigInt(p)**29n),den=ps.sub(1),log=S.logRational(p),stable=log.square().mul(ps.div(den.pow(2))),split=log.square().mul(ps).div(den.pow(2));
 ck('arithmetic:exact_ratio_retains_complete_weight',split.lo.le(stable.lo)&&stable.hi.le(split.hi)&&pos(stable.lo)&&width(stable).le(O.decimal('1e-100')));
 no('arithmetic:split_huge_weight_gives_narrow_enclosure',pos(width(split).sub('1/100000')));
 fixtures.push({kind:'fixed_grid_precision_control',p,s:29,exact_weight:ps.div(den.pow(2)),stable,split,scope:'both formulas enclose the same quantity; only the unsplit rational weight has useful precision'});
}
const archive=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-prime-observation-closure/sources/SE49_RESULT.json'),'utf8'));
const old=O.primeBank(archive),cross=S.sourceJet(23),newJet=S.sourceJet(29),held=S.sourceJet(31);
ck('source:independent_old_value_and_derivative_agree',overlap(cross.q,archive.values[20])&&overlap(cross.derivative,I.of(archive.L[20][20]).neg()));
ck('source:new_complete_jets_have_narrow_enclosures',[newJet,held].every(j=>width(j.q).le(O.decimal('1e-98'))&&width(j.derivative).le(O.decimal('1e-98'))));
ck('source:Euler_cutoff_all_primes_and_nonzero_tails',[cross,newJet,held].every(j=>j.components.euler.cutoff===4096&&j.components.euler.retained_prime_count===564&&j.components.euler.last_retained_prime===4093&&pos(j.components.euler.tail0)&&pos(j.components.euler.tail1)));
ck('source:positive_real_source_signs',[cross,newJet,held].every(j=>j.q.positive()&&j.derivative.hi.le(0)&&j.components.euler.log_derivative.hi.le(0)&&j.components.euler.log_derivative_prime.positive()));
ck('source:Gamma_remainder_is_explicit_and_retained',[cross,newJet,held].every(j=>j.components.psi.remainder.lo.le(j.components.psi.first_omitted)&&j.components.psi.first_omitted.le(j.components.psi.remainder.hi)&&!j.components.psi.first_omitted.zero()));
const acceptance=W.acceptance(old,newJet.t,newJet.q,newJet.derivative),enlarged=W.append(old,newJet),L=O.intervalLDL(enlarged.L),J=O.intervalLDL(enlarged.J);
ck('actual:prime_29_new_mixed_direction',newJet.s===29&&newJet.t.eq('3249/4')&&!old.nodes.some(t=>new F(t).eq(newJet.t))&&enlarged.nodes.length===9&&enlarged.L.every(row=>row.length===9)&&enlarged.J.every(row=>row.length===9));
ck('actual:both_ninth_order_forms_positive',L.positive&&J.positive&&L.d.length===9&&J.d.length===9);
ck('actual:direct_and_factored_innovations_agree',acceptance.accepted&&overlap(acceptance.alpha,acceptance.rawAlpha)&&overlap(acceptance.beta,acceptance.rawBeta)&&overlap(acceptance.rawAlpha,L.d[8])&&overlap(acceptance.rawBeta,J.d[8]));
ck('actual:conservative_new_innovation_floors',pos(acceptance.rawAlpha.lo.sub(O.decimal('1.36e-18')))&&pos(acceptance.rawBeta.lo.sub(O.decimal('6.68e-15'))));
ck('actual:sharp_value_and_slope_coordinates',inside(acceptance.valueFraction,'140/1000','141/1000')&&inside(acceptance.theta,'142/1000','143/1000'));
ck('actual:normalized_source_slope_in_its_cone',acceptance.normalizedResidual.positive()&&pos(acceptance.normalizedDerivative.lo.sub(acceptance.normalizedLower.hi))&&pos(acceptance.normalizedUpper.lo.sub(acceptance.normalizedDerivative.hi))&&acceptance.normalizedDerivative.hi.le(0));
const fit=O.response(enlarged.L,enlarged.J,enlarged.values,held.t),error=held.q.sub(fit.value),roots=[...old.primes,29].map(p=>new F(2*p-1,2)),analytic=held.w.add('5/2').pow(2).div(held.w).mul(O.schwarzProduct(roots,held.w,1)),oldFit=O.response(old.L,old.J,old.values,held.t),oldError=held.q.sub(oldFit.value);
ck('actual:held_out_31_not_in_fit',held.s===31&&held.t.eq('3721/4')&&!enlarged.nodes.some(t=>new F(t).eq(held.t)));
ck('actual:held_out_error_below_3e_minus_17',error.positive()&&pos(O.decimal('3e-17').sub(error.absUpper())));
ck('actual:held_out_error_inside_independent_analytic_bound',error.absUpper().le(analytic));
ck('actual:new_fit_improves_same_held_out_point',fit.value.lo.le(held.q.hi)&&oldFit.value.hi.le(fit.value.lo)&&error.absUpper().le(oldError.lo));
no('scope:old_finite_archive_proves_RH',archive.rh_proved===false&&archive.all_order_positivity_certified===false&&archive.global_shift_frontier_unchanged==='1079/2000');
fixtures.push({kind:'complete_source_prime_update',constants,cross_check_jet:cross,new_jet:newJet,held_out_jet:held,primes:[...old.primes,29],bank:enlarged,acceptance,L_pivots:L.d,J_pivots:J.d,held_out:{s:31,t:held.t,actual_q:held.q,prime_bank_q:fit.value,actual_minus_finite:error,absolute_error_upper:error.absUpper(),finite_analytic_bound:analytic,old_eight_prime_error:oldError},scope:'new ninth prime direction beyond the archived s<=26 samples; not the full old V24 plus one direction'});

{
 const T=register(['x','d','e','b','A','C'],'residual_innovation_scalars'),[x,d,e,b,a,c]=['x','d','e','b','A','C'].map(t=>T.word([t])),one=T.one(),d2=d.times(d),alpha=e.scale(-1).plus(b.times(d).scale(2)).minus(a.times(d2)),beta=d.plus(x.times(e)).minus(x.times(b).times(d).scale(2)).minus(x.times(x).times(c).times(d2)),kappa=x.times(a.plus(x.times(c)));
 prove(T,'two_sign_budget_factorization',x.times(alpha).plus(beta),d.times(one.minus(kappa.times(d))));
 prove(T,'slope_window_width_after_x_clear',x.times(b.times(d).scale(2).minus(a.times(d2))).minus(x.times(b.times(d).scale(2)).minus(d).plus(x.times(x).times(c).times(d2))),d.times(one.minus(kappa.times(d))));
}
{
 const T=register(['g','R','Rp','b','A','C','x'],'normalized_source_residual'),[g,r,rp,b,a,c,x]=['g','R','Rp','b','A','C','x'].map(t=>T.word([t])),delta=g.times(r),dp=g.times(rp.plus(b.times(r).scale(2)));
 prove(T,'normalized_alpha',dp.scale(-1).plus(b.times(delta).scale(2)).minus(a.times(delta).times(delta)),g.times(rp).scale(-1).minus(a.times(g).times(g).times(r).times(r)));
 prove(T,'normalized_beta',delta.plus(x.times(dp)).minus(x.times(b).times(delta).scale(2)).minus(x.times(x).times(c).times(delta).times(delta)),g.times(r).plus(x.times(g).times(rp)).minus(x.times(x).times(c).times(g).times(g).times(r).times(r)));
}
const payload={protocol:'PIW1_PRIME_INNOVATION_WINDOW_V1',status:'PASS_PIW1_PRIME_INNOVATION_WINDOW',source_parent:'ab7c44cabf73b23898038deb0f649aa1d5de3669',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,
scope:{value_window:'sharp existence interval for a derivative admitting both actual source forms',slope:'actual derivative is independently required; interior value alone is insufficient',covariance:'old-source congruence with u and a co-transported',source:'new complete xi jets in the Euler half-plane; explicit Gamma, pole, prime and all-tail bounds',old_bank:'eight unchanged prime directions selected from complete SE49 archive',new_accepted_s:29,held_out_s:31,full_old_V24_plus_new_direction_certified:false,all_order_positivity_proved:false,rh_proved:false,global_eta_improved:false,physical_mass_or_quantum_gravity_closed:false,finite_witnesses_prove_universal_theorems:false,arithmetic:'native rational F and 110-place outward intervals; Python 3.12 I/O only'},
summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length,actual_prime_observations:9,source_jets_evaluated:3}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
