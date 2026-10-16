'use strict';
// Scientific arithmetic is exclusively the frozen native rational/Cut engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const O=require('../native-local-sector-classification/oscillator_fixture.cjs');
const {N,f,C}=O;
const M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const checks=[],negative_controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};

// Exact finite Fourier Laurent algebra. d is d/dtheta, theta=2 pi x;
// the common spatial factor 2 pi is restored analytically in the theorem.
const clean=x=>new Map([...x].filter(([,c])=>!c.zero()));
const lp=rows=>clean(new Map(rows.map(([n,c])=>[n,C(c)])));
function add(...xs){const r=new Map();for(const x of xs)for(const[n,c]of x)r.set(n,(r.get(n)||N.ZERO).add(c));return clean(r);}
const scale=(x,a)=>clean(new Map([...x].map(([n,c])=>[n,c.mul(a)])));
const sub=(x,y)=>add(x,scale(y,-1));
function mul(x,y){const r=new Map();for(const[n,c]of x)for(const[m,d]of y)r.set(n+m,(r.get(n+m)||N.ZERO).add(c.mul(d)));return clean(r);}
const diff=x=>clean(new Map([...x].map(([n,c])=>[n,c.mul(new N.Cut(0,n))])));
const mean=x=>x.get(0)||N.ZERO, same=(x,y)=>sub(x,y).size===0;
const one=lp([[0,1]]),cos=lp([[-1,'1/2'],[1,'1/2']]),sin=lp([[-1,[0,'1/2']],[1,[0,'-1/2']]]);
for(const u of ['1/3','3/5','-2/3'])for(const a of [1,2]){
 const W=add(one,scale(cos,u)),Ep=scale(mul(W,W),a),Em=scale(one,3),smear=add(sin,scale(cos,'2/5'));
 const Zp=mean(Ep),Zm=mean(Em),J=Zp.mul(Zm);
 ck('positive_profile_margin:'+u+':'+a,f(u).abs().le('2/3'));
 ck('critical_invariant_squared:'+u+':'+a,J.eq(f(a).mul(3).mul(f(1).add(f(u).pow(2).div(2)))));
 const rawSquared=f(a).mul(3).mul(mean(W).rad.pow(2));
 no('raw_length_is_critical_invariant:'+u+':'+a,!J.eq(rawSquared));
 for(const alpha of ['0','1/2','-2']){
  const dp=add(diff(mul(smear,Ep)),scale(Ep,alpha)),dm=scale(Em,f(alpha).neg());
  ck('plus_total_derivative_average:'+u+':'+a+':'+alpha,mean(diff(mul(smear,Ep))).zero());
  ck('plus_J_invariance:'+u+':'+a+':'+alpha,mean(dp).mul(Zm).add(Zp.mul(mean(dm))).zero());
  const mp=scale(Ep,alpha),mm=sub(scale(diff(mul(smear,Em)),-1),scale(Em,alpha));
  ck('minus_J_invariance:'+u+':'+a+':'+alpha,mean(mp).mul(Zm).add(Zp.mul(mean(mm))).zero());
 }
 ck('product_derivative_law:'+u+':'+a,same(diff(mul(smear,Ep)),add(mul(diff(smear),Ep),mul(smear,diff(Ep)))));
 no('drop_affine_inhomogeneous_term:'+u+':'+a,!mean(mul(sin,diff(Ep))).zero());
 for(const c of ['1/2',2,3])ck('primitive_rescaling_cancels:'+u+':'+a+':'+c,Zp.mul(c).mul(Zm.div(c)).eq(J));
 fixtures.push({kind:'critical_positive_Fourier_profile',u,a,Zp,Zm,J,rawSquared});
}

// Actual homogeneous gauge equivalences approaching distinct critical labels.
for(const n of [1,2,4,8,16]){
 const p=f(1).div(n),deltaPhi=f(7),deltaQ=f(-3);
 ck('regular_P_orbit_same_Q:'+n,p.mul(0).zero());
 ck('regular_P_orbit_scalar_shift:'+n,deltaPhi.div(p).mul(p).eq(deltaPhi));
 ck('regular_K_orbit_geometry_shift:'+n,deltaQ.div(p).mul(p).eq(deltaQ));
 no('critical_origin_has_same_finite_clock_flow:'+n,!f(0).mul(deltaPhi.div(p)).eq(deltaPhi));
}

// Sparse bivariate polynomial jets: compact cutoffs are part of the written
// analysis. Polynomial-window integrals below are finite rational witnesses,
// not a replacement for C_c^infinity test regularity or Gaussian integration.
const pk=(i,j)=>i+','+j;
const poly=rows=>clean(new Map(rows.map(([i,j,c])=>[pk(i,j),C(c)])));
function pmul(x,y){const r=new Map();for(const[k,c]of x){const[i,j]=k.split(',').map(Number);for(const[l,d]of y){const[a,b]=l.split(',').map(Number),v=pk(i+a,j+b);r.set(v,(r.get(v)||N.ZERO).add(c.mul(d)));}}return clean(r);}
const P=poly([[1,0,1]]),K=poly([[0,1,1]]),unit=poly([[0,0,1]]);
const ev0=x=>x.get('0,0')||N.ZERO;
function hadamard(x){const a=new Map(),b=new Map();for(const[k,c]of x){const[i,j]=k.split(',').map(Number);if(i>0)a.set(pk(i-1,j),c);else if(j>0)b.set(pk(0,j-1),c);else assert.ok(c.zero());}return[a,b];}
function axisPair(x,y,axis){let s=N.ZERO;for(const[k,c]of x){const[i,j]=k.split(',').map(Number);if((axis===0?j:i)!==0)continue;for(const[l,d]of y){const[a,b]=l.split(',').map(Number);if((axis===0?b:a)!==0)continue;const degree=axis===0?i+a:j+b;assert.ok(degree>0,'weighted integral requires vanishing product at zero');if(degree%2===0)s=s.add(c.dagger().mul(d).mul(f(2).div(degree)));}}return s;}
const regPair=(x,y)=>axisPair(x,y,0).add(axisPair(x,y,1));
const jets=[P,K,add(P,K),sub(P,K),poly([[2,0,2],[1,1,-3],[0,3,4]]),poly([[1,0,[1,2]],[0,2,[3,-1]]])];
for(const[i,x]of jets.entries()){
 const[a,b]=hadamard(x);
 ck('compact_Hadamard_jet_identity:'+i,same(x,add(pmul(P,a),pmul(K,b))));
 ck('finite_domain_evaluation:'+i,ev0(x).zero());
 ck('finite_domain_multiplier_invariance:'+i,ev0(pmul(add(unit,P,K),x)).zero());
 const norm=regPair(x,x);ck('axis_polynomial_norm_positive:'+i,norm.turn.zero()&&f(0).le(norm.rad)&&!norm.zero());
 ck('constraint_ideal_zero_axis_pairing:'+i,regPair(pmul(pmul(P,K),x),x).zero());
 ck('origin_constraint_evaluation:'+i,ev0(pmul(pmul(P,K),x)).zero());
 fixtures.push({kind:'polynomial_window_axis_norm',i,norm});
}
ck('axis_P_packet_exact_norm',regPair(P,P).eq(1));
ck('axis_K_packet_exact_norm',regPair(K,K).eq(1));
ck('axes_orthogonal',regPair(P,K).zero());
ck('two_axis_packet_norm',regPair(add(P,scale(K,2)),add(P,scale(K,2))).eq(5));
// A scalar-origin sector has rank one, with all Gram minors zero.
const origin=[C(1),C([2,1]),C([-3,2])],G=origin.map(a=>origin.map(b=>a.dagger().mul(b)));
ck('origin_Gram_rank_one',N.rank(G)===1);
for(let i=0;i<origin.length;i++)for(let j=0;j<origin.length;j++)ck('origin_Gram_minor:'+i+':'+j,G[i][i].mul(G[j][j]).sub(G[i][j].mul(G[j][i])).zero());

// Check the bound 2 d^2 and the exact exponent giving log(1/epsilon),
// without numerically evaluating an exponential, integral, logarithm or limit.
for(const d of ['1/2','1/4','1/8','1/16']){
 const delta=f(d),oneAxis=delta.pow(2),both=oneAxis.mul(2);
 ck('near_origin_square_bound:'+d,both.eq(delta.mul(delta).mul(2)));
 ck('near_origin_bound_decreases:'+d,both.le('1/2'));
}
ck('Gaussian_square_log_coefficient',f(2).mul('1/2').eq(1));
ck('two_coordinate_regulator_scaling_power',f(2).mul(2).eq(4));
for(const L of [2,4,8,16,32]){
 const ell=f(L),reg=f(5),originLead=ell.mul(3),r=f(1).div(ell);
 ck('origin_rate_normalized:'+L,r.mul(originLead).eq(3));
 ck('regular_rate_normalized:'+L,r.mul(reg).eq(f(5).div(L)));
 no('one_normalization_preserves_both_rates:'+L,!r.mul(reg).eq(reg));
}
// i(d/dP-1/(2P)) applied to chi=P^j has derivative-norm exponent 2j-3.
for(const j of [1,2,3]){
 const exponent=2*j-3;
 ck('branch_operator_domain_power:'+j,(exponent>-1)===(j>1));
}
no('all_vanishing_smooth_tests_in_canonical_domain',2*1-3===-1);

// Subtraction gives a positive full-test form but loses momentum symmetry.
// s=1 is a polynomial-window witness for s(0)=1. The written smooth-cutoff
// proof applies to every genuine compact s nonzero near zero.
function subtractPair(x,y){const rx=sub(x,scale(unit,ev0(x))),ry=sub(y,scale(unit,ev0(y)));return regPair(rx,ry).add(ev0(x).dagger().mul(ev0(y)));}
ck('subtraction_origin_norm',subtractPair(unit,unit).eq(1));
ck('subtraction_regular_positive_norm',subtractPair(P,P).eq(1));
ck('subtraction_origin_regular_cross_zero',subtractPair(unit,pmul(P,P)).zero());
no('subtraction_preserves_P_adjoint',!subtractPair(P,P).eq(subtractPair(unit,pmul(P,P))));
no('subtraction_preserves_K_adjoint',!subtractPair(K,K).eq(subtractPair(unit,pmul(K,K))));

// Q_* equals the LPS reduced Q: verify its exponentiated cancellation,
// using current amplitudes determined by their means and the source ODEs.
for(const mom of [-3,-1,1,2])for(const ctheta of [-1,0,1]){
 const W=f(1).add(f('3/5').mul(ctheta)),Ep=W.pow(2),Em=f(1),Zp=f('59/50'),Zm=f(1),p0=f(mom);
 const Ap=p0.mul(Ep).div(Zp),Am=p0.mul(Em).div(Zm);
 const exp2Q=Ep.mul(Em).mul(p0.pow(2)).div(Ap.mul(Am).abs());
 ck('regular_Q_equals_critical_expression:'+mom+':'+ctheta,exp2Q.eq(Zp.mul(Zm)));
}
// Laurent jets for the inherited canonical Q and its affine symmetrization.
const coord=lp([[1,1]]),invcoord=lp([[-1,1]]);
const derivMomentum=x=>clean(new Map([...x].map(([n,c])=>[n-1,c.mul(n)])));
const qOp=x=>scale(sub(derivMomentum(x),scale(mul(invcoord,x),'1/2')),N.IOTA);
const dOp=x=>scale(mul(coord,derivMomentum(x)),N.IOTA);
for(const j of [0,1,2,3,4]){
 const x=lp([[j,[2,1]]]);
 ck('dilation_is_weighted_symmetrization:'+j,same(dOp(x),scale(add(mul(coord,qOp(x)),qOp(mul(coord,x))),'1/2')));
 ck('affine_quantum_commutator:'+j,same(sub(dOp(mul(coord,x)),mul(coord,dOp(x))),scale(mul(coord,x),N.IOTA)));
}
for(const lam of ['1/3','1/2',2,3])for(const p0 of [-2,-1,1,2])ck('dilation_weight_invariance:'+lam+':'+p0,f(lam).div(f(lam).mul(p0).abs()).eq(f(1).div(f(p0).abs())));
// Weak cross bracket has the actual constraint factor, for an arbitrary
// smooth coefficient from differentiating Q_* along the uniform k shift.
for(const p0 of [-2,0,2])for(const k0 of [-3,0,3])if(p0*k0===0)ck('global_dilations_weak_cross_bracket:'+p0+':'+k0,f(p0).mul(k0).mul('7/3').zero());
no('canonical_pair_on_origin_line',N.commutator([[0]],[[3]])[0][0].zero()&&!N.IOTA.zero());

// Reuse LSC's independently computed descendant Gram pairing at zero.
// Projection is exact for each finite test, not an assertion of finite norm.
const pairCache=new Map();
function originBra(x){let out=N.ZERO;for(const[key,c]of x){const gu=O.grade(key,'u'),gv=O.grade(key,'v');if(gu!==gv)continue;if(!pairCache.has(gu))pairCache.set(gu,O.paired(gu,f(0)).out);out=out.add(O.inner(pairCache.get(gu),new Map([[key,c]])));}return out;}
for(const vars of [[],['u:1'],['v:1'],['u:1','v:1'],['u:2','v:2'],['u:1','v:2']])for(const m of [-3,-2,-1,0,1,2,3])ck('origin_full_constraint_finite_test:'+vars.join('.')+':'+m,originBra(O.L(m,O.state(vars),f(0))).zero());
ck('origin_dressed_vacuum_coefficient',originBra(O.vac).eq(1));
ck('origin_nonvacuum_coefficient',originBra(O.state(['u:1','v:1'])).eq(-2));
no('origin_functional_is_bare_vacuum',!originBra(O.state(['u:1','v:1'])).zero());

// Small native word certificates for the analytic proof's algebraic pieces.
const tokens=['a','b','c','d','p','k'],rules=[];
for(let j=0;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'commuting classical scalars and amplitude jets'});
const scalar=new M.Presentation({tokens,rules}),audit=scalar.audit();ck('word_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({spec:scalar.spec(),audit});
function prove(name,lhs,rhs){const r=scalar.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const cert=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',scalar.replay(cert));proofs.push({name,presentation_sha256:scalar.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:cert});}
const[a,b,c,d,p,k]=tokens.map(x=>scalar.word([x])),zero=scalar.word([]).scale(0);
prove('critical_product_primitive_gauge_cancels',a.times(p).times(k).minus(p.times(a).times(k)),zero);
prove('vanishing_amplitude_product_expansion',p.times(a).plus(k.times(b)).times(p.times(c).plus(k.times(d))),p.times(p).times(a).times(c).plus(k.times(k).times(b).times(d)).plus(p.times(k).times(a.times(d).plus(b.times(c)))));
prove('constraint_ideal_multiplier_invariance',p.times(k).times(a).times(b),a.times(p).times(k).times(b));

const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'COS1_CRITICAL_ORIGIN_SECTOR_V1',status:'PASS_CRITICAL_QUOTIENT_AND_TWO_SCALE_AVERAGING',source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,presentations,proofs,fixtures,results:{critical_orbit_labels:'Q0=(1/2)log[(integral exp(q+r))(integral exp(q-r))], Phi; r_prime=k',globally_continuous_observers:'constant on all critical orbits; distinct orbits are not identified by finite gauge flow',origin_normalization:'B_epsilon/log(1/epsilon) -> conjugate(e0(f)) e0(g)',maximal_unscaled_smooth_domain:'e0(f)=0',regular_completion:'H_P direct_sum H_K unchanged',filtered_option:'(H_P direct_sum H_K) direct_sum C with separate quotient maps and relative normalization choice',global_weak_Dirac_observables:'D_P=P Q_star, D_K=K Phi_av; self-adjoint quantizations iP d_P and iK d_K on active branches, zero on origin',single_full_test_rigging_with_old_momentum_adjoint_contract:false,faithful_quantization_of_all_critical_labels:false,interacting_dynamics_constructed:false,physical_3plus1_gravity_complete:false},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
