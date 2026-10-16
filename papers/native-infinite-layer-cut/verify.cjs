'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const W=require('../native-prime-innovation-window/window_adapter.cjs');
const S=require('./source_layers.cjs'),L=require('./layer_adapter.cjs');
const {F}=N,{Interval:I}=O,argv=process.argv.slice(2);assert.ok(argv.length===1&&['--check','--write'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}},pos=x=>!F.of(x).le(0);
const overlap=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.hi)&&b.lo.le(a.hi);};
const real=M=>N.matrix(M).map(row=>row.map(x=>{assert.ok(x.turn.zero());return x.rad;}));
const det=M=>W.determinant(M),eq=N.equal,mm=N.mul,adj=N.dagger;
const dot=(a,b)=>a.reduce((z,x,i)=>z.add(F.of(x).mul(b[i])),new F(0));
const fall=(n,r)=>{let z=1n;for(let j=0;j<r;j++)z*=BigInt(n-j);return z;};

// An exact four-mode source exercises the universal algebra independently.
const atoms=['1/5','4/5','9/25','16/25'].map(F.of),weights=['1/8','2/8','3/8','2/8'].map(F.of),roots=['2/5','2/5','12/25','12/25'].map(F.of);
const moments=Array.from({length:18},(_,n)=>dot(weights,atoms.map(y=>y.pow(n))));
const exact=[];
for(let order=1;order<=4;order++){
 const f=L.exactForms(moments,order),boxed=L.forms(moments,order),id=N.identity(order),e=N.matrix(Array.from({length:order},(_,j)=>[j===0?1:0]));
 ck('toy:both_confluent_forms_positive:'+order,boxed.factors.A.positive&&boxed.factors.B.positive&&boxed.factors.H.positive);
 ck('toy:acquired_self_adjoint_contraction:'+order,eq(mm(adj(f.Y),f.H),f.A)&&eq(mm(f.H,N.sub(id,f.Y)),f.B));
 ck('toy:determinant_ratios:'+order,det(f.Y).eq(det(f.A).div(det(f.H)))&&det(N.sub(id,f.Y)).eq(det(f.B).div(det(f.H))));
 let v=e;const recovered=[];for(let k=0;k<2*order;k++){recovered.push(mm(adj(e),mm(f.H,v))[0][0].rad);v=mm(f.Y,v);}
 ck('toy:all_2N_coefficients_recovered:'+order,recovered.every((x,k)=>x.eq(moments[k])));
 exact.push({order,...f,recovered});
}
const toy=exact[3],V=N.matrix(atoms.map(y=>Array.from({length:4},(_,j)=>y.pow(j)))),diag=x=>N.matrix(x.map((z,i)=>x.map((_,j)=>i===j?z:0))),Sin=mm(N.inverse(V),mm(diag(roots),V)),wave=L.cutWave(toy.Y,Sin,toy.H),G=L.block(toy.H,N.zeros(4),N.zeros(4),toy.H),id8=N.identity(8);
ck('toy:source_square_root',eq(mm(Sin,Sin),mm(toy.Y,N.sub(N.identity(4),toy.Y)))&&eq(mm(adj(Sin),toy.H),mm(toy.H,Sin)));
ck('toy:positive_projection',eq(mm(wave.cut,wave.cut),wave.cut)&&eq(mm(adj(wave.cut),G),mm(G,wave.cut)));
ck('toy:two_reflections',eq(N.power(wave.reflection,2),id8)&&eq(N.power(wave.fixedReflection,2),id8));
ck('toy:wave_preserves_acquired_pairing',eq(mm(adj(wave.wave),mm(G,wave.wave)),G));
ck('toy:wave_determinant_one',det(wave.wave).eq(1));
const D=N.sub(N.scale(toy.Y,2),N.identity(4)),Db=L.block(D,N.zeros(4),N.zeros(4),D);
ck('toy:wave_recurrence',eq(N.add(N.sub(N.power(wave.wave,2),N.scale(mm(Db,wave.wave),2)),id8),N.zeros(8)));
const psi=N.matrix([[1],[0],[0],[0],[0],[0],[0],[0]]),observer=N.matrix([[1,2,0,0,0,0,0,0]]),inverse=N.inverse(wave.wave),initial=mm(observer,psi),energy=mm(adj(psi),mm(G,psi));
for(const depth of [1,3,8]){const state=mm(N.power(wave.wave,depth),psi),carried=mm(observer,N.power(inverse,depth));ck('toy:carried_observer_and_norm:'+depth,eq(mm(carried,state),initial)&&eq(mm(adj(state),mm(G,state)),energy));}
no('wave:fixed_observer_automatically_invariant',!eq(mm(observer,mm(wave.wave,psi)),initial));
no('layer_operator:is_already_a_projection',!eq(N.power(toy.Y,2),toy.Y));
no('wave:nonfixed_state_norm_convergence',!N.isZero(mm(N.sub(wave.wave,id8),psi)));
let avg=N.zeros(8);for(let n=0;n<5;n++)avg=N.add(avg,N.power(wave.wave,n));avg=N.scale(avg,'1/5');
ck('toy:finite_wave_average_bound',O.intervalLDL(real(N.sub(N.scale(G,'1/4'),mm(adj(avg),mm(G,avg))))).positive);
const shear=N.matrix([[1,1],[0,1]]);no('determinant_one:implies_norm_preservation',det(shear).eq(1)&&!eq(mm(adj(shear),shear),N.identity(2)));
no('positive_determinant:implies_positive_form',det([[-1,0],[0,-1]]).eq(1)&&!O.intervalLDL([[-1,0],[0,-1]]).positive);
fixtures.push({kind:'exact_four_mode_layer_cut',atoms,weights,moments,finite_banks:exact,square_root:Sin,cut:wave.cut,wave:wave.wave,pairing:G,scope:'declared rational four-mode source, not actual xi; dimension does not grow after four independent modes'});

const binary=[];
for(const depth of [1,5,15]){
 const ds=Array.from({length:depth+1},(_,j)=>L.exactDifference(moments,j,depth-j)),ws=ds.map((d,j)=>d.mul(L.choose(depth,j)));
 ck('toy:binomial_record_normalization:'+depth,ws.every(pos)&&ws.reduce((a,b)=>a.add(b),new F(0)).eq(1));
 ck('toy:binomial_record_recovers_falling_moments:'+depth,Array.from({length:depth+1},(_,r)=>ws.reduce((sum,w,j)=>sum.add(w.mul(new F(fall(j,r),fall(depth,r)))),new F(0)).eq(moments[r])).every(Boolean));
 binary.push({depth,branch_readings:ds,multiplicity_weights:ws});
}
fixtures.push({kind:'exact_combinatorial_recovery',records:binary});
const finite=L.forms(moments,4),changed=moments.slice();changed[8]=new F(0);changed[9]=new F(-1);
no('finite_layers:force_unobserved_next_signs',L.forms(changed,4).factors.A.positive&&L.forms(changed,4).factors.B.positive&&!L.forms(changed,5).factors.A.positive);
no('finite_four_modes:have_five_independent_directions',!L.forms(moments,5).factors.H.positive);
no('positive_decreasing_layers:suffice',[1,'1/2','1/8','1/64'].map(F.of).every(pos)&&!L.forms([1,'1/2','1/8','1/64'],2).factors.H.positive);
no('guard:incomplete_layer_pair',rejects(()=>L.forms([1,'1/2','1/3'],2)));
no('guard:nonunit_normalization',rejects(()=>L.forms([2,1],1)));
no('guard:missing_difference_layer',rejects(()=>L.difference([1,'1/2'],1,1)));
no('guard:outside_positive_source_domain',rejects(()=>L.response(moments,4,-1)));
no('guard:wrong_square_root',rejects(()=>L.cutWave(toy.Y,N.identity(4),toy.H)));
no('guard:uncarried_metric',rejects(()=>L.cutWave(toy.Y,Sin)));
no('guard:unsupported_even_source',rejects(()=>S.sourceLayers(28,15)));
no('guard:nondecreasing_tail_integrand',rejects(()=>S.eulerLayers(3,15,2)));

// Formal identities use declared functional-calculus relations, not examples.
const T=new P.Presentation({tokens:['Y','S','R'],rules:[
 {id:'root_square',lhs:['S','S'],rhs:[[['Y'],1],[['Y','Y'],-1]],source:'S is the positive square root of Y(I-Y)'},
 {id:'root_commutes',lhs:['S','Y'],rhs:[[['Y','S'],1]],source:'functional calculus of the same self-adjoint Y'},
 {id:'turn_Y',lhs:['R','Y'],rhs:[[['Y','R'],1]],source:'two-channel turn acts on a separate slot'},
 {id:'turn_S',lhs:['R','S'],rhs:[[['S','R'],1]],source:'two-channel turn acts on a separate slot'},
 {id:'turn_square',lhs:['R','R'],rhs:[[[],-1]],source:'retained native quarter-turn square'}
]});
const audit=T.audit();ck('words:confluent_layer_cut_algebra',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name:'layer_cut_wave',spec:T.spec(),audit});
const one=T.one(),y=T.word(['Y']),s=T.word(['S']),r=T.word(['R']),d=y.scale(2).minus(one),w=d.plus(s.times(r).scale(2));
function prove(name,left,right){const z=T.reduce(left.minus(right),{witness:true});ck('words:'+name,z.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(z.certificate));ck('words:replay:'+name,T.replay(certificate));proofs.push({name,presentation_hash:T.hash(),lhs:left.toJSON(),rhs:right.toJSON(),certificate});}
prove('cut_upper_diagonal',y.times(y).plus(s.times(s)),y);
prove('cut_off_diagonal',y.times(s).plus(s.times(one.minus(y))),s);
prove('cut_lower_diagonal',s.times(s).plus(one.minus(y).times(one.minus(y))),one.minus(y));
prove('circle_identity',d.times(d).plus(s.times(s).scale(4)),one);
prove('wave_unitarity',w.times(d.minus(s.times(r).scale(2))),one);
prove('wave_recurrence',w.times(w).minus(d.times(w).scale(2)).plus(one),T.zero());

// Actual higher layers: freshly evaluated, with every analytic tail retained.
const parent=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-prime-innovation-window/CERTIFICATE.json'),'utf8'));
const old=parent.fixtures.find(x=>x.kind==='complete_source_prime_update');assert.ok(old);
const jet=S.sourceLayers(29,15),forms=L.forms(jet.moments,8),floorA=O.decimal('1.51e-11'),floorB=O.decimal('6.38e-11');
ck('actual:anchor_and_complete_layer_count',jet.s===29&&jet.a.eq('3249/4')&&jet.q_layers.length===16&&jet.moments.length===16);
ck('actual:archived_value_and_first_derivative_overlap',overlap(jet.q0,old.new_jet.q)&&overlap(jet.q_layers[1].div(jet.a),old.new_jet.derivative));
ck('actual:archived_Q_and_Qprime_overlap',overlap(jet.Q_layers[0],old.new_jet.Q)&&overlap(jet.Q_layers[1],old.new_jet.Qprime));
ck('actual:Hurwitz_agrees_with_half_integer_trigamma',overlap(jet.components.hurwitz[0].value,old.new_jet.components.psi1));
ck('actual:archived_Euler_low_layers_overlap',overlap(jet.components.euler.coefficients[0],old.new_jet.components.euler.log_derivative)&&overlap(jet.components.euler.coefficients[1],old.new_jet.components.euler.log_derivative_prime));
ck('actual:all_Hurwitz_remainders_retained',jet.components.hurwitz.length===15&&jet.components.hurwitz.every(z=>pos(z.remainder)&&z.shift===512&&z.terms===32));
ck('actual:all_Euler_tails_retained',jet.components.euler.retained_prime_count===564&&jet.components.euler.last_retained_prime===4093&&jet.components.euler.tails.length===16&&jet.components.euler.tails.every(pos));
ck('actual:source_layers_narrow',jet.moments.every(z=>z.hi.sub(z.lo).le(O.decimal('1e-70'))));
ck('actual:normalization_for_unconditional_analytic_bound',jet.q0.positive()&&jet.q0.hi.le('1/2')&&jet.moments[0].lo.eq(1)&&jet.moments[0].hi.eq(1));
ck('actual:both_eighth_order_confluent_forms_positive',forms.factors.A.positive&&forms.factors.B.positive&&forms.factors.H.positive);
ck('actual:conservative_layer_pivot_floors',forms.factors.A.d.every(z=>pos(z.lo.sub(floorA)))&&forms.factors.B.d.every(z=>pos(z.lo.sub(floorB))));
const lower=O.intervalLDL(forms.A.map((row,i)=>row.map((x,j)=>x.sub(forms.H[i][j].mul('1/10000'))))),upper=O.intervalLDL(forms.H.map((row,i)=>row.map((x,j)=>x.mul('81/100').sub(forms.A[i][j]))));
ck('actual:strict_finite_spectral_interval',lower.positive&&upper.positive);
ck('actual:finite_wave_average_constant',pos(new F('25/4').mul('19/100').sub(1)));
const triangle=[];for(let r=0;r<=15;r++)for(let k=0;k+r<=15;k++)triangle.push({r,k,value:L.difference(jet.moments,r,k)});
ck('actual:all_136_binomial_readings_positive',triangle.length===136&&triangle.every(z=>z.value.positive()));
const binomialWeights=Array.from({length:16},(_,j)=>L.difference(jet.moments,j,15-j).mul(L.choose(15,j))),total=binomialWeights.reduce((s,w)=>s.add(w),new I(0));
ck('actual:depth_15_record_normalization_enclosed',binomialWeights.every(z=>z.positive())&&total.contains(1));
const held=old.held_out_jet,t=F.of(held.t),lambda=t.div(jet.a).sub(1),fit=L.response(jet.moments,8,lambda),qfit=fit.value.mul(jet.q0),error=I.of(held.q).sub(qfit),wheld=F.of(held.w),analytic=wheld.add(jet.w).pow(2).div(wheld).mul(wheld.sub(jet.w).div(wheld.add(jet.w).sub(2)).pow(16));
ck('actual:held_out_point_distinct_from_anchor',held.s===31&&lambda.eq('472/3249')&&!t.eq(jet.a));
ck('actual:held_out_positive_error_below_6e_minus_28',error.positive()&&pos(O.decimal('6e-28').sub(error.absUpper())));
ck('actual:independent_confluent_analytic_bound',analytic.eq('456261632/757403880336151185370351664061')&&error.absUpper().le(analytic)&&analytic.le(O.decimal('6.025e-22')));
no('actual:finite_layer_execution_proves_RH',parent.scope.rh_proved===false&&parent.scope.all_order_positivity_proved===false);
fixtures.push({kind:'complete_actual_single_point_layers',jet,forms,finite_spectral_lower_pivots:lower.d,finite_spectral_upper_pivots:upper.d,binomial_triangle:triangle,depth15_weights:binomialWeights,weight_sum:total,held_out:{s:31,t,lambda,actual_q:held.q,layer_response:qfit,error,analytic_bound:analytic},scope:'eight confluent directions at one actual anchor; not an added prime direction or all-order source certificate'});

const payload={protocol:'ILC1_INFINITE_LAYER_CUT_V1',status:'PASS_ILC1_INFINITE_LAYER_CUT',source_parent:'200c2697d14323b46b0087774cbcbdea0505bc85',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,
scope:{all_order_single_point_RH_criterion_proved:true,all_order_actual_inequalities_proved:false,rh_proved:false,finite_anchor_original_s:29,highest_actual_derivative:15,actual_confluent_dimension:8,selected_cut_dilation_dimension:16,unfitted_original_s:31,lambda:'dimensionless local displacement (t-a)/a; not identified with thermodynamic response lambda',wave:'unitary two-reflection transport; no Maxwell or physical frequency identification',infinite_depth:'carried readings invariant; nonfixed unitary states do not converge in norm',finite_uniform_mass_gap_claimed:false,global_eta_improved:false,physical_mass_or_quantum_gravity_closed:false,finite_witnesses_prove_universal_theorems:false,arithmetic:'frozen native rational engine and inherited outward intervals; Python 3.12 I/O only'},
summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length,actual_derivative_layers:16,actual_confluent_directions:8,positive_binomial_readings:triangle.length}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
