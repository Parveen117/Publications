'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const A=require('./observation_adapter.cjs'),O=require('../native-prime-observer-sectors/observer_adapter.cjs');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const {F,Cut}=N,{Interval:IV}=A,argv=process.argv.slice(2);
assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],presentations=[],proofs=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}};
const eq=N.equal,mm=N.mul,ff=F.of,positive=x=>!ff(x).le(0);
const same=(x,y)=>JSON.stringify(x)===JSON.stringify(y),overlap=(x,y)=>{x=IV.of(x);y=IV.of(y);return x.lo.le(y.hi)&&y.lo.le(x.hi);};
function register(tokens,name){const rules=[];for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'declared real scalar coordinate algebra'});const S=new P.Presentation({tokens,rules}),audit=S.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});return S;}
function prove(S,name,lhs,rhs){const r=S.reduce(lhs.minus(rhs),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,S.replay(certificate));proofs.push({name,presentation_hash:S.hash(),lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}

// The metric and lattice are the already declared TET/POS adapted frame.
const points=[[2,0],[0,3],[5,0]],c=A.centre(points),H0=N.matrix([[2,1],[1,5]]);
ck('centre:prime_content_inputs',same(points.map(x=>O.content(x).toString()),['2','3','5']));
ck('centre:actual_source_system',eq(c.H,H0)&&eq(c.A,[[-2,26],[12,6]])&&eq(c.b,[[37],[42]])&&c.determinant.eq(-324));
ck('centre:exact_solution',eq(c.x,[['145/54'],['44/27']])&&eq(mm(c.A,c.x),c.b));
ck('centre:three_equal_positive_readings',c.radii.every(x=>x.eq(c.radii[0]))&&positive(c.radii[0].rad));
no('centre:automatic_prime_content_closure',!O.integral(c.x));
const contrast=N.matrix([[-1,1,0],[-1,0,1]]),mean=[new Cut('1/3'),new Cut('1/3'),new Cut('1/3')],readings=c.radii.map(x=>[x]),complete=contrast.concat([mean]);
ck('reading:contrast_kernel_exactly_common_line',N.rank(contrast)===2&&N.isZero(mm(contrast,[[1],[1],[1]])));
ck('reading:mean_and_contrasts_recover_all',eq(mm(N.inverse(complete),mm(complete,readings)),readings));
no('reading:zero_contrast_is_zero_reading',N.isZero(mm(contrast,readings))&&!N.isZero(readings));
{
 const B=N.matrix([[1,1],[0,1]]),Bi=N.inverse(B),H=mm(N.dagger(Bi),mm(H0,Bi)),moved=c.points.map(x=>mm(B,x).map(r=>r[0])),transported=A.centre(moved,H),fixed=A.centre(moved,H0),endpoint=O.rationalEndpoint(),end=A.centre(c.points.map(x=>mm(endpoint.U,x).map(r=>r[0])),endpoint.H);
 ck('centre:common_frame_covariance',eq(transported.x,mm(B,c.x))&&transported.radii.every((x,i)=>x.eq(c.radii[i])));
 ck('centre:actual_TET_endpoint_covariance',eq(endpoint.H0,H0)&&eq(end.x,mm(endpoint.U,c.x))&&end.radii.every((x,i)=>x.eq(c.radii[i])));
 no('centre:fixed_metric_under_nonisometric_shear',!eq(fixed.x,mm(B,c.x)));
 fixtures.push({kind:'source_metric_centre_and_transport',source:c,shear:B,transported,actual_TET_endpoint:{parameter:endpoint.t,U:endpoint.U,H:endpoint.H,centre:end.x},fixed_metric_centre:fixed.x});
}
no('guard:collinear_centre',rejects(()=>A.centre([[2,0],[3,0],[5,0]])));
no('guard:nonpositive_metric',rejects(()=>A.centre(points,[[1,0],[0,-1]])));
no('guard:nonsymmetric_metric',rejects(()=>A.centre(points,[[1,1],[0,1]])));
no('guard:complex_reading_point',rejects(()=>A.centre([[new Cut(2,1),0],[0,3],[5,0]])));
{
 const modular=Object.fromEntries([2,3,5,7,9,25,35,49].map(n=>[n,A.residues(c.A,c.b,n)]));
 ck('local:good_moduli_unique', [5,7,25,35,49].every(n=>same(modular[n],[A.uniqueResidue(c.A,c.b,n)])));
 ck('local:prime_depths_compatible',[5,7].every(p=>same(A.uniqueResidue(c.A,c.b,p*p).map(x=>(BigInt(x)%BigInt(p)).toString()),A.uniqueResidue(c.A,c.b,p))));
 ck('local:explicit_good_prime_and_CRT_values',same(modular[5],[['0','2']])&&same(modular[7],[['1','5']])&&same(modular[35],[['15','12']]));
 ck('local:inherited_CRT_map_agrees',same(O.crtInverse(modular[5][0],modular[7][0],5,7).map(String),modular[35][0]));
 ck('local:bad_prime_three_has_three_first_layer_solutions',same(modular[3],[['0','2'],['1','0'],['2','1']]));
 no('local:every_prime_has_a_centre',modular[2].length===0);
 no('local:first_layer_forces_lift',modular[3].length===3&&modular[9].length===0);
 ck('local:mod_nine_obstruction_identity',eq(N.add(N.scale([c.A[0]],6),[c.A[1]]),[[0,162]])&&c.b[0][0].mul(6).add(c.b[1][0]).eq(264)&&264n%9n!==0n);
 no('local:nonunit_determinant_unique_inverse',rejects(()=>A.uniqueResidue(c.A,c.b,3)));
 no('local:nonunit_rational_denominator',rejects(()=>A.rationalMod('1/3',9)));
 fixtures.push({kind:'exact_modular_centres',determinant:c.determinant,modular,scope:'finite witnesses of the written adjugate and inverse-limit theorem'});
}

// All interval endpoints are native exact fractions, rounded outwards only.
ck('interval:exact_scientific_decimal',A.decimal('-1.25e-7').eq('-1/8000000')&&A.decimal('2.50E2').eq(250));
ck('interval:positive_and_negative_directed_rounding',[new F(1,3),new F(-1,3),new F(1,7)].every(x=>new IV(x).contains(x)));
ck('interval:arithmetic_contains_exact_results',new IV('1/3').add('1/7').contains('10/21')&&new IV('-1/3').mul('2/7').contains('-2/21')&&new IV('1/3').div('-2/7').contains('-7/6'));
ck('interval:square_across_zero',same(new IV(-2,1).square().toJSON(),['0','4']));
no('interval:tiny_negative_rounded_positive',!new IV('-1e-150').positive());
no('guard:division_through_zero',rejects(()=>new IV(1).div([-1,1])));
no('guard:reversed_interval',rejects(()=>new IV(2,1)));
no('guard:floating_input',rejects(()=>A.decimal(0.1)));
no('interval:indefinite_form_certified',!A.intervalLDL([[1,2],[2,1]]).positive);
no('guard:asymmetric_interval_form',rejects(()=>A.intervalLDL([[1,0],[1,1]])));
no('guard:response_dimension_mismatch',rejects(()=>A.response([[1]],[[1]],[1,2],1)));
no('guard:response_nonpositive_point',rejects(()=>A.response([[1]],[[1]],[1],0)));

// Independent exact rational Stieltjes fixture checks the inherited pencil
// identities. It is explicitly not the actual completed xi response.
const toyMu=[2,5,11].map(ff),toyWeights=[1,2,3].map(ff),nodes=['25/4','81/4','169/4'].map(ff);
const toyQ=t=>toyMu.reduce((s,m,i)=>s.add(toyWeights[i].div(ff(t).add(m))),new F(0));
const toyD=t=>toyMu.reduce((s,m,i)=>s.sub(toyWeights[i].div(ff(t).add(m).pow(2))),new F(0));
function pencil(ts,q,d){const u=ts.map(q),L=ts.map((t,i)=>ts.map((x,j)=>i===j?d(t).neg():u[i].sub(u[j]).div(x.sub(t)))),J=ts.map((t,i)=>ts.map((x,j)=>i===j?u[i].add(t.mul(d(t))):x.mul(u[j]).sub(t.mul(u[i])).div(x.sub(t))));return {u,L,J};}
function exactResponse(bank,z){const col=bank.u.map(x=>[x]),inv=N.inverse(N.add(bank.J,N.scale(bank.L,z))),v=mm(inv,col);return {q:mm(N.dagger(col),v)[0][0].rad,d:mm(N.dagger(v),mm(bank.L,v))[0][0].rad.neg()};}
const toy=pencil(nodes,toyQ,toyD),T=nodes.map((x,i)=>nodes.map((_,j)=>i===j?x:0)),ones=N.matrix([[1],[1],[1]]),uCol=toy.u.map(x=>[x]);
ck('toy:both_pencils_positive',A.intervalLDL(toy.L).positive&&A.intervalLDL(toy.J).positive);
ck('toy:both_displacement_identities',eq(N.add(toy.J,mm(toy.L,T)),mm(uCol,N.dagger(ones)))&&eq(N.add(toy.J,mm(T,toy.L)),mm(ones,N.dagger(uCol))));
ck('toy:values_and_derivatives_interpolate',nodes.every(t=>{const r=exactResponse(toy,t);return r.q.eq(toyQ(t))&&r.d.eq(toyD(t));}));
ck('toy:full_three_pole_response_recovered',exactResponse(toy,'2809/4').q.eq(toyQ('2809/4'))&&A.response(toy.L,toy.J,toy.u,'2809/4').value.contains(toyQ('2809/4')));
{
 const old=pencil(nodes.slice(0,2),toyQ,toyD),ell=toy.L.slice(0,2).map(r=>[r[2]]),j=toy.J.slice(0,2).map(r=>[r[2]]),quad=(M,v)=>mm(N.dagger(v),mm(N.inverse(M),v))[0][0].rad;
 const alpha=toy.L[2][2].sub(quad(old.L,ell)),beta=toy.J[2][2].sub(quad(old.J,j));
 ck('toy:both_actual_update_innovations_positive',positive(alpha)&&positive(beta));
 const z=new F('2809/4'),Pc=N.add(old.J,N.scale(old.L,z)),k=N.add(j,N.scale(ell,z)),delta=toy.J[2][2].add(z.mul(toy.L[2][2])).sub(quad(Pc,k)),e=toy.u[2].sub(mm(N.dagger(k),mm(N.inverse(Pc),old.u.map(x=>[x])))[0][0].rad),increment=exactResponse(toy,z).q.sub(exactResponse(old,z).q);
 ck('toy:response_increment_and_lower_bound',increment.eq(e.pow(2).div(delta))&&beta.add(z.mul(alpha)).le(delta));
 const bad=toy.L.map(r=>r.slice());bad[2][2]=new F(0);
 no('toy:new_prime_forces_positive_innovation',!A.intervalLDL(bad).positive);
 fixtures.push({kind:'explicit_toy_Stieltjes_pencil',poles:toyMu.map(x=>x.neg()),weights:toyWeights,nodes,...toy,alpha,beta,held_out_response:exactResponse(toy,z),increment,scope:'algebra fixture only; not actual xi or a zeta-zero catalogue'});
}
{
 const eps=new F('1/100000000'),maskedQ=t=>{const a=t.add(4);return toyQ(t).add(eps.mul(2).mul(a).div(a.pow(2).add(1)));},maskedD=t=>{const a=t.add(4);return toyD(t).add(eps.mul(2).mul(new F(1).sub(a.pow(2))).div(a.pow(2).add(1).pow(2)));},masked=pencil(nodes,maskedQ,maskedD),l=A.intervalLDL(masked.L),j=A.intervalLDL(masked.J);
 ck('masking:finite_positive_data_with_off_axis_poles',l.positive&&j.positive);
 no('masking:finite_positive_bank_forces_all_poles_on_negative_axis',positive(eps)&&l.positive&&j.positive);
 fixtures.push({kind:'toy_finite_sign_masking_countermodel',epsilon:eps,additional_poles:[new Cut(-4,1),new Cut(-4,-1)],residue_at_each_added_pole:eps,nodes,L:masked.L,J:masked.J,L_pivots:l.d,J_pivots:j.d,scope:'a rational countermodel to the finite-to-global inference; not a zero of actual xi'});
}

// Source provenance is separately bound by SHA-256 pins and Git blob hashes.
const source=JSON.parse(fs.readFileSync(path.join(__dirname,'sources/SE49_RESULT.json'),'utf8')),bank=A.primeBank(source),l=A.intervalLDL(bank.L),j=A.intervalLDL(bank.J),ts=bank.nodes.map(ff),us=bank.values.map(IV.of);
ck('source:correct_prime_subbank',same(bank.primes,[3,5,7,11,13,17,19,23])&&same(bank.indices,[0,2,4,8,10,14,16,20]));
ck('source:both_outward_LDL_certificates',l.positive&&j.positive&&l.d.length===8&&j.d.length===8);
ck('source:conservative_pivot_floors',positive(l.minimumLower.sub('1/10000000000000000000'))&&positive(j.minimumLower.sub('1/1000000000000000')));
ck('source:common_positive_anchor_below_half',us[0].positive()&&us[0].hi.le('1/2'));
ck('source:all_retained_displacement_identities_enclosed',ts.every((t,i)=>ts.every((x,k)=>overlap(IV.of(bank.J[i][k]).add(IV.of(bank.L[i][k]).mul(x)),us[i])&&overlap(IV.of(bank.J[i][k]).add(IV.of(bank.L[i][k]).mul(t)),us[k]))));
ck('source:all_off_diagonal_divided_differences_consistent',ts.every((t,i)=>ts.every((x,k)=>i===k||(overlap(us[i].sub(us[k]).div(x.sub(t)),bank.L[i][k])&&overlap(us[k].mul(x).sub(us[i].mul(t)).div(x.sub(t)),bank.J[i][k])))));
no('source:finite_archive_claims_RH',source.rh_proved===false&&source.all_order_positivity_certified===false&&source.global_eta_zero_certified===false);
ck('source:global_shift_frontier_retained',source.global_shift_frontier_unchanged==='1079/2000'&&source.global_shift_floor==='1/492000');
const held=source.held_out_response_readings[0],w=new F(held.held_out_w),z=w.pow(2),actual=IV.of(held.actual_q),fit=A.response(bank.L,bank.J,bank.values,z),error=actual.sub(fit.value),ws=bank.primes.map(p=>new F(2*p-1,2)),B=A.schwarzProduct(ws,w,1),analytic=w.add('5/2').pow(2).div(w).mul(B);
ck('source:held_out_point_not_in_fit',held.held_out_s===27&&z.eq('2809/4')&&!ts.some(t=>t.eq(z)));
ck('source:positive_held_out_error_below_1e_minus_15',error.positive()&&positive(new F('1/1000000000000000').sub(error.absUpper())));
ck('source:exact_finite_analytic_bound',analytic.eq('103938208563200/38778237372401894277'));
ck('source:independent_error_inside_analytic_bound',error.absUpper().le(analytic));
fixtures.push({kind:'actual_SE49_prime_bank_replay',primes:bank.primes,source_indices_zero_based:bank.indices,nodes:bank.nodes,L_pivots:l.d,J_pivots:j.d,minimum_L_pivot_lower:l.minimumLower,minimum_J_pivot_lower:j.minimumLower,held_out:{s:27,w,t:z,actual_q:actual,prime_bank_q:fit.value,actual_minus_prime_bank:error,absolute_error_upper:error.absUpper(),B_N_1:B,finite_analytic_bound:analytic},scope:'principal selection of the existing complete V24 boxes, not an enlarged source space or fresh tail recomputation'});

// Analytic fixtures witness factors and obstructions; the infinite-rate and
// meromorphic-identity arguments remain the written theorems, not finite tests.
for(const a of [0,1]){
 const test=new Cut(2,'3/2'),products=ws.map((_,i)=>A.schwarzProduct(ws.slice(0,i+1),test,a));
 ck('density:strict_factor_contraction:'+a,products.every((x,i)=>positive(x)&&x.le(i?products[i-1]:1)&&!x.eq(i?products[i-1]:1)));
}
ck('density:exact_zero_at_a_retained_probe',A.schwarzProduct(ws,'9/2',0).zero());
no('density:uniform_contraction_on_boundary',ws.every(b=>new Cut(0,'3/2').sub(b).norm2().eq(new Cut(0,'3/2').add(b).norm2())));
const exponential=Array.from({length:12},(_,i)=>new F(2n**BigInt(i+1))),expProduct=A.schwarzProduct(exponential,'1/8');
ck('density:geometric_deficit_inequalities',exponential.every((b,i)=>new F(1).sub(new F(1,2n**BigInt(i+2))).le(A.schwarzProduct([b],'1/8'))));
no('density:observation_count_alone_forces_decay',new F('1/2').le(expProduct));
no('density:duplicate_value_and_derivative_add_independent_rank',N.rank([[1,2],[1,2]])===1);
no('guard:Schwarz_point_outside_domain',rejects(()=>A.schwarzProduct(ws,1,1)));
fixtures.push({kind:'analytic_density_controls',geometric_nodes:exponential,test_w:'1/8',finite_product:expProduct,universal_geometric_lower_bound:'1/2',boundary_test:new Cut(0,'3/2'),scope:'the prime-reciprocal proof and compact-uniform estimates are in POC-T4'});

{
 const S=register(['x','y','a','b','c','d','h','k','l'],'quadratic_reading_coordinates'),[x,y,a,b,c0,d,h,k,l0]=['x','y','a','b','c','d','h','k','l'].map(t=>S.word([t]));
 const Q=(u,v)=>h.times(u).times(u).plus(k.times(u).times(v).scale(2)).plus(l0.times(v).times(v));
 const lhs=Q(x.minus(c0),y.minus(d)).minus(Q(x.minus(a),y.minus(b))),rhs=a.minus(c0).times(h.times(x).plus(k.times(y))).scale(2).plus(b.minus(d).times(k.times(x).plus(l0.times(y))).scale(2)).plus(Q(c0,d)).minus(Q(a,b));
 prove(S,'quadratic_difference_cancels_common_source',lhs,rhs);
}
{
 const S=register(['r','b','a','t'],'half_plane_factor'),[r,b,a,t]=['r','b','a','t'].map(x=>S.word([x])),sq=x=>x.times(x);
 prove(S,'Schwarz_denominator_difference',sq(r.plus(b).minus(a.scale(2))).plus(sq(t)).minus(sq(r.minus(b)).plus(sq(t))),r.minus(a).times(b.minus(a)).scale(4));
}
{
 const S=register(['A','h','x','y','d'],'scalar_Schur_completion'),[a,h,x,y,d]=['A','h','x','y','d'].map(t=>S.word([t])),sq=x=>x.times(x);
 prove(S,'Schur_congruence_polynomial',a.times(sq(x.plus(h.times(y)))).plus(d.minus(a.times(sq(h))).times(sq(y))),a.times(sq(x)).plus(a.times(h).times(x).times(y).scale(2)).plus(d.times(sq(y))));
}
{
 const S=register(['x','t','u','r'],'divided_difference_interpolation'),[x,t,u,r]=['x','t','u','r'].map(q=>S.word([q]));
 prove(S,'cross_pencil_interpolation_after_denominator_clear',x.times(r).minus(t.times(u)).plus(x.times(u.minus(r))),x.minus(t).times(u));
}

const payload={protocol:'POC1_PRIME_OBSERVATION_CLOSURE_V1',status:'PASS_POC1_PRIME_OBSERVATION_CLOSURE',source_parent:'4b4a1c4c3b2d51b6f202a4b88c41b05f606f66e9',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,
 scope:{centre:'chosen squared-distance readings in the retained TET metric; not a canonical embedding of primes',local:'good primes give compatible unique centres; bad-prime lifting must be tested',actual_RH_source:'complete pinned SE49 enclosures, sampled at eight odd-prime original s values',prime_density:'written Euler-product and Schwarz-factor proof; analytic uniqueness density, not topological density',all_order_criterion:'RH iff both actual source pencils are positive at every odd-prime bank order',all_order_positivity_proved:false,rh_proved:false,new_source_directions_beyond_V24:false,original_special_function_tails_recomputed:false,global_eta_improved:false,finite_model_poles_identified_as_actual_zeros:false,physical_mass_or_spacetime_quantization_proved:false,numerical_checks_prove_written_analytic_theorems:false,rounding:'110 decimal places, directed native rational bounds; no binary floating scientific arithmetic'},
 summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length,actual_prime_observations:bank.primes.length}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
