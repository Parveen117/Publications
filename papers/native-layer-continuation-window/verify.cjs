'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const O=require('../native-prime-observation-closure/observation_adapter.cjs');
const W=require('../native-prime-innovation-window/window_adapter.cjs');
const L=require('../native-infinite-layer-cut/layer_adapter.cjs');
const S=require('../native-infinite-layer-cut/source_layers.cjs');
const C=require('./continuation_adapter.cjs');
const {F}=N,{Interval:I}=O,argv=process.argv.slice(2);
assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],negative_controls=[],fixtures=[],proofs=[],presentations=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const pos=x=>!F.of(x).le(0),overlap=(a,b)=>{a=I.of(a);b=I.of(b);return a.lo.le(b.hi)&&b.lo.le(a.hi);};
const rejects=f=>{try{f();return false;}catch{return true;}},mm=N.mul,adj=N.dagger,eq=N.equal;
const uniform=Array.from({length:12},(_,j)=>new F(1,j+1)),exact=[];
for(let n=1;n<=4;n++){
 const a=C.exactWindow(uniform,n,uniform[2*n]),p=uniform[2*n].sub(a.lower).div(a.width),theta=uniform[2*n+1].sub(a.oddLower).div(a.oddWidth),next=C.exactWindow(uniform,n+1,uniform[2*n+2]);
 ck('uniform:canonical_positions:'+n,p.eq(new F(n,2*n+1))&&theta.eq('1/2'));
 ck('uniform:odd_width_identity:'+n,a.oddWidth.eq(a.x.sub(a.lower).mul(a.upper.sub(a.x)).div(a.width)));
 ck('uniform:next_even_width:'+n,next.width.eq(a.width.mul(p).mul(new F(1).sub(p)).mul(theta).mul(new F(1).sub(theta))));
 ck('uniform:width_contraction:'+n,next.width.le(a.width.div(16)));
 const f=L.exactForms(uniform,n),Q=Array.from({length:2*n},(_,i)=>Array.from({length:2*n},(_,j)=>{
  if(j<n)return i%n===j?1:0;
  if(j<2*n-1){const k=j-n;return (i<n&&i===k?1:0)-(i%n===k+1?1:0);}
  return i===n-1?1:0;
 })),gram=mm(adj(Q),mm(L.block(f.A,N.zeros(n),N.zeros(n),f.B),Q));
 const expected=Array.from({length:2*n},(_,i)=>Array.from({length:2*n},(_,j)=>{
  if(i<n&&j<n)return uniform[i+j];
  if(i===2*n-1&&j===2*n-1)return uniform[2*n-1];
  if(i===2*n-1||j===2*n-1){const k=i===2*n-1?j:i;return k<n?uniform[n+k]:uniform[k].sub(uniform[k+1]);}
  if(i<n||j<n)return 0;
  return uniform[i+j-2*n+1].sub(uniform[i+j-2*n+2]);
 }));
 ck('uniform:general_Gram_congruence:'+n,eq(gram,N.matrix(expected)));
 const invsum=N.inverse(f.A)[n-1][n-1].rad.add(N.inverse(f.B)[n-1][n-1].rad);
 ck('uniform:reciprocal_window_width:'+n,invsum.mul(a.width).eq(1));
 const points=['-1/3','0','1/4','1/2','1','4/3'].map(F.of);
 ck('uniform:quadratic_identity_inside_and_outside:'+n,points.every(t=>{const z=C.exactWindow(uniform,n,a.lower.add(a.width.mul(t)));return z.oddWidth.eq(a.width.mul(t).mul(new F(1).sub(t)));}));
 const accepted=C.acceptance(uniform,n,uniform[2*n],uniform[2*n+1]);
 ck('uniform:outward_acceptance:'+n,accepted.accepted&&accepted.evenFraction.contains(p)&&accepted.oddFraction.contains(theta));
 const ext=L.forms(uniform,n+1);
 ck('uniform:Schur_pivots_match:'+n,overlap(accepted.oddLower,ext.factors.A.d[n])&&overlap(accepted.oddUpper,ext.factors.B.d[n]));
 const H=Array.from({length:n+1},(_,i)=>Array.from({length:n+1},(_,j)=>uniform[i+j])),g=Array.from({length:n+1},(_,j)=>[j===n?new F(1):N.mul(N.inverse(f.H),N.matrix(uniform.slice(n,2*n).map(z=>[z])))[j][0].rad.neg()]);
 ck('uniform:even_residual_norm:'+n,mm(adj(g),mm(H,g))[0][0].rad.eq(uniform[2*n].sub(a.lower)));
 exact.push({order:n,...a,p,theta,next_even_width:next.width});
}
fixtures.push({kind:'exact_uniform_moment_windows',moments:uniform,windows:exact,scope:'declared uniform measure on [0,1], not the actual xi source'});
const basic=[1,'1/2'],x=new F('3/10'),window=C.exactWindow(basic,1,x);
ck('one_layer:explicit_next_window',window.lower.eq('1/4')&&window.upper.eq('1/2')&&window.oddLower.eq('9/50')&&window.oddUpper.eq('11/50'));
for(const [name,xx,yy] of [['below_even','1/5','1/10'],['above_even','3/5','1/2'],['below_odd','3/10','17/100'],['above_odd','3/10','23/100'],['lower_endpoint','1/4','1/8'],['upper_endpoint','1/2','1/2']])no('next_layer:'+name,!C.acceptance(basic,1,xx,yy).accepted);
const bad=[1,'1/2','1/8'];
no('finite_binomial_triangle:implies_positive_moment_form',[[0,0],[0,1],[0,2],[1,0],[1,1],[2,0]].every(([r,k])=>L.exactDifference(bad,r,k).le(0)===false)&&!O.intervalLDL([[1,'1/2'],['1/2','1/8']]).positive);
no('positive_decreasing_layers:force_odd_window',!C.acceptance(basic,1,'3/10','17/100').accepted);
no('guard:missing_old_layer',rejects(()=>C.evenWindow([1],1)));
no('guard:nonunit_normalization',rejects(()=>C.evenWindow([2,1],1)));
no('guard:old_nonpositive_pair',rejects(()=>C.evenWindow([1,2],1)));
no('chosen_continuation:is_actual_source',C.acceptance(basic,1,'3/10','1/5').accepted&&!new F('3/10').eq('1/3'));
const lower=C.acceptance(basic,1,'3/10','17/100'),upper=C.acceptance(basic,1,'3/10','23/100');
no('one_positive_innovation:suffices',lower.oddUpper.positive()&&!lower.oddLower.positive()&&upper.oddLower.positive()&&!upper.oddUpper.positive());

const T=new P.Presentation({tokens:['p','t'],rules:[{id:'commuting_scalar_positions',lhs:['t','p'],rhs:[[['p','t'],1]],source:'real scalar canonical positions'}]});
const audit=T.audit();ck('words:scalar_window_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name:'scalar_continuation_windows',spec:T.spec(),audit});
const one=T.one(),p=T.word(['p']),t=T.word(['t']),pq=p.times(one.minus(p)),tq=t.times(one.minus(t)),ph=p.minus(one.scale('1/2')),th=t.minus(one.scale('1/2'));
function prove(name,left,right){const z=T.reduce(left.minus(right),{witness:true});ck('words:'+name,z.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(z.certificate));ck('words:replay:'+name,T.replay(certificate));proofs.push({name,presentation_hash:T.hash(),lhs:left.toJSON(),rhs:right.toJSON(),certificate});}
prove('factorized_odd_width',p.minus(p.times(p)),pq);
prove('one_step_quarter_bound_identity',one.scale('1/4').minus(pq),ph.times(ph));
prove('two_step_sixteenth_bound_identity',one.scale('1/16').minus(pq.times(tq)),ph.times(ph).scale('1/4').plus(pq.times(th.times(th))));
prove('complementary_innovations',pq.times(t).plus(pq.times(one.minus(t))),pq);

const parent=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-infinite-layer-cut/CERTIFICATE.json'),'utf8'));
const old=parent.fixtures.find(x=>x.kind==='complete_actual_single_point_layers');assert.ok(old);
const jet=S.sourceLayers(29,17),c=jet.moments,a=C.acceptance(c,8,c[16],c[17]),forms=L.forms(c,9),next=C.evenWindow(c,9);
ck('actual:complete_18_layers',jet.s===29&&jet.a.eq('3249/4')&&c.length===18&&jet.q_layers.length===18);
ck('actual:all_16_retained_layers_overlap',c.slice(0,16).every((z,k)=>overlap(z,old.jet.moments[k])));
ck('actual:complete_Hurwitz_and_Euler_tails',jet.components.hurwitz.length===17&&jet.components.hurwitz.every(z=>pos(z.remainder)&&z.shift===512&&z.terms===32)&&jet.components.euler.tails.length===18&&jet.components.euler.tails.every(pos)&&jet.components.euler.retained_prime_count===564);
ck('actual:source_window_accepted',a.accepted);
const inside=(z,lo,hi)=>pos(z.lo.sub(O.decimal(lo)))&&pos(O.decimal(hi).sub(z.hi));
ck('actual:canonical_position_bounds',inside(a.evenFraction,'0.170','0.171')&&inside(a.oddFraction,'0.145','0.146'));
ck('actual:window_width_bounds',inside(a.even.width,'1.224e-11','1.225e-11')&&inside(a.odd.width,'1.727e-12','1.728e-12')&&inside(next.width,'2.144e-13','2.145e-13'));
ck('actual:source_intervals_narrow',c.every(z=>z.hi.sub(z.lo).le(O.decimal('1e-66'))));
ck('actual:both_ninth_order_forms_positive',forms.factors.A.positive&&forms.factors.B.positive&&forms.factors.H.positive);
ck('actual:independent_LDL_innovations_overlap',overlap(a.oddLower,forms.factors.A.d[8])&&overlap(a.oddUpper,forms.factors.B.d[8]));
ck('actual:new_pivot_floors',pos(a.oddLower.lo.sub(O.decimal('2.509e-13')))&&pos(a.oddUpper.lo.sub(O.decimal('1.476e-12'))));
ck('actual:factorized_and_direct_width_overlap',overlap(a.factorizedWidth,a.odd.width));
const recurrence=a.even.width.mul(a.evenFraction).mul(new I(1).sub(a.evenFraction)).mul(a.oddFraction).mul(new I(1).sub(a.oddFraction));
ck('actual:next_even_width_recurrence_overlap',overlap(next.width,recurrence));
ck('actual:strict_sixteenth_contraction',pos(a.even.width.lo.div(16).sub(next.width.hi)));
const triangle=[];for(let r=0;r<=17;r++)for(let k=0;k+r<=17;k++)triangle.push({r,k,value:L.difference(c,r,k)});
ck('actual:171_binomial_readings_positive',triangle.length===171&&triangle.every(z=>z.value.positive()));
const lambda=F.of(old.held_out.lambda),fit=L.response(c,9,lambda).value.mul(jet.q0),error=I.of(old.held_out.actual_q).sub(fit),oldfit=I.of(old.held_out.layer_response),wheld=new F('61/2');
const analytic=wheld.add(jet.w).pow(2).div(wheld).mul(wheld.sub(jet.w).div(wheld.add(jet.w).sub(2)).pow(18));
ck('actual:unfitted_31_positive_improvement',fit.sub(oldfit).positive()&&error.positive()&&error.hi.le(I.of(old.held_out.error).lo));
ck('actual:independent_unconditional_error_bound',error.absUpper().le(analytic));
ck('actual:reported_error_bounds',error.hi.le(O.decimal('1.412e-31'))&&analytic.le(O.decimal('7.417e-25')));
no('actual:finite_execution_proves_RH',parent.scope.rh_proved===false&&parent.scope.all_order_actual_inequalities_proved===false);
no('actual:infinite_uniform_width_lower_bound',next.width.hi.le(a.even.width.lo.div(16)));
fixtures.push({kind:'actual_two_layer_continuation',jet,acceptance:a,forms,next_unexecuted_even_window:next,recurrence_width:recurrence,binomial_triangle:triangle,held_out:{s:31,lambda,actual_q:old.held_out.actual_q,layer_response:fit,error,analytic_bound:analytic},scope:'actual layers 0..17 and nine positive confluent directions; c18 is bounded conditionally for continuation but not evaluated'});
const payload={protocol:'LCW1_LAYER_CONTINUATION_WINDOW_V1',status:'PASS_LCW1_LAYER_CONTINUATION_WINDOW',source_parent:'1e62ad12a5a2e4c61643ae1f636ca694a3abd0a0',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,presentations,proofs,fixtures,scope:{sharp_finite_two_layer_windows_proved:true,exact_width_recurrence_proved:true,all_order_actual_source_window_membership_proved:false,rh_proved:false,source_derivatives:18,confluent_directions:9,cut_dilation_dimension:18,unexecuted_next_derivative:18,physical_mass_or_Maxwell_identification:false,global_eta_improved:false,arithmetic:'frozen native rational/rewrite engines and inherited outward intervals; Python 3.12 I/O only'},summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:negative_controls.length,actual_derivative_layers:18,actual_confluent_directions:9,positive_binomial_readings:triangle.length}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
