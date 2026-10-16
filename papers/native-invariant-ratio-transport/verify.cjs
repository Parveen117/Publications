'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),M=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const f=N.F.of,mm=(...a)=>a.reduce(N.mul),sc=N.scale,eq=N.equal,zero=N.isZero,cm=N.commutator;
const checks=[],negative_controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const tr=a=>a.reduce((s,row,i)=>s.add(row[i]),N.ZERO),real=z=>{assert.ok(z.turn.zero());return z.rad;},ex=(rho,A)=>real(tr(mm(rho,A)));
function unit(n,i,j){const a=N.zeros(n);a[i][j]=N.ONE;return a;}
const I=N.identity(2),X=N.matrix([[0,1],[1,0]]),Z=N.matrix([[1,0],[0,-1]]),R=N.matrix([[0,-1],[1,0]]),Y=sc(R,N.IOTA),A=sc(N.add(I,Z),'1/2'),B=sc(N.sub(I,Z),'1/2');
const rho=sc(N.add(I,N.add(sc(Z,'3/5'),sc(Y,'4/5'))),'1/2');
const dual=(h,a)=>sc(cm(h,a),N.IOTA),num=(h,a,b)=>ex(rho,dual(h,a)).mul(ex(rho,b)).sub(ex(rho,a).mul(ex(rho,dual(h,b))));
ck('ratio:positive_pure_source',eq(mm(rho,rho),rho)&&tr(rho).eq(1));
ck('ratio:initial',ex(rho,A).div(ex(rho,B)).eq(4));
ck('ratio:commuting_transport',num(Z,A,B).zero());
no('ratio:dimensionless_implies_conserved',!num(X,A,B).zero());
ck('ratio:noninvariant_derivative',num(X,A,B).div(ex(rho,B).pow(2)).eq(20));
ck('ratio:proportional_exception',num(X,sc(B,2),B).zero());
no('ratio:proportional_requires_separate_conservation',!zero(dual(X,B)));
const frame=N.matrix([['3/5','-4/5'],['4/5','3/5']]),rhof=mm(frame,rho,N.dagger(frame)),af=mm(frame,A,N.dagger(frame)),bf=mm(frame,B,N.dagger(frame));
ck('ratio:simultaneous_frame_covariance',ex(rhof,af).div(ex(rhof,bf)).eq(4));
const hb=[I,X,Y,Z],row=O=>N.matrix([hb.map(e=>tr(mm(O,e)))]),outer=(a,b)=>mm(N.dagger(a),b);
function polynomialDefect(h,a,b){const ar=row(a),br=row(b),ur=row(dual(h,a)),vr=row(dual(h,b));return N.sub(N.add(outer(ur,br),outer(br,ur)),N.add(outer(ar,vr),outer(vr,ar)));}
ck('ratio:all_state_polynomial_zero',zero(polynomialDefect(Z,A,B)));
no('ratio:all_state_polynomial_false',!zero(polynomialDefect(X,A,B)));
ck('ratio:trace_rate_obstruction',tr(mm(A,dual(X,A))).zero());

function rotor(radius){const n=2*radius+1,id=N.identity(n),u=N.zeros(n),l=N.zeros(n),s=N.zeros(n),t=N.zeros(n);
 for(let j=0;j<n;j++){l[j][j]=N.Cut.of(j-radius);s[n-1-j][j]=N.ONE;t[j][j]=N.Cut.of((j-radius)%2===0?1:-1);if(j+1<n)u[j+1][j]=N.ONE;}
 const ud=N.dagger(u),cos=sc(N.add(u,ud),'1/2'),l2=mm(l,l),hg=N.sub(N.add(l2,sc(id,8)),N.add(mm(u,u),mm(ud,ud))),v=N.sub(sc(id,'13/10'),sc(cos,'3/10'));
 return{radius,n,id,u,ud,l,l2,s,t,cos,hg,v,plus:sc(N.add(id,s),'1/2'),minus:sc(N.sub(id,s),'1/2')};}
const commutants=[];
for(const radius of [1,2,3]){
 const q=rotor(radius),n=q.n,basis=Array.from({length:n*n},(_,k)=>unit(n,Math.floor(k/n),k%n));
 const columns=basis.map(e=>[...N.flatten(cm(q.l2,e)),...N.flatten(cm(q.cos,e))]);
 const equations=N.matrix(Array.from({length:2*n*n},(_,i)=>columns.map(col=>col[i]))),rank=N.rank(equations);
 ck('commutant:dimension_two:'+radius,rank===n*n-2);
 ck('reflection:selfadjoint_involution:'+radius,eq(N.dagger(q.s),q.s)&&eq(mm(q.s,q.s),q.id));
 ck('reflection:shift_inverse:'+radius,eq(mm(q.s,q.u,q.s),q.ud));
 ck('reflection:momentum_sign:'+radius,eq(mm(q.s,q.l,q.s),sc(q.l,-1)));
 ck('reflection:geometry_generator:'+radius,zero(cm(q.s,q.hg)));
 ck('reflection:lapse:'+radius,zero(cm(q.s,q.v)));
 ck('reflection:projectors:'+radius,eq(mm(q.plus,q.plus),q.plus)&&eq(mm(q.minus,q.minus),q.minus)&&zero(mm(q.plus,q.minus)));
 ck('half_turn:uncoupled_symmetry:'+radius,zero(cm(q.t,q.hg)));
 no('half_turn:coupled_symmetry:'+radius,!zero(cm(q.t,q.v)));
 commutants.push({radius,carrier:n,equation_rank:rank,commutant_dimension:n*n-rank,reflection:q.s});
}
results.commutant_fixtures=commutants;
const tiny=rotor(1),larger=rotor(2);
ck('boundary:accidental_small_section_conservation',zero(cm(tiny.hg,tiny.cos)));
no('boundary:small_section_invariant_proves_infinite_invariant',!zero(cm(larger.hg,larger.cos)));
const q=rotor(3),h=N.sub(q.hg,sc(q.v,2)),ket=N.zeros(q.n,1);ket[q.radius][0]=N.ONE;ket[q.radius+1][0]=N.ONE;ket[q.radius-1][0]=N.ONE.neg();
const rho0=sc(mm(ket,N.dagger(ket)),'1/3'),u=N.cayley(sc(h,N.IOTA.neg()),'1/10'),rho1=mm(u,rho0,N.dagger(u));
ck('transport:unitary',eq(mm(N.dagger(u),u),q.id));
ck('transport:commutes_reflection',zero(cm(u,q.s)));
ck('transport:normalized_source',tr(rho0).eq(1)&&eq(mm(rho0,rho0),rho0));
for(const [label,state]of [['before',rho0],['after',rho1]]){
 ck('transport:even_weight:'+label,ex(state,q.plus).eq('1/3'));
 ck('transport:odd_weight:'+label,ex(state,q.minus).eq('2/3'));
 ck('transport:invariant_ratio:'+label,ex(state,q.plus).div(ex(state,q.minus)).eq('1/2'));
}
no('transport:invariant_ratio_implies_stationary_state',!eq(rho0,rho1));
no('transport:reflection_sector_is_scalar_lapse',!eq(mm(q.v,q.plus),sc(q.plus,'13/10')));
results.transport={scope:'Symmetric finite Fourier section and exact rational Cayley unitary; the infinite AGD conservation is proved in the theorem.',initial_density:rho0,H:h,U:u,weights:['1/3','2/3'],ratio:'1/2',state_changed:!eq(rho0,rho1)};

const even0=unit(q.n,q.radius,q.radius),k1=N.zeros(q.n,1);k1[q.radius-1][0]=N.ONE;k1[q.radius+1][0]=N.ONE;
const even1=sc(mm(k1,N.dagger(k1)),'1/2');
for(const [name,state,variance]of [['zero',even0,'9/200'],['cosine',even1,'27/400']]){
 ck('incomplete:even_sector:'+name,ex(state,q.plus).eq(1));
 ck('incomplete:same_lapse_mean:'+name,ex(state,q.v).eq('13/10'));
 ck('incomplete:different_variance:'+name,ex(state,mm(q.v,q.v)).sub(ex(state,q.v).pow(2)).eq(variance));
}
no('incomplete:invariant_record_is_complete_state',!eq(even0,even1));
results.record={even_weight_both:'1',mean_lapse_both:'13/10',lapse_variances:['9/200','27/400']};

// A full exact finite channel, including its erasure outcome.
const PL=N.matrix([[1,0,0],[0,1,0],[0,0,0]]),E=unit(3,2,2),A0=unit(3,0,0),Id3=N.identity(3);
function kraus(c,d){return[N.add(sc(PL,c),E),sc(unit(3,2,0),d),sc(unit(3,2,1),d)];}
const apply=(ks,rho)=>ks.reduce((s,k)=>N.add(s,mm(k,rho,N.dagger(k))),N.zeros(3));
const adjoint=(ks,a)=>ks.reduce((s,k)=>N.add(s,mm(N.dagger(k),a,k)),N.zeros(3));
function blockMap(rho,c){const out=N.zeros(3),qq=f(c).pow(2);for(let i=0;i<3;i++)for(let j=0;j<3;j++)out[i][j]=rho[i][j].mul(i<2&&j<2?qq:(i===2&&j===2?1:c));out[2][2]=out[2][2].add(rho[0][0].add(rho[1][1]).mul(f(1).sub(qq)));return out;}
const ks=kraus('3/5','4/5'),kr=kraus('5/13','12/13');
ck('channel:completeness',eq(adjoint(ks,Id3),Id3));ck('channel:second_completeness',eq(adjoint(kr,Id3),Id3));
for(let i=0;i<3;i++)for(let j=0;j<3;j++){
 const e=unit(3,i,j);
 ck('channel:block_identity:'+i+':'+j,eq(apply(ks,e),blockMap(e,'3/5')));
 ck('channel:semigroup:'+i+':'+j,eq(apply(ks,apply(kr,e)),blockMap(e,'3/13')));
 if(i<2&&j<2)ck('channel:common_dual_multiplier:'+i+':'+j,eq(adjoint(ks,e),sc(e,'9/25')));
}
const logicalKet=N.matrix([[2],[N.IOTA],[0]]),logical=sc(mm(logicalKet,N.dagger(logicalKet)),'1/5'),out=apply(ks,logical);
ck('channel:initial_quantum_coherence',!logical[0][1].zero());
ck('channel:logical_population',ex(out,A0).eq('36/125'));
ck('channel:survival',ex(out,PL).eq('9/25'));
ck('channel:conditional_ratio',ex(out,A0).div(ex(out,PL)).eq('4/5'));
ck('channel:conditional_full_state',eq(sc(mm(PL,out,PL),'25/9'),logical));
ck('channel:erasure_record',ex(out,E).eq('16/25'));
no('channel:invariant_ratio_retains_survival',!ex(out,PL).eq(ex(logical,PL)));
results.attenuation={q:'9/25',input:logical,output:out,population:'36/125',survival:'9/25',conditional_ratio:'4/5',erasure:'16/25'};

const n=f('36/125'),d=f('9/25'),rate=f('-2/3'),ndot=n.mul(rate),ddot=d.mul(rate);
ck('chart:equal_relative_rates',ndot.div(n).eq(ddot.div(d)));
ck('chart:ratio_derivative_zero',ndot.mul(d).sub(n.mul(ddot)).zero());
ck('chart:equal_log_scale_slope_one',ddot.div(d).div(ndot.div(n)).eq(1));
ck('chart:rescaled_log_slope_two',ddot.div(d).mul(2).div(ndot.div(n)).eq(2));
no('chart:dimensionless_forces_plotted_45_degrees',!f(2).eq(1));
no('chart:invariant_ratio_forces_raw_diagonal',!n.div(d).eq(1));
ck('chart:common_nonsemigroup_ratio',f('4/5').mul('1/3').div(f('1/3')).eq('4/5'));
no('chart:diagonal_forces_constant_exponential',!f('1/3').eq(f('1/2').pow(2)));
no('domain:zero_denominator',(()=>{try{f(1).div(0);return false;}catch{return true;}})());
results.chart={normalized_log_slope:'1',after_Y_scale_2:'2',raw_ratio:'4/5',nonsemigroup_q2:'1/3',nonsemigroup_q1_squared:'1/4'};

function register(name,spec){const s=new M.Presentation(spec),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w:Object.fromEntries(spec.tokens.map(t=>{const n=typeof t==='string'?t:t.name;return[n,s.word([n])];}))};}
function prove(name,lhs,rhs){const s=lhs.system,r=s.reduce(lhs.minus(rhs),{witness:true});ck(name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const source=JSON.parse(fs.readFileSync(path.join(__dirname,'../native-action-geometry-dynamics/CERTIFICATE.json'),'utf8')),spec=source.presentations.find(x=>x.name==='circle_momentum').spec;
 spec.tokens.push('s');spec.rules.push({id:'ss',lhs:['s','s'],rhs:[[[],1]],source:'Reflection squared is identity'},
 {id:'su',lhs:['s','u'],rhs:[[['v','s'],1]],source:'Reflection reverses circle shift'},
 {id:'sv',lhs:['s','v'],rhs:[[['u','s'],1]],source:'Reflection reverses inverse shift'},
 {id:'sh',lhs:['s','h'],rhs:[[['h','s'],1]],source:'Reflection leaves central action scale'},
 {id:'sp',lhs:['s','p'],rhs:[[['p','s'],-1]],source:'Reflection reverses angular momentum'});
 const{s:sys,w:{u,v,h,p,s}}=register('circle_with_reflection',spec),one=sys.one(),hg=p.times(p).plus(one.scale(8)).minus(u.times(u)).minus(v.times(v)),lapse=one.scale('13/10').minus(u.plus(v).scale('3/20')),plus=one.plus(s).scale('1/2'),minus=one.minus(s).scale('1/2');
 prove('reflection:kinetic_commutator',s.times(p).times(p).minus(p.times(p).times(s)),sys.zero());
 prove('reflection:geometry_commutator',s.times(hg).minus(hg.times(s)),sys.zero());
 prove('reflection:lapse_commutator',s.times(lapse).minus(lapse.times(s)),sys.zero());
 prove('reflection:positive_sector',plus.times(plus),plus);
 prove('reflection:negative_sector',minus.times(minus),minus);
 prove('reflection:orthogonal_sectors',plus.times(minus),sys.zero());
 prove('reflection:complete_sectors',plus.plus(minus),one);
 prove('reflection:transport_sector_generator',plus.times(hg).minus(hg.times(plus)),sys.zero());
}
{
 const tokens=['n','d','k'],rules=[];for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Scalar reading coefficients commute'});
 const{s,w:{n,d,k}}=register('common_scalar_rate',{tokens,rules});
 prove('ratio:common_rate_cancellation',k.times(n).times(d).minus(n.times(k).times(d)),s.zero());
}
const chosen=proofs.find(x=>x.certificate.steps.length),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['991','0'];
no('integrity:corrupted_word_proof',(()=>{try{systems.get(chosen.certificate.presentation).replay(bad);return false;}catch{return true;}})());
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),payload={protocol:'IRT1_DIMENSIONLESS_INVARIANT_RATIO_TRANSPORT_V1',status:'PASS_NATIVE_INVARIANT_RATIOS_REFLECTION_AND_CALIBRATED_DIAGONAL',scope:'Finite native certificates and written all-state ratio and infinite AGD common-commutant proofs. The 45-degree statement is restricted to a calibrated log chart; no physical quantum-gravity completion.',source_pins_sha256:hash(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls,results,summary:{exact_checks:checks.length,negative_controls:negative_controls.length,native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',sha=hash(bytes),file=path.join(__dirname,'CERTIFICATE.json');
if(process.argv.includes('--write')){fs.writeFileSync(file,bytes);fs.writeFileSync(path.join(__dirname,'EXPECTED.sha256'),sha+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(path.join(__dirname,'EXPECTED.sha256'),'utf8'),sha+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:sha,AGD_ratio:'1/2',commutant_dimension:2,attenuation_ratio:'4/5'}));
