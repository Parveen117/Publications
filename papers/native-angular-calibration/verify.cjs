'use strict';
// NAC-1. Native exact arithmetic and the original serialized word compiler.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const root=__dirname,args=process.argv.slice(2),nativeRoot=args[args.indexOf('--native-root')+1],write=args.includes('--write');
assert.ok(args.includes('--native-root')&&nativeRoot);assert.ok(write!==args.includes('--check'));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex'),pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const[n,h]of Object.entries(pins.native_engine))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,n))),h,n);
for(const[n,h]of Object.entries(pins.packet_files))assert.equal(hash(fs.readFileSync(path.join(root,n))),h,n);
for(const[n,h]of Object.entries(pins.predecessor_files))assert.equal(hash(fs.readFileSync(path.resolve(root,'../..',n))),h,n);
const N=require(path.join(nativeRoot,'native_operator.cjs')),M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const f=N.F.of,mm=(...a)=>a.reduce(N.mul),tr=a=>a.reduce((s,r,i)=>s.add(r[i]),N.ZERO),real=z=>{assert.ok(z.turn.zero());return z.rad;},sc=a=>real(tr(a)).div(2);
const I=N.identity(2),R=N.matrix([[0,-1],[1,0]]),J=N.scale(N.matrix([[0,1],[1,0]]),N.IOTA),L=mm(R,J);
const norm2=a=>sc(mm(N.dagger(a),a)),chi=h=>f(1).sub(sc(h)).div(2),conj=(g,a)=>mm(g,a,N.dagger(g));
const checks=[],controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,value,detail={})=>{assert.ok(value,name);checks.push({name,passed:true,...detail});};
const no=(name,value,detail={})=>{assert.ok(value,name);controls.push({name,rejected:true,...detail});};
function system(name,tokens,rules=[],centralCount=0){
 for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)if(i<centralCount)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Declared real central coefficients'});
 const s=new M.Presentation({tokens,rules}),w=Object.fromEntries(tokens.map(t=>[t,s.word([t])])),audit=s.audit();
 ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w};
}
function prove(name,lhs,rhs){const s=lhs.system,out=s.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}

// Universal polynomial identities, not interpolation through test data.
{
 const{s,w:{m,a,b,c,d,A,B}}=system('finite_two_direction_loop',['m','a','b','c','d','A','B'],[
  {id:'AA',lhs:['A','A'],rhs:[[[],-1]],source:'Unit circular direction'},
  {id:'BB',lhs:['B','B'],rhs:[[[],-1]],source:'Unit circular direction'},
  {id:'BA',lhs:['B','A'],rhs:[[['A','B'],-1],[['m'],-2]],source:'Native direction overlap'},
  {id:'bb',lhs:['b','b'],rhs:[[[],1],[['a','a'],-1]],source:'First circular phase'},
  {id:'dd',lhs:['d','d'],rhs:[[[],1],[['c','c'],-1]],source:'Second circular phase'}
 ],5),one=s.one();
 const U=a.plus(b.times(A)),Ud=a.minus(b.times(A)),V=c.plus(d.times(B)),Vd=c.minus(d.times(B));
 const H=Vd.times(Ud).times(V).times(U),Hd=Ud.times(Vd).times(U).times(V),K=one.minus(m.times(m)),q=b.times(b).times(d).times(d).times(K),br=A.times(B).minus(B.times(A));
 prove('unitarity:first_transport',Ud.times(U),one);
 prove('unitarity:second_transport',Vd.times(V),one);
 prove('directions:commutator_strength',br.times(br).scale(-1),K.scale(4));
 prove('loop:exact_symmetric_part',H.plus(Hd),one.scale(2).minus(q.scale(4)));
 prove('loop:transport_commutator',U.times(V).minus(V.times(U)),b.times(d).times(br));
 prove('one_axis:quadratic_relation',U.times(U).minus(a.times(U).scale(2)).plus(one),s.zero());
 prove('one_axis:cubic_recurrence',U.times(U).times(U),a.times(a).scale(4).minus(one).times(U).minus(a.scale(2)));
}
{
 const{s,w:{t,u,k}}=system('angular_polynomials',['t','u','k'],[],3),one=s.one();
 const P=one.minus(t.times(t)),Q=t.scale(2),D=one.plus(t.times(t));
 prove('cayley:unit_pair',P.times(P).plus(Q.times(Q)),D.times(D));
 prove('cayley:cosine_derivative_numerator',t.scale(-2).times(D).minus(P.times(t.scale(2))),Q.scale(-2));
 prove('cayley:sine_derivative_numerator',D.scale(2).minus(Q.times(t.scale(2))),P.scale(2));
 const Pu=one.minus(u.times(u)),Qu=u.scale(2),Du=one.plus(u.times(u)),E=one.minus(t.times(u)),F=t.plus(u);
 prove('addition:cosine_numerator',P.times(Pu).minus(Q.times(Qu)),E.times(E).minus(F.times(F)));
 prove('addition:sine_numerator',P.times(Qu).plus(Q.times(Pu)),E.times(F).scale(2));
 prove('addition:denominator',D.times(Du),E.times(E).plus(F.times(F)));
 prove('defect:finite_doubling',one.minus(one.minus(k).times(one.minus(k)).scale(2).minus(one)),k.scale(4).minus(k.times(k).scale(2)));
}
{
 const{s,w:{k,b,d,e,g}}=system('calibration_cross_multiplication',['k','b','d','e','g'],[],5);
 const q0=k.times(b).times(b).times(d).times(d),q1=k.times(e).times(e).times(g).times(g);
 prove('calibration:cross_angle_prediction',q0.times(e).times(e).times(g).times(g),q1.times(b).times(b).times(d).times(d));
}

const pair=t=>{t=f(t);const D=f(1).add(t.pow(2));return{c:f(1).sub(t.pow(2)).div(D),s:t.mul(2).div(D)};};
const transport=(a,p)=>N.add(N.scale(I,p.c),N.scale(a,p.s));
const axes=[{name:'perpendicular',mu:f(0),B:J},{name:'tilted',mu:f('3/5'),B:N.add(N.scale(R,'3/5'),N.scale(J,'4/5'))},{name:'opposite_tilt',mu:f('-3/5'),B:N.add(N.scale(R,'-3/5'),N.scale(J,'4/5'))},{name:'coincident',mu:f(1),B:R}];
const steps=[['1/2','1/3'],['1','1'],['-1/2','3/4']],samples=[];
const frame=N.cayley(L,'2/5');
for(const axis of axes){
 const{B,mu}=axis,K=f(1).sub(mu.pow(2)),br=N.commutator(R,B);
 ck('axis:'+axis.name+':contract',N.equal(mm(B,B),N.scale(I,-1))&&N.equal(N.dagger(B),N.scale(B,-1))&&N.equal(N.add(mm(R,B),mm(B,R)),N.scale(I,mu.mul(-2))));
 ck('axis:'+axis.name+':strength',norm2(br).eq(K.mul(4)));
 for(const [t,u]of steps){
  const a=pair(t),b=pair(u),U=transport(R,a),V=transport(B,b),H=mm(N.dagger(V),N.dagger(U),V,U),value=chi(H),pred=K.mul(a.s.pow(2)).mul(b.s.pow(2)),name=axis.name+':'+t+':'+u;
  ck(name+':Cayley_transport',N.equal(U,N.cayley(R,t))&&N.equal(V,N.cayley(B,u)));
  ck(name+':unitary_loop',N.equal(mm(N.dagger(H),H),I));
  ck(name+':finite_angle_law',value.eq(pred)&&norm2(N.sub(I,H)).div(4).eq(value));
  ck(name+':calibration_recovery',value.div(a.s.pow(2).mul(b.s.pow(2))).eq(K));
  const Ug=conj(frame,U),Vg=conj(frame,V),Hg=mm(N.dagger(Vg),N.dagger(Ug),Vg,Ug);
  ck(name+':simultaneous_frame_invariance',N.equal(Hg,conj(frame,H))&&chi(Hg).eq(value));
  samples.push({axis:axis.name,mu,parameters:[t,u],first_phase:a,second_phase:b,holonomy:H,return_defect:value,calibration:K});
 }
}
const perpendicular=samples.find(x=>x.axis==='perpendicular'),tilted=samples.find(x=>x.axis==='tilted'),opposite=samples.find(x=>x.axis==='opposite_tilt'),parallel=samples.find(x=>x.axis==='coincident');
ck('benchmarks:exact_values',perpendicular.return_defect.eq('144/625')&&tilted.return_defect.eq('2304/15625')&&parallel.return_defect.zero());
no('omission:angles_alone_determine_loop',!perpendicular.return_defect.eq(tilted.return_defect)&&!tilted.return_defect.eq(parallel.return_defect));
no('omission:scalar_defect_determines_oriented_holonomy',tilted.return_defect.eq(opposite.return_defect)&&!N.equal(tilted.holonomy,opposite.holonomy));
no('omission:scalar_defect_determines_loop_orientation',chi(N.dagger(perpendicular.holonomy)).eq(perpendicular.return_defect)&&!N.equal(perpendicular.holonomy,N.dagger(perpendicular.holonomy)));
const central=N.scale(I,-1),V=N.cayley(J,'1/3'),centralLoop=mm(N.dagger(V),central,V,central);
no('inference:zero_finite_return_implies_commuting_generators',N.equal(centralLoop,I)&&!N.isZero(N.commutator(R,J)));
no('division:zero_sine_calibration_rejected',(()=>{try{f(0).div(0);return false;}catch{return true;}})());
const antipodal=samples.find(x=>x.axis==='perpendicular'&&x.parameters[0]==='1');
ck('loop:antipode_is_regular',N.equal(antipodal.holonomy,N.scale(I,-1))&&antipodal.return_defect.eq(1));
no('inference:antipodal_Cayley_chart_pole_is_a_transport_singularity',(()=>{try{N.inverse(N.add(I,antipodal.holonomy));return false;}catch{return antipodal.return_defect.eq(1);}})());
const Fconstant=N.commutator(R,J),Fvariable=N.sub(Fconstant,N.scale(L,2));
no('omission:instantaneous_axes_determine_variable_field_curvature',!N.isZero(Fconstant)&&N.isZero(Fvariable));
results.loops={samples,calibration_formula:'K=chi(theta0,phi0)/(sin(theta0)^2 sin(phi0)^2)=1-mu^2 for a nonzero-sine reference and fixed unit directions.',finite_prediction:'chi(theta,phi)=K sin(theta)^2 sin(phi)^2',angular_area_correction:'chi/A_ang^2=sinc(theta)^2 sinc(phi)^2; real completion and limiting proof are in THEOREM.md.',variable_field_control:{instantaneous_axes:[R,J],constant_curvature:Fconstant,derivative_partial_y_Ax:N.scale(L,2),resulting_curvature:Fvariable}};

// Native recurrence independently compared against exact matrix powers.
const U=N.cayley(R,'1/2'),c=sc(U),kappa=f(1).sub(c),multiples=[];
let prev=f(1),current=c;
for(let n=0;n<=8;n++){
 let poly;if(n===0)poly=prev;else if(n===1)poly=current;else{const next=c.mul(2).mul(current).sub(prev);prev=current;current=next;poly=next;}
 const actual=N.power(U,n),defect=f(1).sub(sc(actual));
 ck('repeat:'+n+':polynomial_prediction',sc(actual).eq(poly));
 multiples.push({multiple:n,holonomy:actual,scalar:poly,compact_defect:defect});
}
ck('repeat:doubled_and_tripled_defects',multiples[2].compact_defect.eq('32/25')&&multiples[3].compact_defect.eq('242/125'));
no('approximation:finite_defect_is_purely_quadratic',!multiples[2].compact_defect.eq(kappa.mul(4)));
ck('repeat:negative_power_defect',f(1).sub(sc(N.dagger(N.power(U,3)))).eq(multiples[3].compact_defect));
results.one_axis={reference:U,reference_defect:kappa,multiples,recurrence:'P0=1; P1=c; P_(n+1)=2 c P_n-P_(n-1); kappa_n=1-P_n(1-kappa).',calibration:'U(theta)=exp((theta/theta0)*L0), L0=theta0*A; theta0 != 0 and the generator increment/lift is retained.'};

// Fractional winding: exact pi/2 rotations, not numerical pi evaluations.
const windingCases=[];
for(const [ratio,quarterPower]of[['1',4],['1/2',2],['3/4',3]]){
 const r=f(ratio),phase=Array.from({length:8},(_,n)=>N.power(R,quarterPower*n));
 let matches=true;for(let i=0;i<8;i++)for(let j=0;j<8;j++)matches=matches&&(N.equal(phase[i],phase[j])===(r.mul(i-j).d===1n));
 ck('lift:'+ratio+':exact_residue_classification',matches);
 windingCases.push({rescaling:r,initial_endpoint:'I',angle_lifts:'2*pi*n',outputs:phase,required_winding_modulus:r.d.toString()});
}
no('omission:endpoint_fixes_fractional_transport',N.equal(N.power(R,0),N.power(R,4))&&!N.equal(N.power(R,0),N.power(R,2)));
results.winding={cases:windingCases,general_statement:'For r=p/q in lowest terms, q>0, equal rescaled transports occur exactly when the winding difference is divisible by q. For irrational r they require equal winding. Written proof uses the inherited 2*pi kernel.'};

const chosen=proofs.find(p=>p.certificate.steps.length),sys=systems.get(chosen.certificate.presentation),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['997','0'];
no('integrity:corrupted_rewrite_rejected',(()=>{try{sys.replay(bad);return false;}catch{return true;}})());
const payload={protocol:'NAC1_NATIVE_ANGULAR_CALIBRATION_V1',status:'PASS_FINITE_ANGLE_TRANSPORT_AND_CURVATURE_CALIBRATION',scope:'Fixed unit circular generators in the existing native sector; universal algebraic finite-loop identities, exact recurrence and matrix fixtures, and operation-specific angle-lift records. Completion and all-angle arguments are written. No physical coupling or complete gravitational dynamics is selected.',source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,exact_checks:checks.length,negative_controls:controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',h=hash(bytes),certificate=path.join(root,'CERTIFICATE.json'),expected=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(certificate,bytes);fs.writeFileSync(expected,h+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(certificate,'utf8'),bytes);assert.equal(fs.readFileSync(expected,'utf8'),h+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,mode:write?'generation':'read-only replay',...payload.summary,sha256:h,perpendicular_defect:perpendicular.return_defect,tilted_defect:tilted.return_defect,doubled_defect:multiples[2].compact_defect},null,2));
