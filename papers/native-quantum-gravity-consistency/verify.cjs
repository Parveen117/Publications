'use strict';
// QGC-1: unchanged native exact arithmetic and proof compiler.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const root=__dirname,args=process.argv.slice(2),nativeRoot=args[args.indexOf('--native-root')+1],write=args.includes('--write');
assert.ok(args.includes('--native-root')&&nativeRoot);assert.ok(write!==args.includes('--check'));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const[n,h]of Object.entries(pins.native_engine))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,n))),h,n);
for(const[n,h]of Object.entries(pins.packet_files))assert.equal(hash(fs.readFileSync(path.join(root,n))),h,n);
for(const[n,h]of Object.entries(pins.predecessor_files))assert.equal(hash(fs.readFileSync(path.resolve(root,'../..',n))),h,n);
const N=require(path.join(nativeRoot,'native_operator.cjs')),M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const f=N.F.of,C=N.Cut.of,mm=(...a)=>a.reduce(N.mul),tr=a=>a.reduce((s,row,i)=>s.add(row[i]),N.ZERO);
const real=z=>{assert.ok(z.turn.zero());return z.rad;},col=a=>N.matrix(a.map(x=>[x]));
const checks=[],controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,value,detail={})=>{assert.ok(value,name);checks.push({name,passed:true,...detail});};
const no=(name,value,detail={})=>{assert.ok(value,name);controls.push({name,rejected:true,...detail});};
const I=N.identity(2),Z=N.matrix([[1,0],[0,-1]]),X=N.matrix([[0,1],[1,0]]),R=N.matrix([[0,-1],[1,0]]),Y=N.scale(R,N.IOTA);
const cut=(z,x,y=0)=>N.scale(N.add(I,N.add(N.scale(Z,z),N.add(N.scale(X,x),N.scale(Y,y)))),'1/2');
const Pp=cut(1,0),Pm=cut(-1,0),p=cut('3/5','4/5'),q=cut('3/5','-4/5');
const moments=rho=>[I,Z,X,Y].map(b=>real(tr(mm(b,rho)))),conj=(u,rho)=>mm(u,rho,N.dagger(u));
const sum=xs=>xs.reduce(N.add,N.zeros(xs[0].length,xs[0][0].length));
const F=p=>{const[n,z,x]=moments(p);assert.ok(n.eq(1));const h=N.add(N.scale(Z,z.div(-10)),N.scale(X,x.div(-47)));return N.commutator(N.commutator(p,h),p);};

// 1. The chosen nonlinear relaxation fails mixture equivalence if promoted
// to an ensemble-independent quantum evolution on an unknown source.
const ea=[{weight:f('4/5'),state:Pp},{weight:f('1/5'),state:Pm}],eb=[{weight:f('1/2'),state:p},{weight:f('1/2'),state:q}];
const average=(ensemble,fn)=>sum(ensemble.map(e=>N.scale(fn(e.state),e.weight)));
const rho=average(ea,x=>x),fa=average(ea,F),fb=average(eb,F);
ck('ensemble:same_density',N.equal(rho,average(eb,x=>x))&&N.equal(rho,N.matrix([['4/5',0],[0,'1/5']])));
for(const[name,source]of[['z_plus',Pp],['z_minus',Pm],['tilted_plus',p],['tilted_minus',q]]){
 ck('ensemble:'+name+':pure',N.equal(mm(source,source),source)&&tr(source).eq(1));
 const[,z,x]=moments(source),g=f('37/470');
 ck('ensemble:'+name+':gradient_formula',N.equal(F(source),N.add(N.scale(Z,g.mul(z).mul(x.pow(2))),N.scale(X,g.mul(z.pow(2)).mul(x).neg()))));
}
ck('ensemble:first_ensemble_stationary',N.isZero(fa));
ck('ensemble:second_ensemble_nonzero',N.equal(fb,N.scale(Z,'888/29375')));
no('ensemble:affine_extension_preserving_all_pure_velocities',!N.equal(fa,fb));
results.ensemble={ensembles:{ea,eb},common_density:rho,average_velocity:{ea:fa,eb:fb},Z_readout_velocity_difference:real(tr(mm(Z,N.sub(fb,fa)))),contract:'Fixed beta=1, fixed source amount N=1, classical randomization and density-matrix mixture equivalence. Refutes promoting EGR-1 selected gradient relaxation to an ensemble-independent affine evolution. Does not invalidate the certified gradient identities or prove signalling without additional operational premises.'};

// 2. A fixed exponential filter has a unitary implementation with a full
// outcome ledger. It is a different operational model, not a cure that keeps F.
const a=f('3/5'),b=f('4/5'),Ks=N.add(Pp,N.scale(Pm,a)),Kf=N.scale(Pm,b),rotation=N.add(N.scale(I,a),N.scale(R,b));
const W=N.add(N.kron(Pp,I),N.kron(Pm,rotation)),E=r=>N.add(conj(Ks,r),conj(Kf,r));
ck('instrument:circle_coefficients',a.pow(2).add(b.pow(2)).eq(1));
ck('instrument:kraus_completeness',N.equal(N.add(mm(N.dagger(Ks),Ks),mm(N.dagger(Kf),Kf)),I));
ck('instrument:joint_unitarity',N.equal(mm(N.dagger(W),W),N.identity(4))&&N.equal(mm(W,N.dagger(W)),N.identity(4)));
const block=(outcome,input)=>N.matrix([0,1].map(i=>[0,1].map(j=>W[2*i+outcome][2*j+input])));
ck('instrument:dilation_success_block',N.equal(block(0,0),Ks));
ck('instrument:dilation_failure_block',N.equal(block(1,0),Kf));
const samples=[p,q,Pp,Pm,cut(0,1),cut(0,0,1),N.scale(I,'1/2'),N.scale(p,7)],ledgers=[];
const frame=N.add(N.scale(I,'3/5'),N.scale(R,'4/5'));
for(let j=0;j<samples.length;j++){
 const input=samples[j],[n,z,x,y]=moments(input),success=conj(Ks,input),failure=conj(Kf,input),ns=n.add(z).add(a.pow(2).mul(n.sub(z))).div(2),nf=b.pow(2).mul(n.sub(z)).div(2);
 const ms=moments(success),mf=moments(failure),out=E(input);
 ck('instrument:'+j+':success_moments',ms.every((v,i)=>v.eq([ns,n.add(z).sub(a.pow(2).mul(n.sub(z))).div(2),a.mul(x),a.mul(y)][i])));
 ck('instrument:'+j+':failure_moments',mf.every((v,i)=>v.eq([nf,nf.neg(),f(0),f(0)][i])));
 ck('instrument:'+j+':total_moments',moments(out).every((v,i)=>v.eq([n,z,a.mul(x),a.mul(y)][i])));
 ck('instrument:'+j+':amount_preserved',tr(out).eq(n));
 const rotatedOut=N.add(conj(conj(frame,Ks),conj(frame,input)),conj(conj(frame,Kf),conj(frame,input)));
 ck('instrument:'+j+':joint_frame_covariance',N.equal(rotatedOut,conj(frame,out)));
 ledgers.push({input,success,failure,success_moments:ms,failure_moments:mf,output:out});
}
ck('instrument:ensemble_independence',N.equal(average(ea,E),average(eb,E))&&N.equal(average(ea,E),E(rho)));
const success=ledgers[0].success,failure=ledgers[0].failure,ps=real(tr(success)),pf=real(tr(failure)),normalizedSuccess=N.scale(success,f(1).div(ps));
ck('instrument:exact_conditional_probability',ps.eq('109/125')&&pf.eq('16/125'));
ck('instrument:exact_conditional_direction',moments(normalizedSuccess).every((v,i)=>v.eq(['1','91/109','60/109','0'][i])));
ck('instrument:conditional_cut_remains_pure',N.equal(mm(normalizedSuccess,normalizedSuccess),normalizedSuccess));
no('instrument:normalization_retains_source_amount',!ps.eq(1)&&!real(tr(mm(Z,normalizedSuccess))).eq(real(tr(mm(Z,success)))));
no('instrument:successful_branch_is_the_unconditional_operation',!N.equal(normalizedSuccess,E(p)));
no('instrument:new_channel_preserves_old_gradient',!N.equal(N.scale(N.sub(conj(Z,p),p),'1/2'),F(p)),{note:'The new semigroup generator differs from the old state-dependent gradient. The ensemble theorem rules out an affine extension retaining all old pure-state velocities.'});
const repeated=[];
for(const m of[1,2,3,4]){
 const survive=N.power(Ks,m),failures=Array.from({length:m},(_,j)=>mm(Kf,N.power(Ks,j))),complete=N.add(mm(N.dagger(survive),survive),sum(failures.map(k=>mm(N.dagger(k),k))));
 ck('instrument:repeat_'+m+':completeness',N.equal(complete,I));
 ck('instrument:repeat_'+m+':exponential_power',N.equal(survive,N.add(Pp,N.scale(Pm,a.pow(m)))));
 const survived=conj(survive,p),all=N.add(survived,sum(failures.map(k=>conj(k,p))));let iterate=p;for(let j=0;j<m;j++)iterate=E(iterate);
 ck('instrument:repeat_'+m+':full_ledger_matches_channel',N.equal(all,iterate));
 ck('instrument:repeat_'+m+':success_probability',tr(survived).eq(f('4/5').add(f('1/5').mul(a.pow(2*m)))));
 repeated.push({m,success_probability:real(tr(survived)),survived,full_ledger:all});
}
const Ewith=(r,c)=>{const diagonal=N.add(mm(Pp,r,Pp),mm(Pm,r,Pm)),cross=N.add(mm(Pp,r,Pm),mm(Pm,r,Pp));return N.add(diagonal,N.scale(cross,c));};
for(let i=0;i<2;i++)for(let j=0;j<2;j++){
 const unit=N.zeros(2);unit[i][j]=N.ONE;
 ck('instrument:semigroup_basis_'+i+j,N.equal(Ewith(Ewith(unit,'3/5'),'5/13'),Ewith(unit,'3/13')));
 ck('instrument:generator_basis_'+i+j,N.equal(N.scale(N.commutator(Z,N.commutator(Z,unit)),'-1/4'),N.scale(N.sub(conj(Z,unit),unit),'1/2')));
}
const remoteTrace=r=>N.matrix([0,1].map(i=>[0,1].map(j=>r[i][j].add(r[i+2][j+2]))));
const remoteE=r=>N.add(conj(N.kron(Ks,I),r),conj(N.kron(Kf,I),r));
for(let i=0;i<4;i++)for(let j=0;j<4;j++){
 const unit=N.zeros(4);unit[i][j]=N.ONE;
 ck('instrument:remote_marginal_matrix_basis_'+i+j,N.equal(remoteTrace(remoteE(unit)),remoteTrace(unit)));
}
results.instrument={a,b,Ks,Kf,joint_unitary:W,ledgers,normalized_success:normalizedSuccess,repeated,continuous_contract:'a=exp(-s), s>=0 is a dimensionless operation parameter. E_s off-diagonal factor exp(-s), generator (Z*rho*Z-rho)/2. Repeated reduced operations use fresh ancillas; their joint unitary dilation is information preserving before discarding registers.',source_readout_contract:'For fixed alpha/gamma, y1=(alpha/gamma)*z/10 and y2=(alpha/gamma)*x/47. Preserve branch amounts in these moments. The unchanged z and contracted x are a readout prediction, not a derived dynamical law of spacetime.'};

// 3. Exact first-response operators and the necessary added noise of a
// commuting, unbiased joint pointer. chi=1 fixtures; scaling is explicit.
const A1=N.scale(Z,'1/10'),A2=N.scale(X,'1/47'),comm=N.commutator(A1,A2);
ck('response:noncommutation',N.equal(comm,N.scale(Y,N.IOTA.mul('1/235'))));
ck('response:squares',N.equal(mm(A1,A1),N.scale(I,'1/100'))&&N.equal(mm(A2,A2),N.scale(I,'1/2209')));
const etaZ=f('3/5'),etaX=f('4/5'),effects=[],kraus=[];
for(const s of[1,-1])for(const t of[1,-1]){
 const projector=cut(etaZ.mul(s),etaX.mul(t)),effect=N.scale(projector,'1/2'),k=N.scale(projector,C(['1/2','1/2']));
 ck('pointer:'+s+'_'+t+':positive_effect_factor',N.equal(mm(projector,projector),projector)&&N.equal(mm(N.dagger(k),k),effect));
 effects.push({s,t,effect,first:f(s).div(6),second:f(t).mul('5/188'),probability_on_p:real(tr(mm(effect,p)))});kraus.push(k);
}
ck('pointer:effects_sum_to_identity',N.equal(sum(effects.map(e=>e.effect)),I));
const firstMoment=sum(effects.map(e=>N.scale(e.effect,e.first))),secondMoment=sum(effects.map(e=>N.scale(e.effect,e.second)));
ck('pointer:unbiased_first_moments',N.equal(firstMoment,A1)&&N.equal(secondMoment,A2));
const N1=N.sub(sum(effects.map(e=>N.scale(e.effect,e.first.pow(2)))),mm(A1,A1)),N2=N.sub(sum(effects.map(e=>N.scale(e.effect,e.second.pow(2)))),mm(A2,A2));
ck('pointer:exact_added_noises',N.equal(N1,N.scale(I,'4/225'))&&N.equal(N2,N.scale(I,'9/35344')));
ck('pointer:norm_bound_saturated',f('4/225').mul('9/35344').eq(f('1/470').pow(2)));
const V=kraus.flat(),diag=values=>N.matrix(values.map((v,i)=>values.map((_,j)=>i===j?v:0))),B1=diag(effects.flatMap(e=>[e.first,e.first])),B2=diag(effects.flatMap(e=>[e.second,e.second]));
ck('pointer:native_isometry',N.equal(mm(N.dagger(V),V),I));
ck('pointer:commuting_extension',N.isZero(N.commutator(B1,B2)));
ck('pointer:compressed_observables',N.equal(mm(N.dagger(V),B1,V),A1)&&N.equal(mm(N.dagger(V),B2,V),A2));
const complement=N.sub(N.identity(8),mm(V,N.dagger(V))),L1=mm(complement,B1,V),L2=mm(complement,B2,V);
ck('pointer:noise_leakage_factors',N.equal(mm(N.dagger(L1),L1),N1)&&N.equal(mm(N.dagger(L2),L2),N2));
ck('pointer:compressed_commutator_identity',N.equal(comm,N.sub(mm(N.dagger(L2),L1),mm(N.dagger(L1),L2))));
const sharpEffect=N.scale(N.add(I,N.add(Z,X)),'1/4'),badProbability=real(tr(mm(sharpEffect,cut('-3/5','-4/5'))));
no('pointer:two_sharp_channels_in_naive_joint_effect',badProbability.eq('-1/10'));
no('pointer:zero_added_noise_for_both_noncommuting_responses',!N.isZero(comm));
results.joint_readout={chi_fixture:'1',response_operators:[A1,A2],commutator:comm,effects,isometry:V,commuting_pointers:[B1,B2],noise_operators:[N1,N2],noise_product:'1/220900',general_bound:'||N1||*||N2|| >= chi^4/220900, for commuting Hermitian joint pointers and an isometry V with V^dagger B_i V=A_i exactly for all source states.',scope:'Standard quantum joint-readout structure specialized to the two certified response coefficients. These are pulled-back source observables, not a construction of an independent quantum spacetime.'};

// 4. Reuse the ten-edge response, derive all volume slopes from the cell
// Gram matrices, and recover the nonlinear source coefficient from two.
const prior=n=>JSON.parse(fs.readFileSync(path.resolve(root,'../',n,'CERTIFICATE.json'),'utf8'));
const egr1=prior('native-exponential-gravity-response'),egr2=prior('native-nonlinear-gravity-continuation'),crc=prior('native-cut-record-completeness');
ck('gravity:predecessor_status',egr2.status==='PASS_EXACT_NONLINEAR_JET_AND_LOCAL_BRANCH_HYPOTHESES'&&crc.status==='PASS_DECLARED_OPERATION_CLOSURE_AND_SOURCE_LAW_WITNESSES');
const g=egr2.results.geometry,A=N.matrix(g.mode_readout),S=N.matrix(g.section),G=N.matrix(g.gauge),K=N.matrix(g.stiffness),eta0=N.matrix(egr2.results.response.plus.selected_flat_displacement),d0=mm(S,N.inverse(K),col(['3/5','4/5']));
const geo=egr1.results.geometry;
function rates(d){const idx=Object.fromEntries(g.interior_edges.map((e,i)=>[e,i])),dv=(a,b)=>{const key=[a,b].sort().join(':');return idx[key]===undefined?f(0):real(d[idx[key]][0]);};return geo.cells.map((cell,i)=>{const gd=N.matrix(cell.slice(1).map(a=>cell.slice(1).map(b=>dv(cell[0],a).add(dv(cell[0],b)).sub(a===b?0:dv(a,b)).div(2))));return real(tr(mm(N.inverse(N.matrix(geo.cell_records[i].gram)),gd))).div(2);});}
const aAll=rates(d0),bAll=rates(mm(G,eta0)),classes=[0,4,5],avec=classes.map(i=>aAll[i]),bvec=classes.map(i=>bAll[i]);
ck('prediction:uncorrected_volume_slopes',avec.every((v,i)=>v.eq(['3/10',0,'64/47'][i])));
ck('prediction:flat_displacement_volume_slopes',bvec.every((v,i)=>v.eq(['-456869/4241280','4410869/1060320','-531631/1060320'][i])));
for(let i=0;i<9;i++)ck('prediction:volume_classes_'+i,aAll[i].eq(avec[i<4?0:i===4?1:2])&&bAll[i].eq(bvec[i<4?0:i===4?1:2]));
const ratio=bvec[0].div(bvec[1]),c1=avec[2].div(avec[0]),c2=bvec[2].sub(c1.mul(bvec[0])).div(bvec[1]),D=avec[0].mul(bvec[1]).sub(avec[1].mul(bvec[0]));
ck('prediction:two_calibrations_nondegenerate',!D.zero());
ck('prediction:heldout_coefficients',c1.eq('640/141')&&c2.eq('-1860931/621932529'));
function recover(v){const c=v[0].sub(ratio.mul(v[1])).div(avec[0]);if(c.zero())return{identified:false,reason:'zero source slope'};const t=v[1].div(bvec[1].mul(c));if(t.zero())return{identified:false,reason:'no finite kappa in this family'};return{identified:true,c,t,kappa:f(1).sub(f(1).div(t)),heldout_residual:v[2].sub(c1.mul(v[0])).sub(c2.mul(v[1]))};}
const predictions=[];
for(const kappa of[0,2,3,-1])for(const c of['1','7/13']){
 const t=f(1).div(f(1).sub(kappa)),d=N.add(d0,N.scale(mm(G,eta0),t)),direct=rates(d),v=avec.map((v,i)=>f(c).mul(v.add(bvec[i].mul(t)))),reconstruction=recover(v);
 ck('prediction:k'+kappa+'_c'+c+':direct_geometric_derivative',classes.every((j,i)=>direct[j].mul(c).eq(v[i])));
 ck('prediction:k'+kappa+'_c'+c+':parameter_recovery',reconstruction.identified&&reconstruction.kappa.eq(kappa)&&reconstruction.c.eq(c));
 ck('prediction:k'+kappa+'_c'+c+':heldout',reconstruction.heldout_residual.zero());
 if((kappa===0||kappa===2)&&c==='1')ck('prediction:k'+kappa+':previous_packet_reproduced',direct.every((v,i)=>v.eq(crc.results.gravity.cases.find(c=>c.kappa===kappa).volume_log_derivative[i])));
 predictions.push({kappa,c,slopes:v,reconstruction});
}
const altered=predictions[0].slopes.map((v,i)=>i===2?v.add('1/100'):v),alteredFit=recover(altered);
no('prediction:fit_two_slopes_validates_third',alteredFit.identified&&alteredFit.kappa.eq(0)&&alteredFit.c.eq(1)&&alteredFit.heldout_residual.eq('1/100'));
no('prediction:zero_signal_identifies_kappa',!recover([f(0),f(0),f(0)]).identified);
no('prediction:zero_cap_slope_identifies_finite_kappa',!recover(avec).identified);
results.gravity_prediction={classes:['outer:four cells','cap:one cell','inner:four cells'],source:['3/5','4/5'],avec,bvec,calibration_determinant:D,outer_to_cap_b_ratio:ratio,heldout_coefficients:[c1,c2],predictions,contract:'v_i=d log(V_i)/dt at t=0, epsilon(0)=0, c=epsilon_prime(0)!=0, fixed source direction, fixed labelled triangulation, kappa!=1, local analytic branch. v=c*(avec+bvec/(1-kappa)). Two slopes calibrate c,kappa; the third and within-class equality are predictions. No finite-range error bound or physical gravitational-unit calibration is supplied by these exact derivatives.'};

// Native polynomial certificates supplement the written analytic proofs.
function system(name,tokens,rules=[],commute=false){if(commute)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Declared commuting real coefficients'});const s=new M.Presentation({tokens,rules}),w=Object.fromEntries(tokens.map(t=>[t,s.word([t])])),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w};}
function prove(name,lhs,rhs){const s=lhs.system,out=s.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{p,a,b,c}}=system('projector_filter',['p','a','b','c'],[{id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'P squared equals P'},{id:'circle',lhs:['b','b'],rhs:[[[],1],[['a','a'],-1]],source:'a squared plus b squared equals one'}],true),one=s.one(),m=one.minus(p),k=p.plus(a.times(m)),l=b.times(m);
 prove('filter:complement',m.times(m),m);
 prove('filter:completeness',k.times(k).plus(l.times(l)),one);
 prove('filter:success_plus_sector',k.times(p),p);
 prove('filter:failure_plus_sector',l.times(p),s.zero());
 prove('filter:composition',k.times(p.plus(c.times(m))),p.plus(a.times(c).times(m)));
}
{
 const{s,w:{a,b,n,z,x,y}}=system('branch_ledger',['a','b','n','z','x','y'],[{id:'circle',lhs:['b','b'],rhs:[[[],1],[['a','a'],-1]],source:'a squared plus b squared equals one'}],true);
 const ns=n.plus(z).plus(a.times(a).times(n.minus(z))).scale('1/2'),nf=b.times(b).times(n.minus(z)).scale('1/2'),zs=n.plus(z).minus(a.times(a).times(n.minus(z))).scale('1/2');
 prove('ledger:amount',ns.plus(nf),n);prove('ledger:z',zs.minus(nf),z);
 prove('ledger:purity_residual',ns.times(ns).minus(zs.times(zs)).minus(a.times(x).times(a.times(x))).minus(a.times(y).times(a.times(y))),a.times(a).times(n.times(n).minus(z.times(z)).minus(x.times(x)).minus(y.times(y))));
 prove('ledger:total_coherence',a.times(x).plus(a.times(y).scale(N.IOTA)),a.times(x.plus(y.scale(N.IOTA))));
}
{
 const{s,w:{p,h}}=system('conditional_gradient',['p','h'],[{id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'P squared equals P'}]);
 const c=p.times(h).minus(h.times(p));
 prove('gradient:fixed_operator_identity',c.times(p).minus(p.times(c)),p.times(h).times(p).scale(2).minus(h.times(p)).minus(p.times(h)));
}
{
 const{s,w:{z,x}}=system('response_algebra',['z','x'],[{id:'z_square',lhs:['z','z'],rhs:[[[],1]],source:'Z squared equals one'},{id:'x_square',lhs:['x','x'],rhs:[[[],1]],source:'X squared equals one'},{id:'anticommute',lhs:['x','z'],rhs:[[['z','x'],-1]],source:'XZ=-ZX'}]);
 const a=z.scale('1/10'),b=x.scale('1/47'),rot=x.times(z);
 prove('response:commutator',a.times(b).minus(b.times(a)),rot.scale('-1/235'));
 prove('response:first_square',a.times(a),s.one().scale('1/100'));
 prove('response:second_square',b.times(b),s.one().scale('1/2209'));
 prove('response:real_rotation_square',rot.times(rot),s.one().scale(-1));
}
{
 const{s,w:{p,b,c}}=system('pointer_compression',['p','b','c'],[{id:'projection',lhs:['p','p'],rhs:[[['p'],1]],source:'P squared equals P'},{id:'commuting_pointers',lhs:['c','b'],rhs:[[['b','c'],1]],source:'Joint pointers commute'}]);
 const one=s.one(),a=p.times(b).times(p),d=p.times(c).times(p),left=a.times(d).minus(d.times(a)),leakbc=p.times(b).times(one.minus(p)).times(c).times(p),leakcb=p.times(c).times(one.minus(p)).times(b).times(p);
 prove('pointer:compressed_commutator_equals_leakage',left,leakcb.minus(leakbc));
 prove('pointer:positive_noise_factor_identity',p.times(b).times(b).times(p).minus(a.times(a)),p.times(b).times(one.minus(p)).times(b).times(p));
}
{
 const{s,w:{c,t}}=system('source_calibration',['c','t'],[],true),v=avec.map((a,i)=>c.times(s.one().scale(a).plus(t.scale(bvec[i]))));
 prove('prediction:heldout_slope',v[2],v[0].scale(c1).plus(v[1].scale(c2)));
 prove('prediction:recover_source_scale',v[0].minus(v[1].scale(ratio)),c.scale(avec[0]));
 prove('prediction:recover_c_times_t',v[1],c.times(t).scale(bvec[1]));
}
const chosen=proofs.find(p=>p.certificate.steps.length),sys=systems.get(chosen.certificate.presentation),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['997','0'];
no('integrity:corrupted_proof_rejected',(()=>{try{sys.replay(bad);return false;}catch{return true;}})());
const payload={protocol:'QGC1_NATIVE_QUANTUM_GRAVITY_CONSISTENCY_V1',status:'PASS_FINITE_QUANTUM_CONSISTENCY_AND_CONDITIONAL_PREDICTIONS',scope:'Exact ensemble obstruction, selected unitary instrument and probability ledger, simultaneous readout noise saturation, and nonlinear coupling reconstruction with a held-out geometric slope. General positivity, norm inequality and local-branch statements have written proofs under declared premises. Does not certify a completed physical quantum-gravity theory or experimental validation.',source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,exact_checks:checks.length,negative_controls:controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',h=hash(bytes),cert=path.join(root,'CERTIFICATE.json'),expected=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(cert,bytes);fs.writeFileSync(expected,h+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cert,'utf8'),bytes);assert.equal(fs.readFileSync(expected,'utf8'),h+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,mode:write?'generation':'read-only replay',...payload.summary,sha256:h,ensemble_velocity_difference:results.ensemble.Z_readout_velocity_difference,success_probability:ps,added_noise_product:results.joint_readout.noise_product,heldout_coefficients:[c1,c2]},null,2));
