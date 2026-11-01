'use strict';
// EGR-1: unchanged native arithmetic and proof compiler, geometric adapters only.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const root=__dirname,args=process.argv.slice(2),nativeRoot=args[args.indexOf('--native-root')+1],write=args.includes('--write');
assert.ok(args.includes('--native-root')&&nativeRoot);assert.ok(write!==args.includes('--check'));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const[n,h]of Object.entries(pins.native_engine))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,n))),h,n);
for(const[n,h]of Object.entries(pins.packet_files))assert.equal(hash(fs.readFileSync(path.join(root,n))),h,n);
for(const[n,h]of Object.entries(pins.predecessor_files))assert.equal(hash(fs.readFileSync(path.resolve(root,'../..',n))),h,n);
const N=require(path.join(nativeRoot,'native_operator.cjs')),M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const f=N.F.of,C=N.Cut.of,real=z=>{assert.ok(z.turn.zero());return z.rad;},mm=(...xs)=>xs.reduce(N.mul);
const checks=[],controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,value,data={})=>{assert.ok(value,name);checks.push({name,passed:true,...data});};
const no=(name,value,data={})=>{assert.ok(value,name);controls.push({name,rejected:true,...data});};
const col=a=>N.matrix(a.map(x=>[x])),cols=vs=>N.dagger(vs.map(v=>v.map(r=>r[0]))),tr=g=>g.reduce((s,r,i)=>s.add(r[i]),N.ZERO);
function det(g){return g.length===1?g[0][0]:g[0].reduce((s,z,j)=>s.add(z.mul(det(g.slice(1).map(r=>r.filter((_,k)=>j!==k)))).mul(j%2?-1:1)),N.ZERO);}
const positive=x=>!x.le(0),positiveGram=g=>g.every((_,k)=>positive(real(det(g.slice(0,k+1).map(r=>r.slice(0,k+1))))));
function comb(xs,k){if(!k)return [[]];return xs.flatMap((x,i)=>comb(xs.slice(i+1),k-1).map(t=>[x,...t]));}
const key=vs=>vs.slice().sort().join(':'),zero=col([0,0,0,0]);
const coords={A:zero,B:col([1,0,0,0]),C:col([0,1,0,0]),D:col([0,0,1,0]),E:col([0,0,0,1]),O:col(['1/5','1/5','1/5','1/5']),P:col(['6/25','6/25','6/25','6/25'])};
const boundaryVertices=['A','B','C','D','E'],first=boundaryVertices.map((_,i)=>['O',...boundaryVertices.filter((_,j)=>j!==i)]),refined=first[0];
const cells=first.slice(1).concat(refined.map((_,i)=>['P',...refined.filter((_,j)=>j!==i)]));
const edgeKeys=[...new Set(cells.flatMap(c=>comb(c,2).map(key)))].sort(),interiorEdges=edgeKeys.filter(e=>e.includes('O')||e.includes('P'));
const qi=Object.fromEntries(edgeKeys.map(e=>{const[a,b]=e.split(':');return[e,N.energy(N.sub(coords[a],coords[b]))];}));
const triangles=[...new Set(cells.flatMap(c=>comb(c,3).map(key)))].sort(),interiorTriangles=triangles.filter(e=>e.includes('O')||e.includes('P'));
const tetra=new Map();for(const cell of cells)for(const face of comb(cell,4)){const k=key(face);tetra.set(k,(tetra.get(k)||0)+1);}
const boundaryTetra=[...tetra].filter(([,n])=>n===1).map(([k])=>k);
ck('complex:cell_count',cells.length===9);ck('complex:interior_edge_count',interiorEdges.length===10);ck('complex:interior_hinge_count',interiorTriangles.length===20);
ck('complex:original_boundary_preserved',JSON.stringify(boundaryTetra.slice().sort())===JSON.stringify(comb(boundaryVertices,4).map(key).sort()));
ck('complex:manifold_face_pairing',[...tetra.values()].every(n=>n===1||n===2));
ck('complex:interior_edge_definition',interiorEdges.every(e=>!boundaryTetra.some(t=>e.split(':').every(v=>t.split(':').includes(v)))));
const area2=(h,q=qi)=>{const[a,b,c]=h.split(':'),x=q[key([a,b])],y=q[key([a,c])],z=q[key([b,c])];return x.mul(y).sub(x.add(y).sub(z).pow(2).div(4)).div(4);};
const area2prime=(h,e)=>{if(!comb(h.split(':'),2).map(key).includes(e))return f(0);const rest=comb(h.split(':'),2).map(key).filter(k=>k!==e);return qi[rest[0]].add(qi[rest[1]]).sub(qi[e]).div(8);};
function affine(g){const out=N.zeros(5);for(let i=1;i<5;i++)for(let j=1;j<5;j++)out[i][j]=g[i-1][j-1];for(let i=1;i<5;i++){out[i][0]=g[i-1].reduce((s,z)=>s.add(z),N.ZERO).neg();out[0][i]=out[i][0];}out[0][0]=g.flat().reduce((s,z)=>s.add(z),N.ZERO);return out;}
function gramFromLengths(cell,q){return N.matrix(cell.slice(1).map(a=>cell.slice(1).map(b=>q[key([cell[0],a])].add(q[key([cell[0],b])]).sub(a===b?0:q[key([a,b])]).div(2))));}
const curvature=N.zeros(interiorTriangles.length,interiorEdges.length),hessian=N.zeros(interiorEdges.length),cellRecords=[];
for(let ci=0;ci<cells.length;ci++){
 const cell=cells[ci],E=cols(cell.slice(1).map(v=>N.sub(coords[v],coords[cell[0]]))),G=mm(N.dagger(E),E),inv=N.inverse(G),B=affine(inv),volume=real(det(E)).abs().div(24);
 ck('cell:'+ci+':full_dimension',N.rank(E)===4&&positiveGram(G));ck('cell:'+ci+':shared_lengths',N.equal(G,gramFromLengths(cell,qi)));
 ck('cell:'+ci+':volume',volume.eq(ci<4?'1/120':'1/600'));
 const derivatives=[];
 for(let ei=0;ei<interiorEdges.length;ei++){
  const edge=interiorEdges[ei];if(!edge.split(':').every(v=>cell.includes(v)))continue;
  const qprime=(a,b)=>key([a,b])===edge?f(1):f(0),Gp=N.matrix(cell.slice(1).map(a=>cell.slice(1).map(b=>qprime(cell[0],a).add(qprime(cell[0],b)).sub(a===b?0:qprime(a,b)).div(2))));
  const invp=N.scale(mm(inv,Gp,inv),-1),Bp=affine(invp);ck('cell:'+ci+':inverse_derivative:'+edge,N.isZero(N.add(mm(Gp,inv),mm(G,invp))));
  let sch=f(0);const angleTerms=[];
  for(let i=0;i<5;i++)for(let j=i+1;j<5;j++){
   const h=key(cell.filter((_,k)=>k!==i&&k!==j)),hi=interiorTriangles.indexOf(h);
   const k=real(Bp[i][j].sub(B[i][j].mul(Bp[i][i].div(B[i][i]).add(Bp[j][j].div(B[j][j]))).mul('1/2')));sch=sch.add(k);
   if(ei===interiorEdges.findIndex(e=>e.split(':').every(v=>cell.includes(v)))){
    const minor=real(B[i][i].mul(B[j][j]).sub(B[i][j].mul(B[i][j])));
    ck('cell:'+ci+':hinge_area:'+h,area2(h).eq(volume.pow(2).mul(144).mul(minor)));
   }
   if(hi>=0)curvature[hi][ei]=curvature[hi][ei].sub(C(volume.mul(12).mul(k)));
   angleTerms.push({hinge:h,area_times_angle_derivative:volume.mul(12).mul(k).toString()});
  }
  ck('cell:'+ci+':Schlafli:'+edge,sch.zero());derivatives.push({edge,angleTerms});
 }
 cellRecords.push({cell,gram:G,volume:volume.toString(),derivatives});
}
ck('complex:total_volume',cellRecords.reduce((s,c)=>s.add(c.volume),f(0)).eq('1/24'));
for(let i=0;i<interiorEdges.length;i++)for(let j=0;j<interiorEdges.length;j++)for(let h=0;h<interiorTriangles.length;h++)hessian[i][j]=hessian[i][j].add(curvature[h][j].mul(area2prime(interiorTriangles[h],interiorEdges[i]).div(area2(interiorTriangles[h]).mul(2))));
const gauge=N.matrix(interiorEdges.map(edge=>{const[a,b]=edge.split(':');return['O','P'].flatMap(v=>[0,1,2,3].map(k=>real(coords[a][k][0].sub(coords[b][k][0])).mul((a===v?1:0)-(b===v?1:0)).mul(2)));}));
const A=N.matrix([interiorEdges.map(k=>k.split(':').includes('O')&&!k.includes('P')?1:0),interiorEdges.map(k=>k.includes('P')?1:(k.includes('O')&&!k.includes('A')?'-1/5':0))]);
const piv=[interiorEdges.indexOf('A:O'),interiorEdges.indexOf('O:P')],section=N.zeros(10,2);piv.forEach((j,i)=>section[j][i]=N.ONE);
const K=mm(N.dagger(section),hessian,section),Kinv=N.inverse(K),L=mm(curvature,section);
ck('geometry:Hessian_symmetry',N.equal(hessian,N.dagger(hessian)));
ck('geometry:eight_vertex_displacement_directions',N.rank(gauge)===8);
ck('geometry:two_curvature_modes',N.rank(hessian)===2&&N.rank(curvature)===2&&N.rank(A)===2);
ck('geometry:Hessian_annihilates_gauge',N.isZero(mm(hessian,gauge)));
ck('geometry:curvature_annihilates_gauge',N.isZero(mm(curvature,gauge)));
ck('geometry:observed_modes_annihilate_gauge',N.isZero(mm(A,gauge)));
ck('geometry:section',N.equal(mm(A,section),N.identity(2)));
ck('geometry:reduced_stiffness',N.equal(K,N.matrix([[10,0],[0,47]])));
ck('geometry:complete_Hessian_factorization',N.equal(hessian,mm(N.dagger(A),K,A)));
ck('geometry:complete_curvature_factorization',N.equal(curvature,mm(L,A)));
ck('geometry:physical_mode_detects_curvature',N.rank(L)===2);
no('geometry:ten_edges_are_not_ten_curvature_modes',N.rank(hessian)!==10);
no('geometry:ungauged_inverse_is_invalid',(()=>{try{N.inverse(hessian);return false;}catch{return true;}})());
const badSource=col(interiorEdges.map((_,i)=>i===0?1:0));no('source:arbitrary_edge_force_can_violate_vertex_constraints',!N.isZero(mm(N.dagger(gauge),badSource)));
const shifted=Object.fromEntries(Object.entries(coords).map(([k,v])=>[k,N.add(v,col(['1/3','-2/5','1/7','2/9']))]));
ck('geometry:translation_invariance',edgeKeys.every(k=>{const[a,b]=k.split(':');return N.energy(N.sub(shifted[a],shifted[b])).eq(qi[k]);}));
const gaugeShift=col([1,2,-1,0,-2,1,0,1]),dq=col([1,2,3,4,5,6,7,8,9,10]);
ck('geometry:mode_readout_is_linearly_gauge_invariant',N.equal(mm(A,N.add(dq,mm(gauge,gaugeShift))),mm(A,dq)));
results.geometry={coordinates:coords,cells,interior_edges:interiorEdges,interior_hinges:interiorTriangles,boundary_tetrahedra:boundaryTetra,edge_squares:qi,cell_records:cellRecords,gauge,gauge_rank:8,hessian,hessian_rank:2,mode_readout:A,section,reduced_stiffness:K,reduced_inverse:Kinv,area_weighted_deficit_derivative:curvature,curvature_from_modes:L,scope:'Exact first variation of deficits and Hessian at this flat nine-simplex geometry. No nonlinear gauge closure or physical graviton count is asserted.'};

function system(name,tokens,rules=[],commute=false){
 if(commute)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'central real coefficients'});
 const S=new M.Presentation({tokens,rules}),w=Object.fromEntries(tokens.map(t=>[t,S.word([t])])),audit=S.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});systems.set(S.hash(),S);return{S,w};
}
function prove(name,lhs,rhs){const S=lhs.system,out=S.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',S.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{S,w:{p,h}}=system('cut_gradient',['p','h'],[{id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'P^2=P'}]);
 const v=M.bracket(p,h),d=M.bracket(v,p);
 prove('flow:projector_tangent',p.times(d).plus(d.times(p)),d);
 prove('flow:minimal_generator',p.times(v).plus(v.times(p)),v);
 prove('flow:gradient_expansion',d,p.times(h).times(p).scale(2).minus(p.times(h)).minus(h.times(p)));
 prove('flow:scalar_shift_has_no_effect',M.bracket(p,h.plus(S.one().scale(997))),v);
}
{
 const{S,w:{r}}=system('exponential_bridge',['r'],[{id:'quarter_turn',lhs:['r','r'],rhs:[[[],-1]],source:'R^2=-I'}]);
 const series=M.formalExp(r,10);let power=S.one(),factor=f(1);
 for(let n=0;n<=10;n++){if(n){power=power.times(r);factor=factor.mul(n);}const expected=n%2?r.scale(n%4===1?1:-1):S.one().scale(n%4===0?1:-1);prove('exponential:coefficient_'+n,series[n],expected.scale(f(1).div(factor)));}
}
{
 const{S,w:{a,p}}=system('exponential_of_cut',['a','p'],[{id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'P^2=P'}],true);
 const E=S.one().plus(a.minus(S.one()).times(p));
 prove('exponential:cut_exponential_left',E.times(p),a.times(p));
 prove('exponential:cut_exponential_right',p.times(E),a.times(p));
 prove('exponential:cut_exponential_cannot_move_its_own_cut',M.bracket(E,p),S.zero());
}
{
 const{S,w:{c,s,b}}=system('source_alignment',['c','s','b'],[],true),gap=f('1/10').sub('1/47');
 const dc=b.times(c).times(s).times(s).scale(gap.mul(2)),ds=b.times(c).times(c).times(s).scale(gap.mul(-2));
 prove('alignment:circle_invariant',c.times(dc).plus(s.times(ds)).scale(2),S.zero());
 prove('alignment:response_score_increases',c.times(dc).scale('1/5').plus(s.times(ds).scale('2/47')),b.times(c).times(c).times(s).times(s).scale(gap.pow(2).mul(4)));
 prove('alignment:ratio_exponential_numerator',ds.times(c).minus(s.times(dc)),b.times(c).times(s).times(c.times(c).plus(s.times(s))).scale(gap.mul(-2)));
 prove('alignment:potential_dissipation',b.times(c.times(dc).scale('1/5').plus(s.times(ds).scale('2/47'))).scale('-1/2'),b.times(b).times(c).times(c).times(s).times(s).scale(gap.pow(2).mul(-2)));
}
const I2=N.identity(2),Z=N.matrix([[1,0],[0,-1]]),X=N.matrix([[0,1],[1,0]]),R=N.matrix([[0,-1],[1,0]]),P0=N.matrix([[1,0],[0,0]]);
const cut=(c,s)=>N.scale(N.add(I2,N.add(N.scale(Z,c),N.scale(X,s))),'1/2');
const read=Pc=>col([real(tr(mm(Pc,Z))),real(tr(mm(Pc,X)))]);
const alpha=f('1/100'),gamma=f(1),gap=f('37/470'),sourceCases=[];
for(const[c,s]of[[1,0],[0,1],['3/5','4/5'],['5/13','12/13']])for(const amount of[1,2,4]){
 const Pc=cut(c,s),rho=N.scale(Pc,amount),r=read(Pc),chi=alpha.mul(amount).div(gamma),beta=alpha.pow(2).mul(f(amount).pow(2)).div(gamma),y=N.scale(mm(Kinv,r),chi),dq0=mm(section,y),j=N.scale(mm(N.dagger(A),r),alpha.mul(amount));
 const Hp=N.scale(N.add(N.scale(Z,real(y[0][0])),N.scale(X,real(y[1][0]))),alpha.mul(amount).neg()),V=N.commutator(Pc,Hp),dP=N.commutator(V,Pc);
 const name='source:'+c+':'+s+':N'+amount;
 ck(name+':quantized_cut',N.equal(mm(Pc,Pc),Pc)&&real(tr(Pc)).eq(1));
 ck(name+':amount_invariant',real(tr(rho)).eq(amount)&&real(tr(mm(rho,rho))).eq(f(amount).pow(2)));
 ck(name+':all_ten_metric_equations',N.equal(N.scale(mm(hessian,dq0),gamma),j));
 ck(name+':gauge_compatible_source',N.isZero(mm(N.dagger(gauge),j)));
 ck(name+':curvature_response',N.equal(mm(curvature,dq0),mm(L,y)));
 ck(name+':force',N.equal(V,N.scale(R,beta.mul(gap).mul(c).mul(s).neg())));
 ck(name+':tangent',N.equal(N.add(mm(Pc,dP),mm(dP,Pc)),dP));
 ck(name+':gradient_dissipation',real(tr(mm(Hp,dP))).eq(N.energy(V).neg()));
 const score=real(mm(N.dagger(y),K,y)[0][0]),rprime=read(dP),scorePrime=real(mm(N.dagger(r),Kinv,rprime)[0][0]).mul(2).mul(chi.pow(2));
 ck(name+':bending_score_monotone',f(0).le(scorePrime));
 const qpert={...qi};interiorEdges.forEach((e,i)=>qpert[e]=qi[e].add(real(dq0[i][0])));
 ck(name+':small_response_has_valid_cells',cells.every(cell=>positiveGram(gramFromLengths(cell,qpert))));
 const G=N.cayley(N.scale(X,N.IOTA),'1/3'),Pg=mm(G,Pc,N.dagger(G)),Zg=mm(G,Z,N.dagger(G)),Xg=mm(G,X,N.dagger(G));
 ck(name+':observer_frame_invariance',real(tr(mm(Pg,Zg))).eq(c)&&real(tr(mm(Pg,Xg))).eq(s));
 sourceCases.push({c:String(c),s:String(s),amount:String(amount),response:y,response_score:score.toString(),score_derivative:scorePrime.toString(),gradient_generator:V});
}
for(const c of['1','0','3/5','5/13']){
 const cases=sourceCases.filter(r=>r.c===c),base=cases.find(r=>r.amount==='1');
 for(const r of cases)ck('strength:linear_response:'+c+':'+r.amount,N.equal(r.response,N.scale(base.response,r.amount))&&f(r.response_score).eq(f(base.response_score).mul(f(r.amount).pow(2))));
}
const sameN=sourceCases.filter(r=>r.amount==='1');no('strength:amount_alone_does_not_fix_bending',sameN[0].response_score!==sameN[1].response_score);
no('flow:larger_scalar_offset_does_not_change_force',N.equal(N.commutator(P0,Z),N.commutator(P0,N.add(Z,N.scale(I2,1000)))));
ck('alignment:exact_rate',gap.mul(2).eq('37/235'));
// Rational snapshots of the exact analytic relaxation, with decay factor 5/17.
const initial=cut('7/25','24/25'),later=cut('119/169','120/169'),Ui=N.add(N.scale(I2,'63/65'),N.scale(R,'-16/65')),Plift=mm(Ui,initial,N.dagger(Ui));
ck('exponential:exact_Cayley_cut_lift',N.equal(mm(Plift,Plift),Plift)&&real(tr(Plift)).eq(1)&&N.equal(Ui,N.cayley(R,'-1/8')));
ck('exponential:exact_unitary_transport',N.equal(mm(N.dagger(Ui),Ui),I2));
ck('relaxation:exact_two_snapshot_transport',N.equal(Plift,later));
const ri=read(initial),rl=read(later),mi=real(ri[1][0].div(ri[0][0])),ml=real(rl[1][0].div(rl[0][0]));
ck('relaxation:exact_exponential_decay_factor',ml.div(mi).eq('5/17'));
ck('relaxation:strict_response_increase',positive(real(mm(N.dagger(rl),Kinv,rl)[0][0]).sub(real(mm(N.dagger(ri),Kinv,ri)[0][0]))));
for(const Pc of[initial,later]){
 const r=read(Pc),Hp=N.scale(N.add(N.scale(Z,real(r[0][0]).div(10)),N.scale(X,real(r[1][0]).div(47))),-1),dP=N.commutator(N.commutator(Pc,Hp),Pc),rp=read(dP);
 const mdot=real(rp[1][0].mul(r[0][0]).sub(r[1][0].mul(rp[0][0])).div(r[0][0].mul(r[0][0])));
 ck('relaxation:ratio_derivative:'+real(r[0][0]),mdot.eq(real(r[1][0].div(r[0][0])).mul(gap).mul(-2)));
}
// Different continuous generators can give the same sampled return; written phase proof provides 2*pi lift.
no('exponential:sampled_identity_does_not_fix_generator',N.equal(N.power(R,4),I2),{interpretation:'exp(0)=exp(2*pi*R)=I in the completed algebra, while the intervening cut paths differ.'});
const va=N.cayley(R,'1/2'),vb=N.cayley(N.scale(Z,N.IOTA),'1/3');
no('exponential:unitary_transports_need_not_be_flat',!N.equal(mm(va,vb,N.dagger(va),N.dagger(vb)),I2));
ck('regularity:finite_response_bound',sourceCases.every(r=>f(r.response_score).le(alpha.mul(r.amount).div(gamma).pow(2).div(10))));
no('singularity:finite_cut_concentration_is_not_a_blowup',N.energy(P0).eq(1)&&real(tr(P0)).eq(1));
results.coupling={alpha:alpha.toString(),gamma:gamma.toString(),cases:sourceCases,response:'y=(alpha*N/gamma)*diag(1/10,1/47)*r(P)',score:'B^2=y^T*K*y=(alpha*N/gamma)^2*(r_1^2/10+r_2^2/47)',scope:'Exactly solved quadratic metric model. Full nonlinear Regge response retains higher-order corrections.'};
results.exponential={generator:'V=[P,H_eff], H_eff=-(alpha^2*N^2/gamma)*(r_1*Z/10+r_2*X/47)',ratio_law:'m(lambda)=m(0)*exp(-(37/235)*(alpha^2*N^2/gamma)*lambda), m=r_2/r_1, r_1 nonzero',cut_transport:'P(lambda)=U(lambda)*P(0)*U(lambda)^dagger',U:'exp((phi(lambda)-phi(0))*R) in this one-generator real sector; ordered exponential generally',monotone:'d B^2/dlambda=4*(alpha*N/gamma)^2*(alpha^2*N^2/gamma)*(37/470)^2*r_1^2*r_2^2 >= 0',model_choice:'First-order gradient relaxation. Not the earlier conservative CID-D2 evolution or physical clock time.',finite_bound:'B^2 <= (alpha*N/gamma)^2/10; no finite-parameter singularity is produced in this regular finite model.'};
const chosen=proofs.find(p=>p.certificate.steps.length),sys=systems.get(chosen.certificate.presentation),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['997','0'];
no('integrity:corrupt_witness_rejected',(()=>{try{sys.replay(bad);return false;}catch{return true;}})());
const payload={protocol:'EGR1_NATIVE_EXPONENTIAL_GRAVITY_RESPONSE_V1',status:'PASS_LINEARIZED_GEOMETRY_AND_SELECTED_CUT_FLOW',scope:'Exact native linearized metric geometry, gauge null directions, polynomial witnesses and finite source checks. Analytic exponential and relaxation proofs retain THEOREM.md hypotheses; not full nonlinear or Lorentzian quantum gravity.',source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,exact_checks:checks.length,negative_controls:controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',h=hash(bytes),file=path.join(root,'CERTIFICATE.json'),expected=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(file,bytes);fs.writeFileSync(expected,h+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(expected,'utf8'),h+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,mode:write?'generation':'read-only replay',...payload.summary,sha256:h,geometry:{cells:cells.length,interior_edges:interiorEdges.length,gauge_rank:N.rank(gauge),curvature_rank:N.rank(curvature),K},response:results.coupling.response,exponential_rate:results.exponential.ratio_law},null,2));
