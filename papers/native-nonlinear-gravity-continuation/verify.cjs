'use strict';
// EGR-2: native arithmetic and compiler unchanged; finite geometry adapters only.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const root=__dirname,args=process.argv.slice(2),nativeRoot=args[args.indexOf('--native-root')+1],write=args.includes('--write');
assert.ok(args.includes('--native-root')&&nativeRoot);assert.ok(write!==args.includes('--check'));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const[n,h]of Object.entries(pins.native_engine))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,n))),h,n);
for(const[n,h]of Object.entries(pins.packet_files))assert.equal(hash(fs.readFileSync(path.join(root,n))),h,n);
for(const[n,h]of Object.entries(pins.predecessor_files))assert.equal(hash(fs.readFileSync(path.resolve(root,'../..',n))),h,n);
const N=require(path.join(nativeRoot,'native_operator.cjs')),M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const inherited=JSON.parse(fs.readFileSync(path.resolve(root,'../native-exponential-gravity-response/CERTIFICATE.json'),'utf8'));
const data=inherited.results.geometry,f=N.F.of,C=N.Cut.of,mm=(...xs)=>xs.reduce(N.mul);
const real=x=>{assert.ok(x.turn.zero());return x.rad;},col=v=>N.matrix(v.map(x=>[x]));
const columns=vs=>N.matrix(vs[0].map((_,i)=>vs.map(c=>c[i][0]))),trace=g=>g.reduce((s,row,i)=>s.add(row[i]),N.ZERO);
const checks=[],controls=[],presentations=[],proofs=[],systems=new Map(),results={};
const ck=(name,passed,detail={})=>{assert.ok(passed,name);checks.push({name,passed:true,...detail});};
const no=(name,rejected,detail={})=>{assert.ok(rejected,name);controls.push({name,rejected:true,...detail});};
const positive=x=>!x.le(0);
function det(g){return g.length===1?g[0][0]:g[0].reduce((s,z,j)=>s.add(z.mul(det(g.slice(1).map(r=>r.filter((_,k)=>j!==k)))).mul(j%2?-1:1)),N.ZERO);}
const positiveGram=g=>g.every((_,k)=>positive(real(det(g.slice(0,k+1).map(r=>r.slice(0,k+1))))));
const cells=data.cells,edges=data.interior_edges,hinges=data.interior_hinges,q=Object.fromEntries(Object.entries(data.edge_squares).map(([k,v])=>[k,f(v)]));
const H=N.matrix(data.hessian),A=N.matrix(data.mode_readout),S=N.matrix(data.section),G=N.matrix(data.gauge),K=N.matrix(data.reduced_stiffness);
const key=v=>v.slice().sort().join(':'),qof=(a,b)=>a===b?f(0):q[key([a,b])],edgeIndex=Object.fromEntries(edges.map((k,i)=>[k,i]));
const sCols=[0,1].map(j=>S.map(row=>[row[j]])),gCols=Array.from({length:8},(_,j)=>G.map(row=>[row[j]]));
const B10=columns([...sCols,...gCols]);
ck('predecessor:passing_packet',inherited.status==='PASS_LINEARIZED_GEOMETRY_AND_SELECTED_CUT_FLOW');
ck('predecessor:complete_edge_coordinates',N.rank(B10)===10);
ck('predecessor:metric_factorization',N.equal(H,mm(N.dagger(A),K,A)));
ck('predecessor:section_and_null_directions',N.equal(mm(A,S),N.identity(2))&&N.isZero(mm(A,G)));
ck('predecessor:stiffness',N.equal(K,N.matrix([[10,0],[0,47]])));

function affine(g){const B=N.zeros(5);for(let i=1;i<5;i++)for(let j=1;j<5;j++)B[i][j]=g[i-1][j-1];for(let i=1;i<5;i++){B[i][0]=g[i-1].reduce((s,z)=>s.add(z),N.ZERO).neg();B[0][i]=B[i][0];}B[0][0]=g.flat().reduce((s,z)=>s.add(z),N.ZERO);return B;}
const gram=(cell,qs)=>N.matrix(cell.slice(1).map(a=>cell.slice(1).map(b=>qs[key([cell[0],a])].add(qs[key([cell[0],b])]).sub(a===b?0:qs[key([a,b])]).div(2))));
function area2(h){const[a,b,c]=h.split(':'),x=qof(a,b),y=qof(a,c),z=qof(b,c);return x.mul(y).sub(x.add(y).sub(z).pow(2).div(4)).div(4);}
const triangleEdges=h=>{const[a,b,c]=h.split(':');return[key([a,b]),key([a,c]),key([b,c])];};
function si(h,e){const es=triangleEdges(h);return es.includes(e)?es.filter(k=>k!==e).reduce((s,k)=>s.add(q[k]),f(0)).sub(q[e]).div(8):f(0);}

// Directional second derivative of the complete ten-edge Regge gradient.
// At q0 all deficits vanish. The identities used are derived in THEOREM.md.
function directional(d,label){
 const dv=k=>edgeIndex[k]===undefined?f(0):real(d[edgeIndex[k]][0]),paird=(a,b)=>a===b?f(0):dv(key([a,b]));
 const cd=hinges.map(()=>f(0)),cp=hinges.map(()=>f(0)),volumeRelative=[];
 let inverseChecks=true,schlafliChecks=true,gramChecks=true;
 for(let ci=0;ci<cells.length;ci++){
  const cell=cells[ci],base=cell[0],Gc=N.matrix(data.cell_records[ci].gram),inv=N.inverse(Gc),B=affine(inv),v=f(data.cell_records[ci].volume);
  gramChecks=gramChecks&&N.equal(Gc,gram(cell,q))&&real(det(Gc)).eq(v.pow(2).mul(576));
  const Gd=N.matrix(cell.slice(1).map(a=>cell.slice(1).map(b=>paird(base,a).add(paird(base,b)).sub(paird(a,b)).div(2))));
  const invp=N.scale(mm(inv,Gd,inv),-1),invpp=N.scale(mm(inv,Gd,inv,Gd,inv),2),Bp=affine(invp),Bpp=affine(invpp),vp=real(trace(mm(inv,Gd))).mul(v).div(2);
  volumeRelative.push(vp.div(v));
  inverseChecks=inverseChecks&&N.isZero(N.add(mm(Gc,invp),mm(Gd,inv)))&&N.isZero(N.add(mm(Gc,invpp),N.scale(mm(Gd,invp),2)));
  let sk=f(0),skp=f(0);
  for(let i=0;i<5;i++)for(let j=i+1;j<5;j++){
   const ell=Bp[i][i].div(B[i][i]).add(Bp[j][j].div(B[j][j]));
   const ellp=Bpp[i][i].div(B[i][i]).sub(Bp[i][i].div(B[i][i]).mul(Bp[i][i].div(B[i][i]))).add(Bpp[j][j].div(B[j][j]).sub(Bp[j][j].div(B[j][j]).mul(Bp[j][j].div(B[j][j]))));
   const k=real(Bp[i][j].sub(B[i][j].mul(ell).mul('1/2'))),kp=real(Bpp[i][j].sub(Bp[i][j].mul(ell).mul('1/2')).sub(B[i][j].mul(ellp).mul('1/2')));
   sk=sk.add(k);skp=skp.add(kp);const h=key(cell.filter((_,a)=>a!==i&&a!==j)),hi=hinges.indexOf(h);
   if(hi>=0){cd[hi]=cd[hi].sub(v.mul(12).mul(k));cp[hi]=cp[hi].sub(vp.mul(k).add(v.mul(kp)).mul(12));}
  }
  schlafliChecks=schlafliChecks&&sk.zero()&&skp.zero();
 }
 const first=edges.map(()=>f(0)),second=edges.map(()=>f(0));
 hinges.forEach((h,hi)=>{
  const s=area2(h),es=triangleEdges(h),sd=es.reduce((sum,e)=>sum.add(si(h,e).mul(dv(e))),f(0));
  edges.forEach((e,ei)=>{
   const a=si(h,e),ad=es.includes(e)?es.filter(k=>k!==e).reduce((sum,k)=>sum.add(dv(k)),f(0)).sub(dv(e)).div(8):f(0);
   first[ei]=first[ei].add(a.div(s.mul(2)).mul(cd[hi]));
   second[ei]=second[ei].add(ad.div(s).sub(a.mul(sd).mul(3).div(s.pow(2).mul(4))).mul(cd[hi]).add(a.div(s.mul(2)).mul(cp[hi])));
  });
 });
 ck(label+':gram_and_volume',gramChecks);
 ck(label+':inverse_derivatives',inverseChecks);
 ck(label+':Schlafli_and_its_derivative',schlafliChecks);
 ck(label+':first_derivative_matches_inherited_H',N.equal(col(first),mm(H,d)));
 return{second:col(second),volumeRelative,area_weighted_deficit_first:cd,derivative_of_area_weighted_deficit_first:cp};
}

// Polarization determines the entire third derivative tensor, not only one source.
const basis=Array.from({length:10},(_,j)=>col(edges.map((_,i)=>i===j?1:0)));
const qBasis=basis.map((b,j)=>directional(b,'tensor:basis_'+j));
const T=Array.from({length:10},()=>Array.from({length:10},()=>Array(10).fill(N.ZERO)));
for(let j=0;j<10;j++)for(let i=0;i<10;i++)T[i][j][j]=qBasis[j].second[i][0];
for(let j=0;j<10;j++)for(let k=j+1;k<10;k++){
 const z=N.scale(N.sub(N.sub(directional(N.add(basis[j],basis[k]),'tensor:pair_'+j+'_'+k).second,qBasis[j].second),qBasis[k].second),'1/2');
 for(let i=0;i<10;i++)T[i][j][k]=T[i][k][j]=z[i][0];
}
ck('tensor:all_third_derivative_permutations',T.every((plane,i)=>plane.every((row,j)=>row.every((v,k)=>v.eq(T[j][i][k])&&v.eq(T[k][j][i])))));
function contract(a,b){return col(T.map(plane=>plane.reduce((sum,row,j)=>sum.add(row.reduce((s,z,k)=>s.add(z.mul(a[j][0]).mul(b[k][0])),N.ZERO)),N.ZERO)));}
function jacobianNull(d){return mm(N.dagger(G),columns(gCols.map(g=>contract(d,g))));}
let flatCubicZero=true;
for(let i=0;i<8;i++)for(let j=i;j<8;j++)flatCubicZero=flatCubicZero&&N.isZero(mm(N.dagger(G),contract(gCols[i],gCols[j])));
ck('flat_family:all_pure_null_cubic_contractions_vanish',flatCubicZero);

const shapeM=(r1,r2)=>N.matrix(Array.from({length:8},(_,i)=>Array.from({length:8},(_,j)=>{
 if(i%4!==j%4)return 0;
 return i<4?(j<4?f(r1).mul(-10).sub(f(r2).mul('2/5')):f(r2).mul(2)):(j<4?f(r2).mul(2):f(r2).mul(-10));
})));
for(const r of[[1,0],[0,1]]){
 const d0=mm(S,N.inverse(K),col(r));
 ck('null_lift:coefficient_'+r.join('_'),N.equal(jacobianNull(d0),shapeM(...r)));
}
results.geometry={predecessor:'EGR-1',cells,interior_edges:edges,interior_hinges:hinges,edge_squares:q,mode_readout:A,section:S,gauge:G,stiffness:K,third_derivative_tensor:T,tensor_index_order:'T[i][j][k] = partial_i partial_j partial_k R at q0; edges in interior_edges order'};
results.general={epsilon:'alpha*N/gamma',action:'R(q) - epsilon*r^T*A*(q-q0)',source_contract:'Fixed linear functional of squared lengths; A is only a linearly gauge-invariant readout at q0.',null_lift:'M(r)=[[-10*r1-2*r2/5,2*r2],[2*r2,-10*r2]] tensor I4',block_determinant:'100*r1*r2',existence:'For fixed r with r1*r2 != 0: a locally unique analytic stationary branch in the scaled-coordinate neighborhood specified in THEOREM.md.',analytic_scope:'Written analytic implicit-function proof; the native certificate checks all finite algebraic hypotheses, not a proof-assistant formalization or a finite convergence radius.'};

const responses=[];
for(const sign of[1,-1]){
 const r=N.scale(col(['3/5','4/5']),sign),z0=mm(N.inverse(K),r),d0=mm(S,z0),qNaive=directional(d0,'source_'+sign+':naive');
 const b=mm(N.dagger(G),qNaive.second),D=N.scale(jacobianNull(d0),2),eta=N.scale(mm(N.inverse(D),b),-1),d=N.add(d0,mm(G,eta));
 const qd=directional(d,'source_'+sign+':corrected'),z2=N.scale(mm(N.inverse(K),N.dagger(S),qd.second),'-1/2'),e=N.mul(S,z2);
 const nullLift=jacobianNull(d),pairBlock=N.matrix([[nullLift[0][0],nullLift[0][4]],[nullLift[4][0],nullLift[4][4]]]);
 ck('source_'+sign+':tensor_matches_direct_geometry',N.equal(contract(d,d),qd.second));
 ck('source_'+sign+':first_order_equations',N.equal(mm(H,d),mm(N.dagger(A),r)));
 ck('source_'+sign+':eight_second_order_compatibility_equations',N.isZero(mm(N.dagger(G),qd.second)));
 ck('source_'+sign+':all_ten_second_order_equations',N.isZero(N.add(mm(H,e),N.scale(qd.second,'1/2'))));
 ck('source_'+sign+':coordinate_projection',N.equal(mm(A,d),z0)&&N.equal(mm(A,e),z2));
 ck('source_'+sign+':consistent_null_selection',N.isZero(N.add(b,mm(D,eta))));
 ck('source_'+sign+':null_lift_unchanged_by_eta',N.equal(nullLift,N.scale(D,'1/2')));
 ck('source_'+sign+':null_lift_has_four_equal_blocks',N.equal(nullLift,shapeM(sign===1?'3/5':'-3/5',sign===1?'4/5':'-4/5')));
 ck('source_'+sign+':block_determinant',real(det(pairBlock)).eq(48));
 ck('source_'+sign+':definiteness',positiveGram(N.scale(pairBlock,-sign)));
 ck('source_'+sign+':all_volume_derivatives_have_source_sign',qd.volumeRelative.every(v=>positive(v.mul(sign))));
 ck('source_'+sign+':all_axis_displacements_are_permutation_symmetric',eta.slice(0,4).every(row=>row[0].eq(eta[0][0]))&&eta.slice(4).every(row=>row[0].eq(eta[4][0])));
 const jdet=real(det(K)).mul(real(det(pairBlock)).pow(4));
 ck('source_'+sign+':IFT_jacobian_nonzero',jdet.eq('2494955520'));
 const scalarCubic=real(mm(N.dagger(d),qd.second)[0][0]),rZ2=real(mm(N.dagger(r),z2)[0][0]);
 ck('source_'+sign+':on_shell_cubic_consistency',scalarCubic.eq(rZ2.mul(-2)));
 const baseTotal=data.cell_records.reduce((sum,c)=>sum.add(c.volume),f(0)),totalPrime=data.cell_records.reduce((sum,c,i)=>sum.add(f(c.volume).mul(qd.volumeRelative[i])),f(0));
 no('source_'+sign+':freezing_flat_null_directions_fails_nonlinearly',!N.isZero(b));
 const naiveSecond=N.scale(mm(S,N.inverse(K),N.dagger(S),qNaive.second),'-1/2');
 no('source_'+sign+':two_selected_equations_leave_full_residual',!N.isZero(N.add(mm(H,naiveSecond),N.scale(qNaive.second,'1/2'))));
 responses.push({sign,readout:r,linear_readout_coefficient:z0,naive_edge_derivative:d0,naive_null_residual:b,desingularized_null_jacobian:D,selected_flat_displacement:eta,edge_derivative:d,gradient_second_derivative:qd.second,readout_second_coefficient:z2,one_second_order_edge_jet:e,jet_scope:'The actual branch has q=q0+epsilon*d+epsilon^2*(e+G*eta1)+O(epsilon^3). eta1 is selected at the next order and is not calculated here; A*G=0 makes the displayed readout coefficient exact.',null_hessian_first_coefficient:nullLift,null_block_determinant:'48',IFT_full_jacobian_determinant:jdet.toString(),volume_log_derivative_at_zero:qd.volumeRelative,total_volume_log_derivative:totalPrime.div(baseTotal).toString(),action_inertia_for_small_positive_epsilon:sign===1?{positive:2,negative:8,zero:0}:{positive:10,negative:0,zero:0},on_shell_action:{epsilon_squared:real(mm(N.dagger(r),z0)[0][0]).mul('-1/2').toString(),epsilon_cubed:rZ2.div(-3).toString()},area_weighted_deficit_first:qd.area_weighted_deficit_first});
}
ck('complement:opposite_first_edge_response',N.equal(responses[1].edge_derivative,N.scale(responses[0].edge_derivative,-1)));
ck('complement:equal_second_readout_coefficient',N.equal(responses[1].readout_second_coefficient,responses[0].readout_second_coefficient));
const cutPlus=N.matrix([['4/5','2/5'],['2/5','1/5']]),cutMinus=N.sub(N.identity(2),cutPlus),rotation=N.matrix([[0,-1],[1,0]]);
ck('complement:both_cuts_are_rank_one',N.equal(mm(cutPlus,cutPlus),cutPlus)&&N.equal(mm(cutMinus,cutMinus),cutMinus)&&real(trace(cutPlus)).eq(1)&&real(trace(cutMinus)).eq(1));
ck('complement:unitary_frame_transport',N.equal(mm(rotation,cutPlus,N.dagger(rotation)),cutMinus)&&N.equal(mm(N.dagger(rotation),rotation),N.identity(2)));
no('source:same_strength_and_spectrum_do_not_fix_stability',responses[0].action_inertia_for_small_positive_epsilon.negative!==responses[1].action_inertia_for_small_positive_epsilon.negative,{invariants:'trace(P)=trace(P^2)=1 for both; Z and X are held fixed, so this is a source-orientation change, not a joint observer-frame change.'});
no('source:positive_quadrant_is_not_a_metric_minimum',positiveGram(N.scale(shapeM('3/5','4/5'),-1)));
no('axes:general_IFT_does_not_cover_channel_axes',N.rank(shapeM(1,0))===4&&N.rank(shapeM(0,1))===4);
results.response={plus:responses[0],minus:responses[1],cut_plus:cutPlus,cut_minus:cutMinus,cut_transport:rotation};

// Envelope derivative of the stationary branch gives the nonlinear cut feedback.
const Z=N.matrix([[1,0],[0,-1]]),X=N.matrix([[0,1],[1,0]]),feedback=[];
for(const res of responses){
 const Pc=res.sign===1?cutPlus:cutMinus,r1=real(res.readout[0][0]),r2=real(res.readout[1][0]),coefficients=[];
 for(const[n,z]of[[1,res.linear_readout_coefficient],[2,res.readout_second_coefficient]]){
  const Hz=N.scale(N.add(N.scale(Z,real(z[0][0])),N.scale(X,real(z[1][0]))),-1),V=N.commutator(Pc,Hz),omega=r1.mul(real(z[1][0])).sub(r2.mul(real(z[0][0]))),dP=N.commutator(V,Pc);
  ck('feedback:'+res.sign+':coefficient_'+n+':generator',N.equal(V,N.scale(rotation,omega)));
  ck('feedback:'+res.sign+':coefficient_'+n+':circle_tangent',N.equal(dP,N.add(N.scale(Z,omega.mul(r2).neg()),N.scale(X,omega.mul(r1)))));
  ck('feedback:'+res.sign+':coefficient_'+n+':dissipation_identity',real(trace(mm(Hz,dP))).eq(N.energy(V).neg()));
  coefficients.push({order:n,omega_over_alphaN:omega.toString(),ratio_derivative_over_alphaN:omega.mul(2).div(r1.pow(2)).toString()});
 }
 ck('feedback:'+res.sign+':leading_exponential_rate',f(coefficients[0].ratio_derivative_over_alphaN).eq(r2.div(r1).mul('-37/235')));
 no('feedback:'+res.sign+':pure_exponential_ratio_law_is_not_exact_nonlinearly',!f(coefficients[1].ratio_derivative_over_alphaN).zero());
 feedback.push({source_sign:res.sign,ratio:r2.div(r1).toString(),coefficients});
}
results.feedback={stationary_envelope:'W(P)=gamma*(R(q(P))-epsilon*r(P)^T*A*(q(P)-q0))',readout:'y(P)=A*(q(P)-q0)',effective_cut_force:'H_eff=-alpha*N*(y1*Z+y2*X)',selected_relaxation:"P'=[[P,H_eff],P]; U'=[P,H_eff]*U",real_generator:'V=omega*R, omega=alpha*N*(r1*y2-r2*y1)',exact_ratio_equation:'m_prime=2*alpha*N*(y2-m*y1)/r1, r1 != 0',transport:'U=exp((phi(lambda)-phi(0))*R), r=(cos(2*phi),sin(2*phi)); exact while the local stationary branch exists',dissipation:'W_prime=-||[P,H_eff]||_F^2',local_coefficients:feedback,scope:'Locally defined, explicitly chosen gradient relaxation. No physical time or global continuation is inferred. The simple scalar exponential ratio law is the leading metric approximation.'};

function system(name,tokens){const rules=[];for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Declared commuting real geometric scalars'});const S=new M.Presentation({tokens,rules}),w=Object.fromEntries(tokens.map(t=>[t,S.word([t])])),audit=S.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});systems.set(S.hash(),S);return{S,w};}
function prove(name,lhs,rhs){const S=lhs.system,out=S.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',S.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
function pdet(g,S){return g.length===1?g[0][0]:g[0].reduce((sum,z,j)=>sum.plus(z.times(pdet(g.slice(1).map(r=>r.filter((_,k)=>j!==k)),S)).scale(j%2?-1:1)),S.zero());}
{
 const{S,w:{x,u,v,w,dx,du,dv,dw}}=system('symmetric_cell_geometry',['x','u','v','w','dx','du','dv','dw']);
 const one=S.one(),zero=S.zero(),a=x.plus(one).minus(u).scale('1/2'),b=w.plus(u).minus(v).scale('1/2');
 const outer=[[x,a,a,a],[a,one,zero,zero],[a,zero,one,zero],[a,zero,zero,one]];
 const cap=Array.from({length:4},(_,i)=>Array.from({length:4},(_,j)=>i===j?v:v.minus(one)));
 const inner=[[w,b,b,b],[b,u,u.minus(one),u.minus(one)],[b,u.minus(one),u,u.minus(one)],[b,u.minus(one),u.minus(one),u]];
 const dout=x.minus(a.times(a).scale(3)),dcap=v.scale(4).minus(one.scale(3)),din=u.scale(3).minus(one.scale(2)).times(w).minus(b.times(b).scale(3));
 prove('volume:outer_Gram_determinant',pdet(outer,S),dout);
 prove('volume:cap_Gram_determinant',pdet(cap,S),dcap);
 prove('volume:inner_Gram_determinant',pdet(inner,S),din);
 const deriv={x:dx,u:du,v:dv,w:dw,dx:zero,du:zero,dv:zero,dw:zero};
 prove('volume:outer_determinant_derivative',M.polynomialDerivative(dout,deriv),dx.minus(a.times(dx.minus(du)).scale(3)));
 prove('volume:cap_determinant_derivative',M.polynomialDerivative(dcap,deriv),dv.scale(4));
 prove('volume:inner_determinant_derivative',M.polynomialDerivative(din,deriv),du.times(w).scale(3).plus(u.scale(3).minus(one.scale(2)).times(dw)).minus(b.times(dw.plus(du).minus(dv)).scale(3)));
}
{
 const{S,w:{r,s,o,p}}=system('source_null_lift',['r','s','o','p']);
 const m11=r.scale(-10).minus(s.scale('2/5')),m12=s.scale(2),m22=s.scale(-10);
 prove('null:general_two_by_two_determinant',m11.times(m22).minus(m12.times(m12)),r.times(s).scale(100));
 const flatReadout1=o.times(o).scale(5),flatReadout2=o.times(o).scale('1/5').plus(p.times(p).scale(5)).minus(o.times(p).scale(2));
 const lift=m11.times(o).times(o).plus(m12.times(o).times(p).scale(2)).plus(m22.times(p).times(p));
 prove('null:flat_family_source_breaks_nonlinear_shift_invariance',lift,r.times(flatReadout1).plus(s.times(flatReadout2)).scale(-2));
}
{
 const{S,w:{s,si,sd,sid,c,cp}}=system('gradient_second_derivative',['s','si','sd','sid','c','cp']);
 // Common denominator 4*s^2: derivative of the gradient, including both area derivatives.
 const fromProduct=s.times(sid).times(c).scale(4).minus(si.times(sd).times(c).scale(2)).plus(s.times(si).times(cp).scale(2)).minus(si.times(sd).times(c));
 const formula=s.times(sid).times(c).scale(4).minus(si.times(sd).times(c).scale(3)).plus(s.times(si).times(cp).scale(2));
 prove('derivative:cleared_denominator_second_gradient',fromProduct,formula);
}

function determinants(x,u,v,w){const a=x.add(1).sub(u).div(2),b=w.add(u).sub(v).div(2);return{outer:x.sub(a.pow(2).mul(3)),cap:v.mul(4).sub(3),inner:u.mul(3).sub(2).mul(w).sub(b.pow(2).mul(3))};}
const baseSym={x:q['A:O'],u:q['B:O'],v:q['B:P'],w:q['O:P']},delta=determinants(baseSym.x,baseSym.u,baseSym.v,baseSym.w);
ck('volume:base_admissibility',positive(baseSym.u.sub('2/3'))&&Object.values(delta).every(positive));
ck('volume:base_determinants',delta.outer.eq('1/25')&&delta.cap.eq('1/625')&&delta.inner.eq('1/625'));
const rateByDet=[],dPlus=responses[0].edge_derivative,dx=real(dPlus[0][0]),du=real(dPlus[1][0]),dv=real(dPlus[2][0]),dw=real(dPlus[9][0]);
const a0=baseSym.x.add(1).sub(baseSym.u).div(2),b0=baseSym.w.add(baseSym.u).sub(baseSym.v).div(2);
rateByDet.push(dx.sub(a0.mul(dx.sub(du)).mul(3)).div(delta.outer.mul(2)));
rateByDet.push(dv.mul(4).div(delta.cap.mul(2)));
rateByDet.push(du.mul(baseSym.w).mul(3).add(baseSym.u.mul(3).sub(2).mul(dw)).sub(b0.mul(dw.add(du).sub(dv)).mul(3)).div(delta.inner.mul(2)));
ck('volume:independent_determinant_rates_match_Gram_trace',rateByDet[0].eq(responses[0].volume_log_derivative_at_zero[0])&&rateByDet[1].eq(responses[0].volume_log_derivative_at_zero[4])&&rateByDet[2].eq(responses[0].volume_log_derivative_at_zero[5]));
ck('volume:published_rates',rateByDet[0].eq('163103/848256')&&rateByDet[1].eq('4410869/1060320')&&rateByDet[2].eq('912209/1060320'));
// At positive amount and fixed receiving weights, I_sigma=n_sigma/v_sigma.
const fixedAmount=f(3),baseVolume=f('1/600'),density2=fixedAmount.pow(2).div(baseVolume.pow(2));
ck('density:determinant_relation',density2.mul(delta.cap).eq(fixedAmount.pow(2).mul(576)));
const collapsed={...q};for(const e of edges)if(e.includes('P')&&e!=='O:P')collapsed[e]=f('3/4');
no('volume:positive_edge_squares_alone_do_not_certify_cells',Object.values(collapsed).every(positive)&&!positiveGram(gram(cells[4],collapsed)));
no('density:initial_contraction_is_not_already_a_degenerate_cell',rateByDet.every(positive)&&Object.values(delta).every(positive),{scope:'Exact data at epsilon=0 establish initial derivative signs only. No finite-epsilon zero, global continuation interval or density blowup is certified.'});
results.volume={symmetric_variables:{x:'q_AO',u:'q_OB=q_OC=q_OD=q_OE',v:'q_PB=q_PC=q_PD=q_PE',w:'q_OP'},base:baseSym,determinants_at_base:delta,volume_squared:'det(Gram)/576',admissibility:'x-3*((x+1-u)/2)^2 > 0; v > 3/4; u > 2/3; (3*u-2)*w-3*((w+u-v)/2)^2 > 0',plus_volume_log_derivative_by_class:rateByDet,minus_volume_log_derivative_by_class:rateByDet.map(v=>v.neg()),density_contract:'For separately assigned fixed n_sigma>0, I_sigma=n_sigma/v_sigma. The source coupling does not specify these weights or a spatial transport law.',density_identity:'I_sigma^2*Delta_sigma=576*n_sigma^2',collapse_scope:'Delta -> 0 at fixed positive assigned amount implies density divergence, if such a limit is reached. The branch theorem only gives a sufficiently small regular neighborhood.'};

const chosen=proofs.find(p=>p.certificate.steps.length),sys=systems.get(chosen.certificate.presentation),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['997','0'];
no('integrity:corrupt_symbolic_witness_rejected',(()=>{try{sys.replay(bad);return false;}catch{return true;}})());
const payload={protocol:'EGR2_NATIVE_NONLINEAR_GRAVITY_CONTINUATION_V1',status:'PASS_EXACT_NONLINEAR_JET_AND_LOCAL_BRANCH_HYPOTHESES',scope:'Exact native third derivatives, full ten-edge second-order solvability, desingularized Jacobian, complementary-cut stability and volume identities. Analytic existence and inertia follow the written proof. No certified finite source interval, global collapse, nonlinear gauge symmetry, quantum measure or complete quantum-gravity theory.',source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,exact_checks:checks.length,negative_controls:controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',h=hash(bytes),file=path.join(root,'CERTIFICATE.json'),expected=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(file,bytes);fs.writeFileSync(expected,h+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(file,'utf8'),bytes);assert.equal(fs.readFileSync(expected,'utf8'),h+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,mode:write?'generation':'read-only replay',...payload.summary,sha256:h,branches:responses.map(r=>({readout:r.readout,IFT_jacobian_determinant:r.IFT_full_jacobian_determinant,inertia:r.action_inertia_for_small_positive_epsilon,volume_log_derivatives:r.volume_log_derivative_at_zero}))},null,2));
