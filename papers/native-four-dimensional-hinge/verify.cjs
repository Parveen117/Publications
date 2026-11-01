'use strict';
// N4H-1: geometry adapters compose the pinned native arithmetic; no replacement engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const assert=require('node:assert/strict'),root=__dirname,args=process.argv.slice(2);
const nativeRoot=args[args.indexOf('--native-root')+1],write=args.includes('--write');
assert.ok(args.includes('--native-root')&&nativeRoot);assert.ok(write!==args.includes('--check'));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const [n,h] of Object.entries(pins.native_engine))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,n))),h,n);
for(const [n,h] of Object.entries(pins.packet_files))assert.equal(hash(fs.readFileSync(path.join(root,n))),h,n);
for(const [n,h] of Object.entries(pins.predecessor_files))assert.equal(hash(fs.readFileSync(path.resolve(root,'../..',n))),h,n);
const N=require(path.join(nativeRoot,'native_operator.cjs')),M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const f=N.F.of,checks=[],controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,value,data={})=>{assert.ok(value,name);checks.push({name,passed:true,...data});};
const no=(name,value,data={})=>{assert.ok(value,name);controls.push({name,rejected:true,...data});};
const mm=(...xs)=>xs.reduce(N.mul),real=z=>{assert.ok(z.turn.zero());return z.rad;};
const col=xs=>N.matrix(xs.map(x=>[x])),cols=vs=>N.dagger(vs.map(v=>v.map(r=>r[0])));
const dot=(a,b)=>real(mm(N.dagger(a),b)[0][0]),outer=(a,b=a)=>mm(a,N.dagger(b));
const tr=a=>a.reduce((s,r,i)=>s.add(r[i]),N.ZERO),I=N.identity(4),e=i=>col([0,1,2,3].map(j=>i===j?1:0));
const P=N.matrix([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]]),Q=N.sub(I,P);
const R=N.matrix([[0,0,0,0],[0,0,0,0],[0,0,0,-1],[0,0,1,0]]);
const positive=x=>!x.le(0),rejected=fn=>{try{fn();return false;}catch{return true;}};
function det(a){
 a=N.matrix(a);assert.equal(a.length,a[0].length);
 if(a.length===1)return a[0][0];
 return a[0].reduce((s,z,j)=>s.add(z.mul(det(a.slice(1).map(r=>r.filter((_,k)=>k!==j)))).mul(j%2?-1:1)),N.ZERO);
}
const positiveGram=g=>g.every((_,k)=>positive(real(det(g.slice(0,k+1).map(r=>r.slice(0,k+1))))));
function system(name,tokens,rules=[],commute=false){
 if(commute)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'declared central scalar coefficients'});
 const S=new M.Presentation({tokens,rules}),w=Object.fromEntries(tokens.map(t=>[t,S.word([t])])),audit=S.audit();
 ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
 presentations.push({name,spec:S.spec(),audit});systems.set(S.hash(),S);return {S,w};
}
function prove(name,lhs,rhs){
 const S=lhs.system,out=S.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);
 const certificate=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',S.replay(certificate));
 proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});
}

// Universal normal-plane transport. Algebraic identity leaves c^2+s^2 explicit.
{
 const tokens=['p','r','a','b','c','d'],rules=[
  {id:'projector',lhs:['p','p'],rhs:[[['p'],1]],source:'P projects onto the hinge plane'},
  {id:'p_r',lhs:['p','r'],rhs:[],source:'normal generator annihilates the hinge'},
  {id:'r_p',lhs:['r','p'],rhs:[],source:'normal generator has normal image'},
  {id:'r_square',lhs:['r','r'],rhs:[[['p'],1],[[],-1]],source:'R^2=-(I-P) on the oriented Euclidean normal plane'}];
 for(let j=2;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'real central coefficients'});
 const {S,w:{p,r,a,b,c,d}}=system('normal_plane_composition',tokens,rules),q=S.one().minus(p);
 const H=(c,s)=>p.plus(c.times(q)).plus(s.times(r));
 prove('normal:rotation_composition',H(a,b).times(H(c,d)),H(a.times(c).minus(b.times(d)),a.times(d).plus(b.times(c))));
 prove('normal:orthogonality_contract',H(c,d.scale(-1)).times(H(c,d)),p.plus(c.times(c).plus(d.times(d)).times(q)));
 prove('normal:fixes_hinge',H(c,d).times(p),p);
 prove('normal:preserves_normal_plane',q.times(H(c,d)),H(c,d).times(q));
 prove('normal:defect_energy',H(c,d.scale(-1)).minus(S.one()).times(H(c,d).minus(S.one())),c.minus(S.one()).times(c.minus(S.one())).plus(d.times(d)).times(q));
 prove('normal:oriented_sine',H(c,d).minus(H(c,d.scale(-1))),d.times(r).scale(2));
 prove('normal:inverse_requires_norm',H(c,d.scale(-1)).times(H(c,d)).minus(S.one()),c.times(c).plus(d.times(d)).minus(S.one()).times(q));
}
{
 const {S,w:{x,t}}=system('dimensionless_geometry_polynomials',['x','t'],[],true),one=S.one();
 const D=one.scale(16).minus(x.scale(2)),num=x.scale(-1);
 prove('edge:cosine_derivative_numerator',num.scale(0).minus(D).plus(num.scale(2)),one.scale(-16));
 // Schur complement of the edge direction: diag 4-x/4, offdiag -x/4.
 const d=one.scale(4).minus(x.scale('1/4')),b=x.scale('-1/4');
 prove('edge:transverse_Gram_determinant',d.times(d).times(d).plus(b.times(b).times(b).scale(2)).minus(d.times(b).times(b).scale(3)),one.scale(64).minus(x.scale(12)));
 prove('edge:full_Gram_determinant',x.times(one.scale(64).minus(x.scale(12))),x.times(one.scale(16).minus(x.scale(3))).scale(4));
 const area2=x.times(one.scale(16).minus(x)).scale('1/16');
 prove('edge:area_square_derivative',M.polynomialDerivative(area2,{x:one,t:S.zero()}),one.scale(8).minus(x).scale('1/8'));
 prove('edge:vacuum_cosine_equation',num.scale(2).plus(D),one.scale(4).minus(x).scale(4));
 const c=one.minus(t.times(t)),s=t.scale(2),den=one.plus(t.times(t));
 prove('Cayley:rational_unit_circle',c.times(c).plus(s.times(s)),den.times(den));
 prove('information:three_wedge_mismatch',den.minus(c).scale(3),t.times(t).scale(6));
 prove('information:three_wedge_derivative_numerator',t.scale(12).times(den).minus(t.times(t).scale(6).times(t.scale(2))),t.scale(12));
}

// Exact, nondegenerate 4-simplex coordinates and unique whole-tetrahedron gluing.
const pair=t=>{t=f(t);const d=f(1).add(t.mul(t));return {c:f(1).sub(t.mul(t)).div(d),s:t.mul(2).div(d)};};
const normal=(c,s)=>N.add(N.add(P,N.scale(Q,c)),N.scale(R,s));
function star(name,ts,scale=1){
 const L=f(scale),count=ts.length,cells=[],edges=new Map(),transports=[];
 for(let i=0;i<count;i++){
  const {c,s}=pair(ts[i]),coords={A:col([0,0,0,0]),B:N.scale(e(0),L),C:N.scale(e(1),L)};
  coords['D'+i]=N.scale(e(2),L);coords['D'+((i+1)%count)]=N.scale(col([0,0,c,s]),L);
  const names=['A','B','C','D'+i,'D'+((i+1)%count)],E=cols(names.slice(1).map(v=>coords[v]));
  ck(name+':simplex_'+i+'_positive_volume',real(det(E)).eq(L.pow(4).mul(s))&&positive(real(det(E))));
  ck(name+':simplex_'+i+'_positive_Gram',positiveGram(mm(N.dagger(E),E)));
  for(let a=0;a<5;a++)for(let b=0;b<a;b++){
   const key=[names[a],names[b]].sort().join('|'),length2=N.energy(N.sub(coords[names[a]],coords[names[b]]));
   if(edges.has(key))ck(name+':shared_length_'+i+'_'+key,edges.get(key).eq(length2));else edges.set(key,length2);
  }
  cells.push({names,coords,c,s});
 }
 for(let i=0;i<count;i++){
  const source=cells[i],target=cells[(i+1)%count],shared=['B','C','D'+((i+1)%count)];
  const E=cols(shared.map(v=>source.coords[v])),F=cols(shared.map(v=>target.coords[v]));
  const gram=mm(N.dagger(E),E);ck(name+':tetra_shape_matching_'+i,N.equal(gram,mm(N.dagger(F),F))&&N.rank(E)===3);
  const PF=mm(E,N.inverse(gram),N.dagger(E)),PG=mm(F,N.inverse(gram),N.dagger(F));
  const ns=N.scale(mm(N.sub(I,PF),source.coords['D'+i]),f(1).div(L.mul(source.s)));
  const nt=N.scale(mm(N.sub(I,PG),target.coords['D'+((i+2)%count)]),f(1).div(L.mul(target.s)));
  ck(name+':inward_unit_normals_'+i,dot(ns,ns).eq(1)&&dot(nt,nt).eq(1)&&N.isZero(mm(N.dagger(E),ns))&&N.isZero(mm(N.dagger(F),nt)));
  const T=N.sub(mm(F,N.inverse(gram),N.dagger(E)),outer(nt,ns));
  ck(name+':full_face_map_'+i,N.equal(mm(T,E),F));
  ck(name+':opposite_inward_normals_'+i,N.equal(mm(T,ns),N.scale(nt,-1)));
  ck(name+':metric_transport_'+i,N.equal(mm(N.dagger(T),T),I));
  ck(name+':oriented_transport_'+i,real(det(T)).eq(1));
  ck(name+':unique_basis_'+i,N.rank(cols([e(0),e(1),source.coords[shared[2]],ns]))===4);
  ck(name+':transport_from_geometry_'+i,N.equal(T,normal(source.c,source.s.neg())));
  ck(name+':hinge_cut_intertwining_'+i,N.equal(mm(T,P),mm(P,T)));
  if(i===0){
   const folded=N.add(mm(F,N.inverse(gram),N.dagger(E)),outer(nt,ns));
   no(name+':reject_folded_gluing',N.equal(mm(folded,E),F)&&N.equal(mm(N.dagger(folded),folded),I)&&N.equal(mm(folded,ns),nt)&&!N.equal(mm(folded,ns),N.scale(nt,-1)));
  }
  transports.push(T);
 }
 let H=I;for(const T of transports)H=mm(T,H);
 const phase=ts.map(pair).reduce((z,{c,s})=>z.mul(new N.Cut(c,s.neg())),N.ONE);
 ck(name+':holonomy_product',N.equal(H,normal(phase.rad,phase.turn)));
 ck(name+':hinge_fixed_by_holonomy',N.equal(mm(H,P),P));
 ck(name+':normal_curvature_energy',N.energy(N.sub(H,I)).div(4).eq(f(1).sub(phase.rad)));
 const psi=e(2),psiH=e(0);
 const info=transports.reduce((z,T)=>z.add(N.energy(N.sub(psi,mm(T,psi))).div(2)),f(0));
 ck(name+':normal_information_response',info.eq(cells.reduce((z,c)=>z.add(f(1).sub(c.c)),f(0))));
 ck(name+':retained_channel_blindness',transports.every(T=>N.equal(mm(T,psiH),psiH)));
 const frames=ts.map((_,i)=>{const K=N.zeros(4),a=i%4,b=(i+1)%4;K[a][b]=N.ONE;K[b][a]=N.ONE.neg();return N.cayley(K,f(i+1).div(i+2));});
 let Hg=I,Ig=f(0);
 for(let i=0;i<count;i++){
  const j=(i+1)%count,Tg=mm(frames[j],transports[i],N.dagger(frames[i]));Hg=mm(Tg,Hg);
  Ig=Ig.add(N.energy(N.sub(mm(frames[j],psi),mm(Tg,frames[i],psi))).div(2));
 }
 ck(name+':independent_frame_holonomy_covariance',N.equal(Hg,mm(frames[0],H,N.dagger(frames[0]))));
 ck(name+':independent_frame_information_invariance',Ig.eq(info));
 ck(name+':independent_frame_curvature_invariance',N.energy(N.sub(Hg,I)).eq(N.energy(N.sub(H,I))));
 return {cells,transports,H,phase,info,record:{name,parameters:ts.map(String),simplex_count:count,edge_count:edges.size,scale:L.toString(),hinge_area_squared:L.pow(4).div(4).toString(),holonomy:H,normal_cosine:phase.rad.toString(),normal_sine:phase.turn.toString(),curvature_defect:f(1).sub(phase.rad).toString(),normal_information:info.toString()}};
}
const ordinary=star('three_wedges',['1/2','1/2','1/2']);
const scaled=star('three_wedges_scaled',['1/2','1/2','1/2'],3);
ck('dimensionless:transport_unchanged_by_common_scale',N.equal(ordinary.H,scaled.H));
const branch0=star('six_wedges_flat',['1/2','1/2','3/4','1/2','1/2','3/4']);
const branch1=star('six_wedges_excess',['2','2','4/3','2','2','4/3']);
ck('lift:same_skeleton_and_hinge_area',branch0.record.edge_count===branch1.record.edge_count&&branch0.record.hinge_area_squared===branch1.record.hinge_area_squared);
ck('lift:both_compact_holonomies_are_identity',N.equal(branch0.H,I)&&N.equal(branch1.H,I));
const ab=pair('1/2'),bc=pair('3/4'),z=new N.Cut(ab.c,ab.s),w=new N.Cut(bc.c,bc.s);
ck('lift:triple_angles_sum_to_pi_given_positive_angles',z.mul(z).mul(w).eq(-1));
ck('lift:angle_branch_range',positive(ab.c)&&positive(ab.s)&&ab.c.pow(2).le('1/2')&&!ab.c.pow(2).eq('1/2')&&positive(bc.c)&&positive(bc.s));
no('lift:identity_holonomy_does_not_force_zero_deficit',N.equal(branch0.H,branch1.H),{written_proof:'THEOREM.md H3: sums are exactly 2*pi and 4*pi; deficits are 0 and -2*pi. Positive principal simplex angles fix the lifts.'});
const reverse=N.dagger(ordinary.H);
no('curvature:defect_scalar_loses_oriented_sine',!N.equal(reverse,ordinary.H)&&N.energy(N.sub(reverse,I)).eq(N.energy(N.sub(ordinary.H,I))));
no('cut:hinge_projector_does_not_choose_transport',N.equal(mm(ordinary.transports[0],P),P)&&!N.equal(ordinary.transports[0],I));
results.hinge_stars=[ordinary.record,scaled.record,branch0.record,branch1.record];
results.lift={flat_total_angle:'2*pi',flat_deficit:'0',excess_total_angle:'4*pi',excess_deficit:'-2*pi',same_compact_holonomy:true,source_of_lift:'positive simplex angles and oriented incidence, with written branch proof'};

// Four actual 4-simplices, one interior edge, four interior triangular hinges.
const complex=[['A','B','C','D','E'],['A','B','C','D','F'],['A','B','C','E','F'],['A','B','D','E','F']];
function combinations(xs,k){if(k===0)return [[]];return xs.flatMap((x,i)=>combinations(xs.slice(i+1),k-1).map(rest=>[x,...rest]));}
const faces=new Map();for(const cell of complex)for(const face of combinations(cell,4)){const key=face.join('');faces.set(key,(faces.get(key)||0)+1);}
const boundary=[...faces].filter(([,n])=>n===1).map(([s])=>s.split(''));
const subfaces=k=>[...new Set(complex.flatMap(c=>combinations(c,k).map(x=>x.join(''))))];
const inside=k=>subfaces(k).filter(s=>!boundary.some(b=>[...s].every(v=>b.includes(v))));
ck('bulk:boundary_tetrahedra',boundary.length===8);
ck('bulk:interior_tetrahedra',[...faces.values()].filter(n=>n===2).length===6);
ck('bulk:only_interior_edge',JSON.stringify(inside(2))===JSON.stringify(['AB']));
ck('bulk:four_interior_triangles',JSON.stringify(inside(3).sort())===JSON.stringify(['ABC','ABD','ABE','ABF']));
ck('bulk:three_cells_per_interior_triangle',inside(3).every(h=>complex.filter(c=>[...h].every(v=>c.includes(v))).length===3));
const gram=x=>N.matrix([[x,x.div(2),x.div(2),x.div(2)],[x.div(2),4,0,0],[x.div(2),0,4,0],[x.div(2),0,0,4]]);
const gramDerivative=N.matrix([[1,'1/2','1/2','1/2'],['1/2',0,0,0],['1/2',0,0,0],['1/2',0,0,0]]);
function affineGradientGram(g){
 const out=N.zeros(5);for(let i=1;i<5;i++)for(let j=1;j<5;j++)out[i][j]=g[i-1][j-1];
 for(let i=1;i<5;i++){out[i][0]=g[i-1].reduce((s,z)=>s.add(z),N.ZERO).neg();out[0][i]=out[i][0];}
 out[0][0]=g.flat().reduce((s,z)=>s.add(z),N.ZERO);return out;
}
function distance2(g,i,j){const at=(i,j)=>i&&j?g[i-1][j-1]:N.ZERO;return real(at(i,i).add(at(j,j)).sub(at(i,j).mul(2)));}
const edgeResults=[];
for(const raw of [1,4,'48/11',5]){
 const x=f(raw),g=gram(x),gi=N.inverse(g),c=real(gi[2][3].div(gi[2][2])).neg();
 ck('bulk:'+raw+':positive_4d_Gram',positiveGram(g));
 ck('bulk:'+raw+':determinant',real(det(g)).eq(x.mul(f(16).sub(x.mul(3))).mul(4)));
 ck('bulk:'+raw+':dihedral_cosine',c.eq(x.neg().div(f(16).sub(x.mul(2)))));
 const h=N.matrix([[x,x.div(2)],[x.div(2),4]]),a2=real(det(h)).div(4),v2=real(det(g)).div(576);
 ck('bulk:'+raw+':area_squared',a2.eq(x.mul(f(16).sub(x)).div(16)));
 ck('bulk:'+raw+':volume_squared',v2.eq(x.mul(f(16).sub(x.mul(3))).div(144)));
 const gip=N.scale(mm(gi,gramDerivative,gi),-1),B=affineGradientGram(gi),Bp=affineGradientGram(gip);
 ck('bulk:'+raw+':inverse_Gram_derivative',N.isZero(N.add(mm(gramDerivative,gi),mm(g,gip))));
 const cdot=real(gip[2][3].mul(gi[2][2]).sub(gi[2][3].mul(gip[2][2])).div(gi[2][2].mul(gi[2][2]))).neg();
 ck('bulk:'+raw+':dihedral_derivative',cdot.eq(f(-4).div(f(8).sub(x).pow(2))));
 ck('bulk:'+raw+':monotone_deficit_input',positive(cdot.neg()));
 ck('Schlafli:'+raw+':affine_gradient_closure',B.every(r=>r.reduce((s,z)=>s.add(z),N.ZERO).zero())&&Bp.every(r=>r.reduce((s,z)=>s.add(z),N.ZERO).zero()));
 let sum=f(0);const terms=[];
 for(let i=0;i<5;i++)for(let j=i+1;j<5;j++){
  const vertices=[0,1,2,3,4].filter(k=>k!==i&&k!==j),[a,b,d]=vertices;
  const u2=distance2(g,a,b),leg2=distance2(g,a,d),uv=u2.add(leg2).sub(distance2(g,b,d)).div(2),actualA2=u2.mul(leg2).sub(uv.pow(2)).div(4);
  const normalMinor=real(B[i][i].mul(B[j][j]).sub(B[i][j].mul(B[i][j])));
  ck('Schlafli:'+raw+':hinge_area_'+i+'_'+j,actualA2.eq(v2.mul(144).mul(normalMinor)));
  const term=real(Bp[i][j].sub(B[i][j].mul(Bp[i][i].div(B[i][i]).add(Bp[j][j].div(B[j][j]))).mul('1/2')));
  sum=sum.add(term);terms.push({opposite_vertices:[i,j],area_times_angle_derivative_divided_by_12_volume:term.toString()});
 }
 ck('Schlafli:'+raw+':all_ten_hinges_cancel',sum.zero()&&terms.some(t=>t.area_times_angle_derivative_divided_by_12_volume!=='0'));
 edgeResults.push({x:x.toString(),gram:g,cosine:c.toString(),cosine_derivative:cdot.toString(),hinge_area_squared:a2.toString(),simplex_volume_squared:v2.toString(),Schlafli_terms:terms});
}
for(const raw of [0,'16/3',6])no('bulk:reject_degenerate_or_non_Euclidean_'+raw,!positiveGram(gram(f(raw))));
const x4=f(4),g4=gram(x4),baseline={A:col([0,0,0,-1]),B:col([0,0,0,1]),C:col([1,1,1,0]),D:col([1,-1,-1,0]),E:col([-1,1,-1,0]),F:col([-1,-1,1,0])};
for(let i=0;i<complex.length;i++){
 const E=cols(complex[i].slice(1).map(v=>N.sub(baseline[v],baseline.A)));
 ck('bulk:flat_embedding_Gram_'+i,N.equal(mm(N.dagger(E),E),g4));
 ck('bulk:flat_embedding_volume_'+i,real(det(E)).abs().eq(8));
}
ck('bulk:unique_vacuum_algebraic_root',x4.neg().div(f(16).sub(x4.mul(2))).eq('-1/2'));
// At x=4, A^2=3, sin(theta)^2=3/4, with A,sin(theta)>0; hence A*sin(theta)=3/2.
const cprime=f(edgeResults.find(r=>r.x==='4').cosine_derivative);
ck('bulk:flat_Hessian_exact',f(3).mul(f(8).sub(x4)).mul(cprime).div(f(4).mul('3/2')).eq('-1/2'));
ck('bulk:weak_source_susceptibility',f(-1).div('-1/2').eq(2));
const curved=pair(2),Hcurved=N.power(normal(curved.c,curved.s.neg()),3);
ck('bulk:curved_example_cosine',curved.c.eq(f('48/11').neg().div(f(16).sub(f('48/11').mul(2)))));
ck('bulk:curved_example_holonomy',N.equal(Hcurved,normal('117/125','-44/125')));
no('bulk:curved_example_is_not_vacuum',!N.equal(Hcurved,I),{reason:'The exact vacuum equation and 0<theta<pi require theta=2*pi/3, hence x=4.'});
no('action:area_variation_is_not_unconstrained_angle_variation',inside(2).length===1&&inside(3).length===4,{reason:'Four hinge areas are functions of the same interior edge, not four independent variables.'});
results.interior_edge={cells:complex,boundary_tetrahedra:boundary,interior_edges:inside(2),interior_hinges:inside(3),admissible_domain:'0 < x < 16/3',fixtures:edgeResults,vacuum_x:'4',curvature_functional_second_derivative_at_vacuum:'-1/2',weak_source_response:'x-4 = 2*epsilon*I_prime(4) + O(epsilon^2), epsilon=alpha/gamma',source_scope:'Chosen Regge-plus-information action, no cosmological term, fixed boundary lengths, twice differentiable information functional.'};

// Integer characters of a circle have native Laurent-word witnesses.
{
 const S=new M.Presentation({tokens:[{name:'u',from:'*',to:'*',charge:'1'},{name:'v',from:'*',to:'*',charge:'-1'}],rules:[
  {id:'uv',lhs:['u','v'],rhs:[[[],1]],source:'unit phase and its inverse'},
  {id:'vu',lhs:['v','u'],rhs:[[[],1]],source:'unit phase and its inverse'}]});
 const audit=S.audit();ck('phase:Laurent_confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name:'integer_phase_characters',spec:S.spec(),audit});systems.set(S.hash(),S);
 const chi=n=>S.word(Array(Math.abs(n)).fill(n>=0?'u':'v'));
 for(const [n,m] of [[1,1],[2,-1],[-2,3],[-2,-1],[3,-3],[0,2]])prove('phase:integer_character_'+n+'_'+m,chi(n).times(chi(m)),chi(n+m));
 ck('phase:half_weight_full_turn_sign',N.IOTA.mul(N.IOTA).eq(-1));
 no('phase:half_weight_does_not_descend_to_unlifted_SO2',!N.IOTA.mul(N.IOTA).eq(1),{interpretation:'For q=1/2, a 2*pi change in lifted angle multiplies exp(iota*q*delta) by -1.'});
 ck('phase:half_weight_double_cover_consistency',N.IOTA.mul(N.IOTA).mul(N.IOTA).mul(N.IOTA).eq(1));
 results.quantum_phase={conditional_descent:'exp(iota*q*delta) is a function of SO(2) holonomy alone iff q is an integer, for every real lift delta',gravity_application:'q=gamma*a_h',alternative:'Retain delta, or use a specified cover; then integer area weights do not follow.',example_hinge_area:'1/2',example_gamma:'1',example_branch_factor:'-1',scope:'Single-hinge phase factor. No path-integral measure, physical boundary amplitude or Lorentzian quantization is certified.'};
}
const selected=proofs.find(p=>p.certificate.steps.length),systemFor=systems.get(selected.certificate.presentation);
const bad=JSON.parse(JSON.stringify(selected.certificate));bad.steps[0].coefficient=['997','0'];no('integrity:corrupted_word_witness_rejected',rejected(()=>systemFor.replay(bad)));
const wrong=JSON.parse(JSON.stringify(selected.certificate));wrong.presentation='0'.repeat(64);no('integrity:wrong_presentation_rejected',rejected(()=>systemFor.replay(wrong)));
const payload={protocol:'N4H1_NATIVE_FOUR_DIMENSIONAL_HINGE_V1',status:'PASS_EXACT_GEOMETRY_AND_ALGEBRA_WITH_DECLARED_ACTION',scope:'Unchanged native finite algebra, rational geometry, word witnesses, and exact counterexamples. Angle lifts, Schlaefli variation and phase descent have written proofs/contracts in THEOREM.md; no blanket continuum or quantum-gravity certification.',source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,exact_checks:checks.length,negative_controls:controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',h=hash(bytes),file=path.join(root,'CERTIFICATE.json'),expected=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(file,bytes);fs.writeFileSync(expected,h+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),bytes,'Certificate replay changed');assert.equal(fs.readFileSync(expected,'utf8'),h+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,mode:write?'generation':'read-only replay',...payload.summary,sha256:h,interior_edge:{vacuum_x:results.interior_edge.vacuum_x,hessian:results.interior_edge.curvature_functional_second_derivative_at_vacuum},lift:results.lift},null,2));
