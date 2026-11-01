'use strict';
// Native arithmetic only. No alternative scalar, polynomial or matrix engine.
const fs=require('node:fs'), path=require('node:path'), crypto=require('node:crypto');
const assert=require('node:assert/strict'), root=__dirname, args=process.argv.slice(2);
const arg=k=>{const i=args.indexOf(k);assert.ok(i>=0&&args[i+1],k);return args[i+1];};
const nativeRoot=path.resolve(arg('--native-root')), upstream=JSON.parse(fs.readFileSync(arg('--upstream-report'),'utf8'));
const write=args.includes('--write');assert.ok(write!==args.includes('--check'),'Choose --write or --check');
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const [file,expected] of Object.entries(pins.native_engine.files))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,file))),expected,file);
const N=require(path.join(nativeRoot,'native_operator.cjs')), M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const {F,Cut,IOTA}=N, f=x=>F.of(x), I=N.identity(2), O=N.zeros(2);
const R=N.matrix([[0,-1],[1,0]]), J=N.scale(N.matrix([[0,1],[1,0]]),IOTA), L=N.mul(R,J);
const basis=[R,J,L], mm=(...xs)=>xs.reduce(N.mul), trace=A=>A[0][0].add(A[1][1]);
const real=x=>{assert.ok(x.turn.zero(),'Real native scalar required');return x.rad;};
const sc=A=>real(trace(A)).div(2), pair=(A,B)=>sc(N.mul(A,B)).neg();
const combo=coeffs=>coeffs.reduce((a,c,i)=>N.add(a,N.scale(basis[i],c)),O);
const checks=[],controls=[],certificates=[],presentations=[],systems=new Map();
const check=(name,ok,detail={})=>{assert.ok(ok,name);checks.push({name,passed:true,...detail});};
const control=(name,ok,detail={})=>{assert.ok(ok,name);controls.push({name,false_inference_rejected:true,...detail});};
const rejects=fn=>{try{fn();return false;}catch{return true;}};
function make(name,tokens,{commutative=false,rules=[],star=null}={}){
  const commute=[];if(commutative)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)commute.push({id:tokens[j]+'_past_'+tokens[i],lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Real central native coefficients'});
  const S=new M.Presentation({tokens,rules:commute.concat(rules)}),w=Object.fromEntries(tokens.map(t=>[t,S.word([t])])),audit=S.audit();
  check(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');
  const images=star?star(S,w):commutative?w:null,dagger=images?M.daggerGate(S,images):null;
  if(dagger)check(name+':dagger',dagger.status==='DAGGER_DESCENDS');
  presentations.push({name,spec:S.spec(),audit,dagger,dagger_images:images?Object.fromEntries(Object.entries(images).map(([k,e])=>[k,e.toJSON()])):null});
  systems.set(S.hash(),S);return {S,w};
}
function prove(name,lhs,rhs){
  const S=lhs.system,r=M.proveEquality(lhs,rhs);check(name,r.status==='EQUAL_IN_DECLARED_QUOTIENT');
  const c=JSON.parse(JSON.stringify(r.certificate));assert.ok(S.fromJSON(c.input).equals(lhs.minus(rhs)));assert.ok(S.fromJSON(c.output).equals(S.zero(lhs.from,lhs.to)));
  check(name+':serialized_replay',S.replay(c));certificates.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate:c});
}
const times=(...xs)=>xs.reduce((a,b)=>a.times(b)),sq=x=>x.times(x);
function sym(S){
  const one=S.one(),zero=S.zero(),constant=A=>A.map(row=>row.map(c=>one.scale(c)));
  const add=(A,B)=>A.map((row,i)=>row.map((x,j)=>x.plus(B[i][j]))),sub=(A,B)=>A.map((row,i)=>row.map((x,j)=>x.minus(B[i][j])));
  const scale=(A,x)=>A.map(row=>row.map(a=>a.times(x))),mul=(A,B)=>A.map(row=>B[0].map((_,j)=>row.reduce((z,x,k)=>z.plus(x.times(B[k][j])),zero)));
  const tr=A=>A[0][0].plus(A[1][1]),eq=(name,A,B)=>{for(let i=0;i<2;i++)for(let j=0;j<2;j++)prove(name+':'+i+j,A[i][j],B[i][j]);};
  const sb=basis.map(constant),vector=c=>c.reduce((A,x,i)=>add(A,scale(sb[i],x)),constant(O));
  return {constant,add,sub,scale,mul,tr,eq,vector,id:constant(I),zero:constant(O)};
}

// R1: exact vertex-frame cancellation, with genuine noncommuting words.
{
  const names=['g0','g1','g2','g3','ux','uy','ut','ur','vx','vy','vt','vr'];
  const tokens=names.flatMap(x=>[x,x+'d']),rules=names.flatMap(x=>[
    {id:x+'_unit_right',lhs:[x,x+'d'],rhs:[[[],1]],source:'Native unit word and its dagger inverse'},
    {id:x+'_unit_left',lhs:[x+'d',x],rhs:[[[],1]],source:'Native unit word and its dagger inverse'}]);
  const {S,w}=make('R1_vertex_word_group',tokens,{rules,star:(_,w)=>Object.fromEntries(names.flatMap(x=>[[x,w[x+'d']],[x+'d',w[x]]]))});
  const edge=(name,s,t)=>times(w['g'+t+'d'],w[name],w['g'+s]);
  const reverse=(name,s,t)=>times(w['g'+s+'d'],w[name+'d'],w['g'+t]);
  const loop=p=>times(w[p+'yd'],w[p+'td'],w[p+'r'],w[p+'x']);
  const invLoop=p=>times(w[p+'xd'],w[p+'rd'],w[p+'t'],w[p+'y']);
  const movedLoop=p=>times(reverse(p+'y',0,3),reverse(p+'t',3,2),edge(p+'r',1,2),edge(p+'x',0,1));
  const movedInvLoop=p=>times(reverse(p+'x',0,1),reverse(p+'r',1,2),edge(p+'t',3,2),edge(p+'y',0,3));
  const conjugate=x=>times(w.g0d,x,w.g0);
  prove('R1:source_loop_covariance',movedLoop('u'),conjugate(loop('u')));
  prove('R1:reference_loop_covariance',movedLoop('v'),conjugate(loop('v')));
  prove('R1:relative_x_covariance',reverse('vx',0,1).times(edge('ux',0,1)),conjugate(w.vxd.times(w.ux)));
  prove('R1:relative_y_covariance',reverse('vy',0,3).times(edge('uy',0,3)),conjugate(w.vyd.times(w.uy)));
  prove('R1:relative_loop_covariance',movedInvLoop('v').times(movedLoop('u')),conjugate(invLoop('v').times(loop('u'))));
}
// Algebraic inverse transport: no dagger premise is needed in this presentation.
{
  const {S,w:{g,gi,q,t}}=make('R2_inverse_conjugation',['g','gi','q','t'],{rules:[
    {id:'g_gi',lhs:['g','gi'],rhs:[[[],1]],source:'Declared two-sided inverse'},
    {id:'gi_g',lhs:['gi','g'],rhs:[[[],1]],source:'Declared two-sided inverse'},
    {id:'q_t',lhs:['q','t'],rhs:[[[],1],[['t'],-1]],source:'(1+q)t=1 on the admitted Cayley chart'},
    {id:'t_q',lhs:['t','q'],rhs:[[[],1],[['t'],-1]],source:'t(1+q)=1 on the admitted Cayley chart'}]});
  const one=S.one(),cq=times(gi,q,g),ct=times(gi,t,g);
  prove('R2:inverse_transforms_left',one.plus(cq).times(ct),one);
  prove('R2:inverse_transforms_right',ct.times(one.plus(cq)),one);
  prove('R2:Cayley_readout_covariance',one.minus(cq).times(ct).scale(2),times(gi,one.minus(q),t,g).scale(2));
}
{
  const {S,w:{s,x,y,z}}=make('R2_native_unit_readout',['s','x','y','z'],{commutative:true,rules:[
    {id:'unit_norm',lhs:['z','z'],rhs:[[[],1],[['s','s'],-1],[['x','x'],-1],[['y','y'],-1]],source:'Native real unit word s^2+x^2+y^2+z^2=1'}]});
  const a=sym(S),v=a.vector([x,y,z]),Q=a.add(a.scale(a.id,s),v),Qd=a.sub(a.scale(a.id,s),v);
  a.eq('R2:unit_word',a.mul(Q,Qd),a.id);
  a.eq('R2:inverse_denominator',a.mul(a.add(a.id,Q),a.add(a.id,Qd)),a.scale(a.id,S.one().plus(s).scale(2)));
  a.eq('R2:Cayley_vector_numerator',a.mul(a.sub(a.id,Q),a.add(a.id,Qd)),a.scale(v,S.one().scale(-2)));
  prove('R2:Cayley_norm_numerator',sq(x).plus(sq(y)).plus(sq(z)).scale(4),S.one().minus(s).times(S.one().plus(s)).scale(4));
  prove('R7:return_defect_native_norm',a.tr(a.mul(a.sub(a.id,Qd),a.sub(a.id,Q))).scale('1/8'),S.one().minus(s).scale('1/2'));
}
{
  const {S,w:{a,b,x,y,z,u,v,w}}=make('R2_pairing_and_area',['a','b','x','y','z','u','v','w'],{commutative:true});
  const m=sym(S),X=m.vector([x,y,z]),Y=m.vector([u,v,w]),dot=x.times(u).plus(y.times(v)).plus(z.times(w));
  const nx=sq(x).plus(sq(y)).plus(sq(z)),ny=sq(u).plus(sq(v)).plus(sq(w));
  prove('R2:native_pairing_coefficients',m.tr(m.mul(X,Y)).scale('-1/2'),dot);
  prove('R2:cross_scalar_reading',m.tr(m.mul(m.sub(m.scale(m.id,a),X),m.add(m.scale(m.id,b),Y))).scale('1/2'),a.times(b).plus(dot));
  prove('R2:Gram_area_minor_identity',nx.times(ny).minus(sq(dot)),sq(x.times(v).minus(y.times(u))).plus(sq(x.times(w).minus(z.times(u)))).plus(sq(y.times(w).minus(z.times(v)))));
}
{
  const {S,w:{a,b,c}}=make('R2_four_reading_denominator',['a','b','c'],{commutative:true});
  const nx=S.one().minus(sq(a)).scale(4),ny=S.one().minus(sq(b)).scale(4),cross=c.minus(a.times(b)).scale(4);
  prove('R2:four_reading_area_denominator',nx.times(ny).minus(sq(cross)),S.one().minus(sq(a)).times(S.one().minus(sq(b))).minus(sq(c.minus(a.times(b)))).scale(16));
}

// R3: both spatial derivatives and the inherited noncommuting word order.
{
  const operators=['Ax','Ay','dxAy','dyAx','Bx','Fa','Fb'],tokens=['e','n',...operators];
  const rules=[{id:'e2',lhs:['e','e'],rhs:[],source:'First comparison jet e^2=0'},
    {id:'n2',lhs:['n','n'],rhs:[],source:'First comparison jet n^2=0'},
    {id:'n_e',lhs:['n','e'],rhs:[[['e','n'],1]],source:'Real central comparison increments'}];
  for(const t of operators)for(const d of ['e','n'])rules.push({id:t+'_past_'+d,lhs:[t,d],rhs:[[[d,t],1]],source:'Central first-order comparison increments'});
  const {S,w}=make('R3_variable_plaquette',tokens,{rules,star:(_,w)=>Object.fromEntries(tokens.map(t=>[t,operators.includes(t)?w[t].scale(-1):w[t]]))});
  const {e,n,Ax,Ay,dxAy,dyAx,Bx,Fa,Fb}=w,one=S.one(),area=e.times(n),Fxy=dxAy.minus(dyAx).plus(M.bracket(Ax,Ay));
  const H=times(one.plus(n.times(Ay)),one.plus(e.times(Ax)).plus(area.times(dyAx)),one.minus(n.times(Ay)).minus(area.times(dxAy)),one.minus(e.times(Ax)));
  prove('R3:full_variable_field_word',H,one.minus(area.times(Fxy)));
  const P=one.plus(area.times(Fb)).times(one.minus(area.times(Fa))),relative=Fa.minus(Fb);
  prove('R3:relative_curvature_jet',P,one.minus(area.times(relative)));
  const inverse=one.scale('1/2').plus(area.times(relative).scale('1/4'));
  prove('R3:relative_loop_inverse',one.plus(P).times(inverse),one);
  prove('R3:relative_loop_Cayley_jet',one.minus(P).times(inverse).scale(2),area.times(relative));
  const edge=one.plus(e.times(Bx)).times(one.minus(e.times(Ax))),edgeInverse=one.scale('1/2').plus(e.times(Ax.minus(Bx)).scale('1/4'));
  prove('R3:relative_edge_inverse',one.plus(edge).times(edgeInverse),one);
  prove('R3:relative_edge_Cayley_jet',one.minus(edge).times(edgeInverse).scale(2),e.times(Ax.minus(Bx)));
}
{
  const {w:{Bx,By,Cx,Cy,D}}=make('R3_connection_difference',['Bx','By','Cx','Cy','D']);
  prove('R3:covariant_reference_difference',D.plus(M.bracket(Bx.plus(Cx),By.plus(Cy))).minus(M.bracket(Bx,By)),D.plus(M.bracket(Bx,Cy)).plus(M.bracket(Cx,By)).plus(M.bracket(Cx,Cy)));
}
{
  const {S,w:{ax,ay,bx,by,cx,cy,da,db,dc}}=make('R4_three_channel_curvature',['ax','ay','bx','by','cx','cy','da','db','dc'],{commutative:true});
  const m=sym(S),X=m.vector([ax,bx,cx]),Y=m.vector([ay,by,cy]);
  const cr=da.plus(bx.times(cy).minus(cx.times(by)).scale(2)),cj=db.plus(cx.times(ay).minus(ax.times(cy)).scale(2));
  const area=ax.times(by).minus(bx.times(ay)),cl=dc.plus(area.scale(2));
  const curvature=m.add(m.vector([da,db,dc]),m.sub(m.mul(X,Y),m.mul(Y,X)));
  m.eq('R4:full_three_channel_coefficients',curvature,m.vector([cr,cj,cl]));
  prove('R4:recover_third_channel_curl',m.tr(m.mul(m.constant(L),curvature)).scale('-1/2').minus(area.scale(2)),dc);
}
{
  const {S,w:{ax,ay,bx,by,u,v}}=make('R5_two_current_calibration',['ax','ay','bx','by','u','v'],{commutative:true});
  const m=sym(S),area=ax.times(by).minus(ay.times(bx)),curvature=m.vector([u,v,area.scale(2)]);
  const gx=sq(ax).plus(sq(bx)),gy=sq(ay).plus(sq(by)),gxy=ax.times(ay).plus(bx.times(by));
  prove('R5:two_current_Gram',gx.times(gy).minus(sq(gxy)),sq(area));
  prove('R5:two_current_curvature_norm',m.tr(m.mul(curvature,curvature)).scale('-1/2'),sq(u).plus(sq(v)).plus(sq(area).scale(4)));
}
{
  const {w:{a,b,c,p,q,r,s}}=make('R5_coordinate_area',['a','b','c','p','q','r','s'],{commutative:true});
  const g11=sq(p).times(a).plus(p.times(r).times(b).scale(2)).plus(sq(r).times(c));
  const g22=sq(q).times(a).plus(q.times(s).times(b).scale(2)).plus(sq(s).times(c));
  const g12=p.times(q).times(a).plus(p.times(s).plus(q.times(r)).times(b)).plus(r.times(s).times(c));
  prove('R5:coordinate_Gram_determinant',g11.times(g22).minus(sq(g12)),sq(p.times(s).minus(q.times(r))).times(a.times(c).minus(sq(b))));
}

// R7: the bounded return defect removes only the loop Cayley pole.
{
  const {S,w:{d}}=make('R7_return_defect_complement',['d'],{commutative:true});
  const chi=S.one().minus(d).scale('1/2');
  prove('R7:information_curvature_complement',S.one().minus(chi).times(S.one().minus(d)),chi.times(S.one().plus(d)));
}

// Finite data decoder and native matrix implementation are independent routes.
function nativeUnit(Q){
  assert.ok(N.equal(mm(N.dagger(Q),Q),I),'Unit transport required');
  const coeff=[sc(Q),...basis.map(B=>pair(B,Q))];
  const rebuilt=N.add(N.scale(I,coeff[0]),combo(coeff.slice(1)));
  assert.ok(N.equal(Q,rebuilt),'Real native unit-word sector required');
  assert.ok(coeff.reduce((sum,x)=>sum.add(x.pow(2)),f(0)).eq(1));return coeff;
}
function readout(Q){
  nativeUnit(Q);const D=N.scale(mm(N.sub(I,Q),N.inverse(N.add(I,Q))),2);
  assert.ok(N.equal(N.dagger(D),N.scale(D,-1))&&sc(D).zero(),'Anti-dagger readout required');return D;
}
function scalarDomain({a,b,c,d},{loopAntipode=false}={}){
  [a,b,c,d]=[a,b,c,d].map(f);
  assert.ok(!a.le(-1)&&!f(1).le(a)&&!b.le(-1)&&!f(1).le(b),'Two nontrivial finite edge charts required');
  assert.ok((loopAntipode?f(-1).le(d):!d.le(-1))&&d.le(1)&&f(-1).le(c)&&c.le(1),'Loop/cross scalar outside native unit domain');
  const delta=f(1).sub(a.pow(2)).mul(f(1).sub(b.pow(2))).sub(c.sub(a.mul(b)).pow(2));
  assert.ok(!delta.le(0),'Positive comparison area required');
  return {a,b,c,d,delta};
}
function scalarDecode(readings){
  const {a,b,d,delta}=scalarDomain(readings);
  return f(1).sub(d).mul(f(1).add(a).pow(2)).mul(f(1).add(b).pow(2)).div(f(16).mul(f(1).add(d)).mul(delta));
}
function informationDecode(readings){
  const {a,b,d,delta}=scalarDomain(readings,{loopAntipode:true});
  return f(1).sub(d).mul(f(1).add(a).pow(2)).mul(f(1).add(b).pow(2)).div(delta.mul(32));
}
const identityLinks=()=>({x:I,y:I,t:I,r:I});
const loop=E=>mm(N.dagger(E.y),N.dagger(E.t),E.r,E.x);
function observeInformation(U,V){
  Object.values(U).concat(Object.values(V)).forEach(nativeUnit);
  const Hu=loop(U),Hv=loop(V),Px=mm(N.dagger(V.x),U.x),Py=mm(N.dagger(V.y),U.y),P=mm(N.dagger(Hv),Hu);
  const Dx=readout(Px),Dy=readout(Py),gxx=pair(Dx,Dx),gxy=pair(Dx,Dy),gyy=pair(Dy,Dy),area=gxx.mul(gyy).sub(gxy.pow(2));
  assert.ok(!area.le(0),'Degenerate comparison area');
  const readings={a:sc(Px),b:sc(Py),c:sc(mm(N.dagger(Px),Py)),d:sc(P)};
  const chi=f(1).sub(readings.d).div(2),information=chi.div(area),shape=area.div(gxx.mul(gyy));
  const scalar=scalarDomain(readings,{loopAntipode:true});
  assert.ok(f(0).le(chi)&&chi.le(1),'Bounded native return defect required');
  assert.ok(chi.eq(sc(mm(N.dagger(N.sub(I,P)),N.sub(I,P))).div(4)),'Native defect norm disagrees');
  assert.ok(information.eq(informationDecode(readings)),'Four-scalar information decoder disagrees with native matrices');
  assert.ok(shape.eq(scalar.delta.div(f(1).sub(scalar.a.pow(2)).mul(f(1).sub(scalar.b.pow(2))))),'Scalar shape decoder disagrees');
  return {Hu,Hv,Px,Py,P,Dx,Dy,gram:{xx:gxx,xy:gxy,yy:gyy},area_squared:area,return_defect:chi,information_density:information,shape,readings};
}
function observe(U,V){
  const o=observeInformation(U,V),curvature=readout(o.P),invariant=pair(curvature,curvature).div(o.area_squared.mul(4));
  const readings=o.readings;
  assert.ok(invariant.eq(scalarDecode(readings)),'Four-scalar decoder disagrees with native matrices');
  assert.ok(o.information_density.eq(f(1).sub(o.return_defect).mul(invariant)),'Return-defect complement identity failed');
  return {...o,curvature,invariant};
}
const edges={x:[0,1],r:[1,2],t:[3,2],y:[0,3]};
const changeFrame=(E,G)=>Object.fromEntries(Object.entries(E).map(([name,Q])=>{const [s,t]=edges[name];return [name,mm(N.dagger(G[t]),Q,G[s])];}));
const frameSets=[
  [I,N.cayley(R,'1/5'),N.cayley(J,'-2/7'),N.cayley(L,'3/8')],
  [mm(N.cayley(R,'1/3'),N.cayley(J,'1/4')),N.cayley(J,'2/5'),N.cayley(L,'-1/6'),mm(N.cayley(L,'2/7'),N.cayley(R,'-1/2'))]
];
const constantLinks=(q,p)=>{const x=N.cayley(R,f(q).neg()),y=N.cayley(J,f(p).neg());return {x,y,t:x,r:y};};
const frozen=[];let frameCases=0;
for(const q of ['1/4','1/3','1/2','2/3','1','2'])for(const p of ['1/4','1/3','1/2','2/3','1','2']){
  const U=constantLinks(q,p),V=identityLinks(),tag=q+','+p;
  if(q==='1'&&p==='1'){control('R6:antipodal_loop_chart_is_rejected',N.equal(loop(U),N.scale(I,-1))&&rejects(()=>observe(U,V)));continue;}
  const o=observe(U,V),qf=f(q),pf=f(p),den=f(1).add(qf.pow(2)).pow(2).mul(f(1).add(pf.pow(2)).pow(2)).sub(qf.pow(2).mul(pf.pow(2)).mul(16));
  check('R6:exact_finite_formula:'+tag,o.invariant.eq(f(1).div(den)));
  check('R2:four_scalar_reconstruction:'+tag,scalarDecode(o.readings).eq(o.invariant));
  for(let j=0;j<frameSets.length;j++){
    const G=frameSets[j],moved=observe(changeFrame(U,G),changeFrame(V,G)),conj=A=>mm(N.dagger(G[0]),A,G[0]);
    check('R1:independent_vertex_frames:'+tag+':'+j,moved.invariant.eq(o.invariant)&&N.equal(moved.P,conj(o.P))&&N.equal(moved.Dx,conj(o.Dx))&&N.equal(moved.Dy,conj(o.Dy))&&N.equal(moved.curvature,conj(o.curvature)));
    frameCases++;
  }
  frozen.push({q:qf,p:pf,readings:o.readings,invariant:o.invariant,area_squared:o.area_squared});
}
const benchmark=observe(constantLinks('1/2','1/3'),identityLinks());
check('R6:benchmark_324_over_481',benchmark.invariant.eq('324/481'));
for(const [h,k] of [['1/7','3'],['2/3','4/5'],['5','11']]){
  const dx=N.scale(benchmark.Dx,f(1).div(h)),dy=N.scale(benchmark.Dy,f(1).div(k)),Fxy=N.scale(benchmark.curvature,f(1).div(f(h).mul(k)));
  const area=pair(dx,dx).mul(pair(dy,dy)).sub(pair(dx,dy).pow(2));
  check('R2:external_spacing_cancels:'+h+','+k,pair(Fxy,Fxy).div(area.mul(4)).eq(benchmark.invariant));
}
// Curved reference is permitted; only the relative value is reconstructed.
const pairedSource=constantLinks('1/2','1/3'),curvedReference=constantLinks('1/5','1/7'),paired=observe(pairedSource,curvedReference);
check('R1:nonflat_reference_admitted',!N.equal(paired.Hv,I));
frameSets.forEach((G,j)=>{check('R1:curved_reference_frame_covariance:'+j,observe(changeFrame(pairedSource,G),changeFrame(curvedReference,G)).invariant.eq(paired.invariant));frameCases++;});
control('R1:changing_the_physical_reference_changes_the_comparison',!paired.invariant.eq(benchmark.invariant));
const oneFrame=[I,N.cayley(R,'1/6'),I,I];
control('R1:omitting_the_reference_frame_change_breaks_the_quotient',!observe(changeFrame(pairedSource,oneFrame),identityLinks()).invariant.eq(benchmark.invariant));

// Realize the same a,b,d but two different cross readings c on full cells.
function prescribedLoop(Px,Py,P){return {x:Px,y:Py,r:I,t:mm(Px,N.dagger(P),N.dagger(Py))};}
const px=N.cayley(R,'-1/2'),py=N.cayley(J,'-1/3'),tiltedPy=N.cayley(combo(['3/5','4/5',0]),'-1/3'),targetLoop=N.cayley(L,'-1/4');
const crossA=observe(prescribedLoop(px,py,targetLoop),identityLinks()),crossB=observe(prescribedLoop(px,tiltedPy,targetLoop),identityLinks());
check('R2:prescribed_full_loops_realized',N.equal(crossA.P,targetLoop)&&N.equal(crossB.P,targetLoop));
control('R2:three_scalars_without_cross_reading_are_insufficient',['a','b','d'].every(k=>crossA.readings[k].eq(crossB.readings[k]))&&!crossA.readings.c.eq(crossB.readings.c)&&crossB.invariant.eq(crossA.invariant.mul('25/16')));

// R7: exact shrinking-cell families; the all-refinements bounds are written proofs.
const concentrationRows=[],collapseRows=[];
check('R7:benchmark_information_324_over_625',benchmark.information_density.eq('324/625'));
for(let n=1;n<=6;n++){
  const t=f(1).div(f(2).pow(n)),t2=t.pow(2),t4=t.pow(4);
  const edgeX=N.cayley(R,t.div(-2)),edgeY=N.cayley(J,t.div(-2));
  const regularLinks=prescribedLoop(edgeX,edgeY,N.cayley(L,t2.neg()));
  const concentratedLinks=prescribedLoop(edgeX,edgeY,N.cayley(L,t.neg()));
  const regular=observe(regularLinks,identityLinks()),concentrated=observe(concentratedLinks,identityLinks());
  check('R7:regular_exact_density:'+n,regular.information_density.eq(f(1).div(f(1).add(t4)))&&regular.area_squared.eq(t4)&&regular.shape.eq(1));
  check('R7:concentration_exact_density:'+n,concentrated.information_density.eq(f(1).div(t2.mul(f(1).add(t2))))&&concentrated.return_defect.eq(t2.div(f(1).add(t2))));
  check('R7:nondegenerate_normalized_comparisons:'+n,N.equal(concentrated.Dx,N.scale(R,t))&&N.equal(concentrated.Dy,N.scale(J,t))&&concentrated.shape.eq(1)&&concentrated.area_squared.eq(t4));
  check('R7:concentration_lower_bound:'+n,f('4/5').div(t2).le(concentrated.information_density));
  check('R7:concentration_curvature_growth:'+n,N.equal(N.scale(concentrated.curvature,f(1).div(t2)),N.scale(L,f(2).div(t))));
  for(let j=0;j<frameSets.length;j++){
    const G=frameSets[j],moved=observeInformation(changeFrame(concentratedLinks,G),changeFrame(identityLinks(),G));
    check('R7:concentration_vertex_frames:'+n+':'+j,moved.information_density.eq(concentrated.information_density)&&moved.return_defect.eq(concentrated.return_defect)&&moved.shape.eq(1));frameCases++;
  }
  const directionR=f(1).sub(t2).div(f(1).add(t2)),directionJ=t.mul(2).div(f(1).add(t2));
  const collapsed=observe(prescribedLoop(edgeX,N.cayley(combo([directionR,directionJ,0]),t.div(-2)),N.cayley(L,t2.neg())),identityLinks());
  check('R7:collapsing_direction_exact_density:'+n,collapsed.shape.eq(directionJ.pow(2))&&collapsed.information_density.eq(f(1).div(f(1).add(t4).mul(directionJ.pow(2))))&&N.equal(N.scale(collapsed.curvature,f(1).div(t2)),N.scale(L,2)));
  concentrationRows.push({n,step:t,regular,concentrated});collapseRows.push({n,step:t,observation:collapsed});
}
const firstConcentration=concentrationRows[0].concentrated,lastConcentration=concentrationRows.at(-1).concentrated;
control('R7:vanishing_loop_defect_does_not_bound_its_density',lastConcentration.return_defect.le(firstConcentration.return_defect)&&!lastConcentration.return_defect.eq(firstConcentration.return_defect)&&firstConcentration.information_density.le(lastConcentration.information_density)&&!firstConcentration.information_density.eq(lastConcentration.information_density));
control('R7:density_alone_does_not_exclude_direction_collapse',collapseRows.at(-1).observation.shape.le(collapseRows[0].observation.shape)&&collapseRows[0].observation.information_density.le(collapseRows.at(-1).observation.information_density));
const antipodalLinks=prescribedLoop(N.cayley(R,'-1/2'),N.cayley(J,'-1/2'),N.scale(I,-1));
const antipodal=observeInformation(antipodalLinks,identityLinks());
check('R7:loop_antipode_has_finite_information_density',antipodal.return_defect.eq(1)&&antipodal.area_squared.eq(1)&&antipodal.information_density.eq(1));
control('R7:Cayley_pole_is_not_information_concentration',rejects(()=>observe(antipodalLinks,identityLinks()))&&antipodal.information_density.eq(1));
frameSets.forEach((G,j)=>{check('R7:antipodal_density_vertex_frames:'+j,observeInformation(changeFrame(antipodalLinks,G),changeFrame(identityLinks(),G)).information_density.eq(1));frameCases++;});
control('R7:information_decoder_rejects_nonunit_loop_scalar',rejects(()=>informationDecode({...benchmark.readings,d:'-3/2'})));

// True variable-field plaquettes, evaluated with native rational Cayley links.
function fieldLinks(field,h,k){h=f(h);k=f(k);return {
  x:N.cayley(field.x(f(0),f(0)),h.div(-2)),y:N.cayley(field.y(f(0),f(0)),k.div(-2)),
  r:N.cayley(field.y(h,f(0)),k.div(-2)),t:N.cayley(field.x(f(0),k),h.div(-2))};}
const fields=[
  {name:'nonclosed_two_coefficient_forms',x:(x,y)=>N.scale(R,f(1).add(y)),y:(x,y)=>N.scale(J,f(1).add(x)),curvature:combo([-1,1,2]),limit:f('3/2')},
  {name:'closed_scalar_phase_forms',x:(x,y)=>N.scale(R,f(1).add(y)),y:(x,y)=>N.add(N.scale(R,x),J),curvature:N.scale(L,2),limit:f(1)}
];
const refinements=[];
for(const field of fields){
  const rows=[];let lastError=null;
  for(const h of ['1/4','1/8','1/16','1/32','1/64']){
    const o=observe(fieldLinks(field,h,h),identityLinks()),Fapprox=N.scale(o.curvature,f(1).div(f(h).pow(2))),error=N.energy(N.sub(Fapprox,field.curvature));
    check('R3:variable_field_orientation:'+field.name+':'+h,error.le(N.energy(N.add(Fapprox,field.curvature))));
    if(lastError!==null)check('R5:specified_refinement_error_decreases:'+field.name+':'+h,error.le(lastError)&&!error.eq(lastError));
    lastError=error;rows.push({step:f(h),invariant:o.invariant,curvature_estimate:Fapprox,squared_native_error:error});
  }
  check('R5:specified_refinement_near_target:'+field.name,rows.at(-1).invariant.sub(field.limit).abs().le('1/10'));
  refinements.push({name:field.name,expected_curvature:field.curvature,expected_limit:field.limit,rows});
}
// Flat endpoint transport with independent base comparisons defeats a frozen-only diagnosis.
const vertices=[I,N.cayley(R,'1/4'),mm(N.cayley(R,'1/4'),N.cayley(J,'1/5')),N.cayley(J,'1/5')];
const pureGauge=Object.fromEntries(Object.entries(edges).map(([e,[s,t]])=>[e,mm(N.dagger(vertices[t]),vertices[s])]));
const flat=observe(pureGauge,identityLinks());
check('R3:endpoint_generated_loop_is_exactly_flat',N.equal(flat.P,I)&&flat.invariant.zero()&&!flat.area_squared.zero());
control('R4:frozen_base_bracket_cannot_detect_derivative_cancellation',!N.isZero(N.commutator(flat.Dx,flat.Dy))&&flat.invariant.zero());
const rawBracket=N.commutator(R,J),cancelled=N.sub(rawBracket,N.scale(L,2));
check('R4:third_channel_curl_reconstruction',pair(L,rawBracket).sub(2).zero()&&pair(L,cancelled).sub(2).eq(-2));
check('R4:comparator_must_transform_with_curvature',frameSets.every(G=>pair(mm(N.dagger(G[0]),L,G[0]),mm(N.dagger(G[0]),rawBracket,G[0])).eq(pair(L,rawBracket))));

control('domain:zero_area_is_rejected',rejects(()=>observe({x:px,y:px,t:px,r:px},identityLinks())));
control('domain:edge_antipode_is_rejected',rejects(()=>observe({x:N.scale(I,-1),y:py,t:px,r:py},identityLinks())));
control('domain:nonunit_transport_is_rejected',rejects(()=>observe({...pairedSource,x:N.scale(pairedSource.x,2)},identityLinks())));
control('domain:outside_real_unit_word_sector_is_rejected',rejects(()=>observe({...pairedSource,x:N.scale(I,IOTA)},identityLinks())));
control('domain:inconsistent_scalar_area_is_rejected',rejects(()=>scalarDecode({a:'3/5',b:'4/5',c:1,d:'1/2'})));
control('domain:loop_antipode_scalar_is_rejected',rejects(()=>scalarDecode({...benchmark.readings,d:-1})));
control('scope:finite_quotient_is_not_its_continuum_calibration',!benchmark.invariant.eq(1));
const resampleStep=f('1/4');
control('R5:pointwise_resampling_is_not_an_exact_finite_frame_change',!N.equal(N.cayley(R,resampleStep.add(resampleStep.pow(2)).neg()),N.cayley(N.scale(R,2),resampleStep.div(-2))));
control('scope:scalar_magnitude_cannot_supply_oriented_curvature',pair(L,rawBracket).eq(2)&&pair(L,N.scale(rawBracket,-1)).eq(-2)&&pair(rawBracket,rawBracket).eq(pair(N.scale(rawBracket,-1),N.scale(rawBracket,-1))));
const c0=certificates.find(x=>x.certificate.steps.length>0),system=systems.get(c0.certificate.presentation);
const bad=JSON.parse(JSON.stringify(c0.certificate));bad.steps[0].coefficient=['997','0'];
control('certificate:changed_rewrite_step_is_rejected',rejects(()=>system.replay(bad)));
const wrong=JSON.parse(JSON.stringify(c0.certificate));wrong.presentation='0'.repeat(64);
control('certificate:wrong_presentation_is_rejected',rejects(()=>system.replay(wrong)));
control('native:floating_point_scalar_is_rejected',rejects(()=>f(0.25)));

const output={schema:'native-relative-curvature.v2',status:'PASS_NATIVE_RELATIVE_CURVATURE',source_commit:pins.source_commit,
  scope:'Exact finite native source/reference curvature and return-defect information observables, algebraic variable-field jets and finite rational families. Smooth C2 limits, all-refinements concentration bounds and continuation obstructions have written proofs in THEOREM.md. No gravitational singularity, dynamics, physical metric identification or quantum-state selection is certified.',
  source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),native_engine:pins.native_engine,upstream,
  presentations,certificates,checks,negative_controls:controls,
  examples:{benchmark,constant_cells:frozen,curved_reference:paired,cross_reading_obstruction:{first:crossA,second:crossB},variable_field_refinements:refinements,flat_endpoint_transport:flat,information_concentration:concentrationRows,direction_collapse:collapseRows,loop_antipode:antipodal},
  summary:{symbolic_certificates:certificates.length,serialized_replays:certificates.length,native_presentations:presentations.length,exact_checks:checks.length,negative_controls:controls.length,constant_cells:frozen.length,vertex_frame_cases:frameCases,variable_field_refinement_cases:refinements.reduce((s,x)=>s+x.rows.length,0),concentration_refinements:concentrationRows.length,direction_collapse_refinements:collapseRows.length}};
const serialized=JSON.stringify(output,null,2)+'\n',digest=hash(serialized),file=path.join(root,'CERTIFICATE.json'),checksum=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(file,serialized);fs.writeFileSync(checksum,digest+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(file,'utf8'),serialized,'Frozen certificate differs from native replay');assert.equal(fs.readFileSync(checksum,'utf8'),digest+'  CERTIFICATE.json\n','Certificate checksum mismatch');}
console.log(JSON.stringify({status:output.status,mode:write?'explicit generation':'read-only replay',...output.summary,benchmark:benchmark.invariant,certificate_sha256:digest},null,2));
