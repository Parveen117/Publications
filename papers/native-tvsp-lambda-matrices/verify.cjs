'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const assert=require('node:assert/strict');
const L=require('./tvsp_adapter.cjs');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const P=require('../native-common-cone-limit/native_engine/paninian_operator.cjs');
const argv=process.argv.slice(2);assert.ok(argv.length===1&&['--write','--check'].includes(argv[0]));
const checks=[],controls=[],proofs=[],presentations=[],fixtures=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);controls.push({name,rejected:true});};
const reject=f=>{try{f();return false;}catch{return true;}};
const f=N.F.of,cut=N.Cut.of,mm=(...a)=>a.reduce(N.mul),eq=N.equal;
const tr=a=>a.reduce((v,r,i)=>v.add(r[i]),N.ZERO);
const sc=(ctx,v)=>L.scalar(ctx,v),proj=L.project,mproj=L.matrixProject;
const mx=(a,b)=>L.matrixMul(a,b),ms=(a,b)=>L.matrixSub(a,b);
const comm=(a,b)=>ms(mx(a,b),mx(b,a));
const coeffMatrix=(a,s,v)=>a.rows.map(row=>row.map(x=>L.coefficient(x,s,v)));
const matrixPoint=L.matrixOuter;
function coefficientSame(a,b){
  const xs=L.entries(a),ys=L.entries(b);
  return xs.length===ys.length&&xs.every(([s,v,c])=>c.eq(L.coefficient(b,s,v)));
}

// Holonomic quotient calculus. Matrix noncommutativity remains ordered.
for(const order of [1,2,3,4]){
  const C=L.createContext(order,'fixture common chart'),s=L.coordinate(C,'s'),v=L.coordinate(C,'v');
  const a=L.add(sc(C,2),L.add(s,L.scale(v,N.IOTA)));
  const b=L.add(sc(C,'3/2'),L.sub(L.scale(s,-2),v));
  const c=L.add(sc(C,'-1/3'),L.mul(s,v));
  ck('jet:associativity:'+order,L.equal(L.mul(L.mul(a,b),c),L.mul(a,L.mul(b,c))));
  ck('jet:inverse:'+order,L.equal(L.mul(a,L.inverse(a)),sc(C,1)));
  ck('jet:dagger:'+order,L.equal(L.dagger(L.mul(a,b)),L.mul(L.dagger(b),L.dagger(a))));
  ck('jet:outer_product:'+order,L.outer(L.mul(a,b)).eq(L.outer(a).mul(L.outer(b))));
  ck('jet:projection_product:'+order,L.equal(proj(L.mul(a,b),order-1),L.mul(proj(a,order-1),proj(b,order-1))));
  const lhs=L.derivative(L.mul(a,b),'s');
  const rhs=L.add(L.mul(L.derivative(a,'s'),proj(b,order-1)),L.mul(proj(a,order-1),L.derivative(b,'s')));
  ck('jet:typed_Leibniz:'+order,L.equal(lhs,rhs));
  const A=L.matrix([[a,b],[sc(C,1),c]]),B=L.matrix([[b,sc(C,0)],[c,a]]),D=L.matrix([[c,a],[b,sc(C,1)]]);
  ck('matrix:associativity:'+order,L.matrixEqual(mx(mx(A,B),D),mx(A,mx(B,D))));
  ck('matrix:dagger_reverses:'+order,L.matrixEqual(L.matrixDagger(mx(A,B)),mx(L.matrixDagger(B),L.matrixDagger(A))));
  ck('matrix:ordinary_recovery:'+order,eq(matrixPoint(mx(A,B)),mm(matrixPoint(A),matrixPoint(B))));
  no('matrix:order_can_be_erased:'+order,!L.matrixEqual(mx(A,B),mx(B,A)));
  const RA=L.regularMatrix(A),RB=L.regularMatrix(B),basis=L.basisIndices(C),d=basis.length;
  ck('regular:dimension:'+order,d===(order+1)*(order+2)/2&&RA.length===2*d);
  ck('regular:product:'+order,eq(L.regularMatrix(mx(A,B)),mm(RA,RB)));
  ck('regular:diagonal_repeats_point:'+order,basis.every((_,k)=>eq(RA.slice(2*k,2*k+2).map(r=>r.slice(2*k,2*k+2)),matrixPoint(A))));
  ck('regular:complete_decoder:'+order,basis.every(([s,v],k)=>eq(RA.slice(2*k,2*k+2).map(r=>r.slice(0,2)),coeffMatrix(A,s,v))));
  no('regular:formal_star_is_positive_adjoint:'+order,!eq(L.regularMatrix(L.matrixDagger(A)),N.dagger(RA)));
  ck('jet:coordinate_nilpotent:'+order,L.outer(s).zero()&&L.equal(L.power(s,order+1),sc(C,0))&&!L.equal(s,sc(C,0)));
  no('jet:zero_depth_determines_derivative:'+order,reject(()=>L.derivative(proj(a,0),'s')));
  no('jet:unknown_coefficient_is_zero:'+order,reject(()=>L.coefficient(a,order+1,0)));
  no('jet:implicit_context_merge:'+order,reject(()=>L.add(a,sc(L.createContext(order,'fixture common chart'),1))));
  no('jet:implicit_depth_merge:'+order,reject(()=>L.mul(a,proj(a,order-1))));
  no('jet:invent_higher_layers:'+order,reject(()=>proj(a,order+1)));
  no('jet:singular_inverse:'+order,reject(()=>L.inverse(s)));
}

// Analytic source input is differentiated before the finite response records.
function geometry(z){
  const depth=z.context.order;assert.ok(depth>=2);
  const inv=L.matrix([[L.div(z.c,z.Delta),L.scale(L.div(z.b,z.Delta),-1)],
    [L.scale(L.div(z.b,z.Delta),-1),L.div(z.a,z.Delta)]]);
  const xs=mx(mproj(inv,depth-1),L.matrixDerivative(z.H,'s'));
  const xv=mx(mproj(inv,depth-1),L.matrixDerivative(z.H,'v'));
  const As=L.matrixScale(xs,'1/2'),Av=L.matrixScale(xv,'1/2');
  const direct=L.matrixAdd(ms(L.matrixDerivative(Av,'s'),L.matrixDerivative(As,'v')),mproj(comm(As,Av),depth-2));
  const F=mproj(L.matrixScale(comm(xs,xv),'-1/4'),depth-2);
  const H=mproj(z.H,depth-2),kappa=L.div(mx(H,F).rows[0][1],proj(z.Delta,depth-2));
  return {inv,xs,xv,As,Av,F,direct,kappa};
}
const parameterValues=['-5','0','1/2','1','5'];
let oldVisible=null,oldPoint=null,oldResponse=null;
for(const mu of parameterValues){
  const z=L.buildFixture(mu,4),C=z.context,s=L.coordinate(C,'s'),v=L.coordinate(C,'v'),m=cut(mu),g=geometry(z);
  ck('source:H_ss:'+mu,L.equal(z.a,L.add(sc(C,2),L.scale(v,m))));
  ck('source:H_sv:'+mu,L.equal(z.b,L.add(sc(C,1),L.scale(s,m))));
  ck('source:H_vv:'+mu,L.equal(z.c,sc(C,5)));
  ck('source:Maxwell:'+mu,L.equal(L.derivative(z.T,'v'),L.scale(L.derivative(z.P,'s'),-1)));
  ck('source:inverse:'+mu,L.matrixEqual(mx(z.H,g.inv),L.identity(C,2)));
  ck('source:response_identity:'+mu,L.equal(L.mul(z.responses.lambdaP,z.c),L.scale(L.mul(z.S,z.Delta),-1)));
  ck('source:response_reciprocal:'+mu,L.equal(L.mul(z.responses.lambdaT,z.Delta),L.mul(z.P,z.a)));
  ck('source:visible_response_entire_germ:'+mu,L.equal(z.r,L.scale(z.V,5)));
  ck('source:point_reading_five:'+mu,L.outer(z.r).eq(5));
  ck('source:deformation_from_response:'+mu,L.coefficient(z.responses.lambdaV,0,1).eq(m.neg()));
  ck('source:all_six_point_responses:'+mu,Object.values(z.responses).every((x,i)=>L.outer(x).eq(['-9/5','-2','-5','-9/2','2','20/9'][i])));
  if(oldVisible){
    ck('source:equal_scalar_layers_across_deformations:'+mu,coefficientSame(z.r,oldVisible));
    ck('source:equal_entire_point_matrix:'+mu,eq(matrixPoint(z.H),oldPoint));
    ck('source:equal_six_point_readings:'+mu,Object.values(z.responses).every((x,i)=>L.outer(x).eq(oldResponse[i])));
  }
  oldVisible=z.r;oldPoint=matrixPoint(z.H);oldResponse=Object.values(z.responses).map(L.outer);
  ck('curvature:actual_connection_derivatives:'+mu,L.matrixEqual(g.direct,g.F));
  const predicted=L.scale(L.inverse(L.power(proj(z.Delta,2),2)),m.mul(m).mul('5/4'));
  ck('curvature:whole_chart_jet_formula:'+mu,L.equal(g.kappa,predicted));
  ck('curvature:point_formula:'+mu,L.outer(g.kappa).eq(m.mul(m).mul('5/324')));
  const F0=N.scale(N.matrix([[5,25],[-10,-5]]),m.mul(m).div(324));
  ck('curvature:oriented_matrix_formula:'+mu,eq(matrixPoint(g.F),F0));
  ck('curvature:scalar_trace_is_blind:'+mu,tr(matrixPoint(g.F)).zero());
  if(!m.zero()){
    no('curvature:positive_commutator_sign:'+mu,!eq(matrixPoint(g.F),N.scale(F0,-1)));
    no('curvature:connection_commutator_alone:'+mu,!eq(matrixPoint(g.F),N.commutator(matrixPoint(g.As),matrixPoint(g.Av))));
    no('curvature:scalar_equals_deformation:'+mu,!L.outer(g.kappa).eq(m));
  }
  // Check actual NT covariant second layer on the visible mechanical response.
  const ds=L.derivative(z.r,'s'),dv=L.derivative(z.r,'v');
  const Gs=mproj(g.As,2),Gv=mproj(g.Av,2);
  const hess=(direction,j)=>{
    const G=direction==='s'?Gs:Gv;
    const second=L.derivative(j===0?ds:dv,direction);
    return L.sub(second,L.add(L.mul(G.rows[0][j],proj(ds,2)),L.mul(G.rows[1][j],proj(dv,2))));
  };
  ck('tower:torsion_free_second_slots:'+mu,L.equal(hess('s',1),hess('v',0)));
  if(!m.zero())no('tower:partial_hessian_is_full_covariant_layer:'+mu,!L.outer(hess('s',0)).zero());
  const tower=L.responseTower(z,3);
  ck('tower:typed_storage_depths:'+mu,tower.every((layer,k)=>layer.derivativeOrder===k&&layer.depth===4-k));
  ck('tower:actual_second_layer:'+mu,L.equal(tower[2].components.lambdaS['00'],L.scale(hess('s',0),-1)));
  for(const [name,response] of Object.entries(z.responses)){
    for(let j=0;j<2;j++){
      const defect=L.sub(tower[3].components[name]['01'+j],tower[3].components[name]['10'+j]);
      const Fc=mproj(g.F,1),d=[proj(L.derivative(response,'s'),1),proj(L.derivative(response,'v'),1)];
      const expected=L.scale(L.add(L.mul(Fc.rows[0][j],d[0]),L.mul(Fc.rows[1][j],d[1])),-1);
      ck('tower:ordered_slot_curvature:'+mu+':'+name+':'+j,L.equal(defect,expected));
    }
  }
  const H1=mproj(z.H,1),reg=L.regularMatrix(H1),H0=matrixPoint(H1),hs=coeffMatrix(H1,1,0),hv=coeffMatrix(H1,0,1);
  ck('one_matrix:exact_first_blocks:'+mu,eq(reg.slice(2,4).map(r=>r.slice(0,2)),hs)&&eq(reg.slice(4,6).map(r=>r.slice(0,2)),hv));
  for(let k=1;k<=6;k++)ck('spectrum:power_trace_blind:'+mu+':'+k,tr(N.power(reg,k)).eq(tr(N.power(H0,k)).mul(3)));
  const characteristic=N.add(N.sub(N.power(reg,2),N.scale(reg,7)),N.scale(N.identity(6),9));
  ck('spectrum:characteristic_polynomial:'+mu,N.isZero(N.power(characteristic,3)));
  fixtures.push({mu,point_H:H0,H_s:hs,H_v:hv,point_responses:oldResponse,reading:L.outer(z.r),curvature:L.outer(g.kappa),F:F0,first_layer_matrix:reg,
    response_lambdaV:L.serialize(z.responses.lambdaV),visible_response:L.serialize(z.r)});
}
const five=fixtures.find(x=>x.mu==='5'),negative=fixtures.find(x=>x.mu==='-5'),flat=fixtures.find(x=>x.mu==='0');
ck('example:mu_five_curvature',five.curvature.eq('125/324'));
no('target:point_or_entire_scalar_germ_determines_geometry',five.reading.eq(flat.reading)&&!five.curvature.eq(flat.curvature));
no('target:curvature_retains_deformation_sign',five.curvature.eq(negative.curvature)&&!eq(five.H_s,negative.H_s));
no('target:flatness_means_all_lambda_responses_zero',flat.curvature.zero()&&flat.point_responses.some(x=>!x.zero()));
no('target:integer_layer_order_quantizes_deformation',fixtures.find(x=>x.mu==='1/2').curvature.eq('5/1296'));
{
  const z=L.buildFixture(0,0);
  ck('tower:zero_order_readout',L.responseTower(z,0).length===1);
  no('tower:unretained_covariant_derivatives',reject(()=>L.responseTower(z,1)));
}
ck('patch:positive_uniform_a',f('7/4').le(2));
ck('patch:positive_uniform_Delta',f(5).mul('7/4').sub(f('5/4').pow(2)).eq('115/16')&&f(0).le('115/16'));

// Fixed channel frames, and a concrete failure of independent entry frames.
{
  const C=L.createContext(2),s=L.coordinate(C,'s'),v=L.coordinate(C,'v');
  const A=L.matrix([[L.add(sc(C,2),s),sc(C,1)],[v,sc(C,5)]]),B=L.matrix([[v,s],[sc(C,1),sc(C,2)]]);
  const G=L.matrix([[sc(C,'3/5'),sc(C,'-4/5')],[sc(C,'4/5'),sc(C,'3/5')]]),gd=L.matrixDagger(G);
  ck('frame:unitary',L.matrixEqual(mx(gd,G),L.identity(C,2)));
  ck('frame:matching_cancel',L.matrixEqual(mx(mx(B,G),mx(gd,A)),mx(B,A)));
  const a=N.matrix([[0,1],[0,0]]),b=N.matrix([[0,0],[1,0]]);
  const phi=x=>[[x[0][0].neg(),x[0][1]],[x[1][0],x[1][1]]];
  no('frame:arbitrary_entrywise_action_composes',!eq(phi(mm(a,b)),mm(phi(a),phi(b))));
}

// Independent source connection selection and differentiation depth controls.
{
  const C=L.createContext(2),s=L.coordinate(C,'s'),one=sc(C,1),zero=sc(C,0);
  const Av=L.matrix([[zero,L.scale(s,-1)],[s,zero]]);
  const F=matrixPoint(L.matrixDerivative(Av,'s'));
  ck('connection:constant_H_compatible_other_transport',eq(N.add(F,N.dagger(F)),N.zeros(2)));
  no('connection:entire_H_jet_fixes_arbitrary_connection',!N.isZero(F));
  const zeroRelation=L.power(s,3);
  ck('depth:truncation_relation',L.equal(zeroRelation,zero));
  no('depth:derivative_is_same_order_derivation',!L.equal(L.scale(L.power(s,2),3),zero));
  const R=L.regularMatrix(L.matrix([[s]]));
  ck('positivity:selfadjoint_formal_coordinate',L.equal(L.dagger(s),s));
  ck('positivity:nonzero_nilpotent_regular_action',!N.isZero(R)&&N.isZero(N.power(R,3)));
  no('positivity:regular_action_is_selfadjoint',!eq(R,N.dagger(R)));
  no('state:operator_products_stay_symmetric_Hessians',!eq(mm(N.matrix([[2,1],[1,5]]),N.matrix([[1,0],[0,2]])),mm(N.matrix([[1,0],[0,2]]),N.matrix([[2,1],[1,5]]))));
}

// Serialized symbolic identities replayed by the unchanged Paninian engine.
function register(S,name){const audit=S.audit();ck('words:confluent:'+name,audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:S.spec(),audit});}
function prove(S,name,a,b){const r=S.reduce(a.minus(b),{witness:true});ck('words:'+name,r.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(r.certificate));ck('words:replay:'+name,S.replay(certificate));proofs.push({name,presentation_hash:S.hash(),lhs:a.toJSON(),rhs:b.toJSON(),certificate});}
{
  const S=new P.Presentation({tokens:['x','y','z'],rules:[]});register(S,'ordered_connection');
  const x=S.word(['x']),y=S.word(['y']),z=S.word(['z']),br=P.bracket;
  const direct=z.minus(x.times(y)).minus(z.minus(y.times(x))).scale('1/2').plus(br(x,y).scale('1/4'));
  prove(S,'selected_curvature_sign_and_factor',direct,br(x,y).scale('-1/4'));
  no('words:wrong_sign',S.reduce(direct.minus(br(x,y).scale('1/4'))).normal.terms.size>0);
}
{
  const tokens=['a','b','m'],rules=[];
  for(let j=1;j<3;j++)for(let i=0;i<j;i++)rules.push({id:'commute_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'central source coefficients'});
  const S=new P.Presentation({tokens,rules});register(S,'whole_source_family');
  const a=S.word(['a']),b=S.word(['b']),m=S.word(['m']),z=S.zero(),one=S.one();
  const sm=(A,B)=>A.map(row=>B[0].map((_,j)=>row.reduce((out,x,k)=>out.plus(x.times(B[k][j])),z)));
  const H=[[a,b],[b,one.scale(5)]],B=[[one.scale(5),b.scale(-1)],[b.scale(-1),a]],ds=[[z,m],[m,z]],dv=[[m,z],[z,z]];
  const X=sm(B,ds),Y=sm(B,dv),XY=sm(X,Y),YX=sm(Y,X);
  const G=[[b.scale(5),one.scale(25)],[a.scale(-5),b.scale(-5)]],HG=sm(H,G),delta=a.scale(5).minus(b.times(b));
  for(let i=0;i<2;i++)for(let j=0;j<2;j++)prove(S,'source_commutator_'+i+j,XY[i][j].minus(YX[i][j]),m.times(m).times(G[i][j]).scale(-1));
  prove(S,'source_curvature_marker',HG[0][1],delta.scale(5));
  prove(S,'source_metric_skew',HG[1][0],HG[0][1].scale(-1));
}
{
  const C=L.createContext(2),S=new P.Presentation(L.presentationSpec(C));register(S,'actual_depth_two_quotient');
  const s=S.word(['s']),v=S.word(['v']);
  prove(S,'formal_coordinates_commute',s.times(v),v.times(s));
  prove(S,'finite_layer_relation',s.times(v).times(s),S.zero());
}
const payload={protocol:'TLM1_TVSP_LAMBDA_MATRICES_V1',status:'PASS_TLM1_TVSP_LAMBDA_MATRICES',
  source_parent:'35d63e4ac0d8cf161c549654ce6b430d68afbde4',checks,negative_controls:controls,presentations,proofs,fixtures,
  scope:{source:'NT six actual responses and selected Hessian strain connection',layers:'source-generated finite jets with typed covariant reconstruction',
    ordinary_recovery:true,regular_matrix_faithful:true,regular_matrix_positive_star_representation:false,
    point_spectrum_determines_geometry:false,scalar_entire_germ_determines_geometry:false,
    physical_space_quantized:false,quantum_gravity_closed:false,universal_claims_from_finite_tests:false},
  summary:{exact_checks:checks.length,symbolic_certificates:proofs.length,negative_controls:controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex');
const cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(argv[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}
else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
