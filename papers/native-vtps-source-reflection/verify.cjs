'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs'),J=require('../native-tvsp-lambda-matrices/tvsp_adapter.cjs'),A=require('./reflection_adapter.cjs');
const {F}=N,args=process.argv.slice(2);assert.ok(args.length===1&&['--write','--check'].includes(args[0]));
const checks=[],negative_controls=[],fixtures=[],ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});},no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const rejects=f=>{try{f();return false;}catch{return true;}},point=J.matrixOuter;
const coeff=(g,a,b)=>g.rows.map(row=>row.map(x=>J.coefficient(x,a,b)));
const sources=[];
for(const raw of [0,1,5,-5]){
 const z=A.source(raw,4),mu=F.of(raw),C=z.context,sc=x=>J.scalar(C,x),g=z.metric,M=z.response;
 const gamma=A.leviCivita(g),Fcurv=A.curvature(gamma),kappa=A.scalarCurvature(z,Fcurv),trans=A.transportedConnection(z);
 ck('source:temperature_inverse:'+raw,J.equal(J.add(J.mul(z.a,z.s),z.v),z.theta));
 ck('source:Helmholtz_P:'+raw,J.equal(J.derivative(z.potential,'s'),J.scale(J.project(z.P,3),-1)));
 ck('source:Helmholtz_S:'+raw,J.equal(J.derivative(z.potential,'v'),J.scale(J.project(z.S,3),-1)));
 const direct=J.matrix([[J.derivative(z.P,'s'),J.derivative(z.P,'v')],[J.derivative(z.S,'s'),J.derivative(z.S,'v')]]);
 ck('Maxwell:actual_response_matrix:'+raw,J.matrixEqual(direct,J.matrixProject(M,3)));
 ck('Maxwell:symmetric_cross_response:'+raw,J.equal(direct.rows[0][1],direct.rows[1][0]));
 ck('metric:pullback:'+raw,J.matrixEqual(g,J.matrixMul(J.matrixMul(J.matrixDagger(z.Q),z.H),z.Q)));
 ck('metric:positive_origin:'+raw,J.outer(g.rows[0][0]).eq('9/2')&&J.outer(g.rows[1][1]).eq('1/2'));
 ck('response:determinant:'+raw,J.equal(A.determinant(M),J.scale(J.inverse(z.a),-5)));
 ck('connection:actual_transport:'+raw,gamma.every((x,i)=>J.matrixEqual(x,trans[i])));
 ck('curvature:tensor_transport:'+raw,J.matrixEqual(Fcurv,J.matrixProject(A.transformedCurvature(z),2)));
 ck('curvature:whole_chart_jet:'+raw,J.equal(kappa,J.project(J.div(sc(mu.pow(2).mul('5/4')),J.power(z.Delta,2)),2)));
 ck('curvature:origin:'+raw,J.outer(kappa).eq(mu.pow(2).mul('5/324')));
 ck('reading:actual_five:'+raw,J.outer(z.reading).eq(5)&&J.equal(J.div(z.reading,z.V),sc(5)));
 ck('metric:first_v_layer:'+raw,N.equal(coeff(g,1,0),N.scale(N.matrix([[3,0],[0,-1]]),mu.div(4))));
 ck('metric:first_theta_layer:'+raw,N.equal(coeff(g,0,1),N.scale(N.matrix([[-1,0],[0,0]]),mu.div(2))));
 const expectedF=N.scale(N.matrix([[0,5],[-45,0]]),mu.pow(2).div(648));
 ck('curvature:explicit_VTPS_matrix:'+raw,N.equal(point(Fcurv),expectedF));
 const naive=['s','v'].map(x=>J.matrixScale(J.matrixMul(J.matrixProject(A.inverse2(g),3),J.matrixDerivative(g,x)),'1/2'));
 ck('control:reselected_diagonal_strain_is_flat:'+raw,J.matrixEqual(A.curvature(naive),J.matrixScale(J.matrixProject(g,2),0)));
 if(!mu.zero()){
  no('inference:diagonal_metric_is_flat:'+raw,!J.outer(kappa).zero());
  no('inference:drop_connection_derivative:'+raw,!gamma.every((x,i)=>J.matrixEqual(x,A.transportedConnection(z,false)[i])));
  no('inference:response_matrix_is_positive_metric:'+raw,J.outer(A.determinant(M)).rad.le(0));
 }
 sources.push({mu,point_metric:point(g),point_response:point(M),metric_v:coeff(g,1,0),metric_theta:coeff(g,0,1),regular_first_metric:J.regularMatrix(J.matrixProject(g,1)),point_connection:gamma.map(point),point_curvature:point(Fcurv),kappa:J.outer(kappa),reading:J.outer(z.reading)});
}
fixtures.push({kind:'actual_VTPS_five_source',sources});
const positive=A.source(5,4),negative=A.source(-5,4);
const centered=z=>J.add(z.potential,J.add(J.add(J.scalar(z.context,10),z.theta),J.scale(z.v,10)));
ck('reflection:centered_potential',A.reflectionEqual(centered(negative),centered(positive)));
ck('reflection:metric_all_retained_layers',A.matrixReflectionEqual(negative.metric,positive.metric));
ck('reflection:response_all_retained_layers',A.matrixReflectionEqual(negative.response,positive.response));
ck('reflection:centered_pressure',A.reflectionEqual(J.sub(negative.P,J.scalar(negative.context,10)),J.sub(positive.P,J.scalar(positive.context,10)),1));
ck('reflection:centered_entropy',A.reflectionEqual(J.sub(negative.S,J.scalar(negative.context,1)),J.sub(positive.S,J.scalar(positive.context,1)),1));
const gp=A.leviCivita(positive.metric),gn=A.leviCivita(negative.metric);
ck('reflection:connection_one_form',gp.every((m,i)=>m.rows.every((row,k)=>row.every((x,j)=>A.reflectionEqual(gn[i].rows[k][j],x,1)))));
const reflection=[];
for(const depth of [1,2]){
 const p=J.matrixProject(positive.metric,depth),n=J.matrixProject(negative.metric,depth),Rp=J.regularMatrix(p),Rn=J.regularMatrix(n),D=A.grade(p.context);
 ck('reflection:involution:'+depth,N.equal(N.mul(D,D),N.identity(D.length)));
 ck('reflection:regular_similarity:'+depth,N.equal(Rn,N.mul(N.mul(D,Rp),D)));
 const x=N.matrix(Array.from({length:D.length},(_,i)=>[i+1])),O=N.matrix([Array.from({length:D.length},(_,i)=>i%2?2:-1)]);
 const before=N.mul(N.mul(O,Rp),x),after=N.mul(N.mul(N.mul(O,D),Rn),N.mul(D,x));
 ck('observer:carried_reading:'+depth,N.equal(before,after));
 no('inference:entrywise_same_two_sides:'+depth,!N.equal(Rp,Rn));
 reflection.push({depth,grade:D,positive:Rp,negative:Rn,before,after});
}
fixtures.push({kind:'graded_matrix_reflection',cases:reflection});
const oldPlus=J.buildFixture(5,1),oldMinus=J.buildFixture(-5,1),oldD=A.grade(oldPlus.context),oldP=J.regularMatrix(oldPlus.H),oldN=J.regularMatrix(oldMinus.H);
ck('inherited_five:original_six_matrix_reflection',N.equal(oldN,N.mul(N.mul(oldD,oldP),oldD)));
fixtures.push({kind:'inherited_six_by_six',grade:oldD,positive:oldP,negative:oldN});
const mirrors=[];
for(const [v,t] of [['1/100','1/50'],['-1/100','1/100'],['1/50','-1/100']]){
 const p=A.exactPoint(5,v,t),n=A.exactPoint(-5,F.of(v).neg(),F.of(t).neg());
 ck('point:metric_isometry:'+v+':'+t,N.equal(p.metric,n.metric)&&N.equal(p.response,n.response));
 ck('point:curvature_match:'+v+':'+t,p.kappa.eq(n.kappa));
 ck('point:centered_responses:'+v+':'+t,p.P.add(n.P).eq(20)&&p.S.add(n.S).eq(2));
 ck('point:calibrated_five:'+v+':'+t,p.reading.div(p.V).eq(5)&&n.reading.div(n.V).eq(5));
 no('inference:uncarried_off_cut_five_is_same:'+v+':'+t,!p.reading.eq(n.reading));
 mirrors.push({positive:p,negative:n});
}
fixtures.push({kind:'exact_off_cut_mirrors',cases:mirrors});
no('guard:singular_temperature_chart',rejects(()=>A.exactPoint(5,'-2/5',0)));
no('guard:nonpositive_Hessian',rejects(()=>A.exactPoint(5,0,2)));
no('guard:insufficient_derivative_record',rejects(()=>A.source(5,1)));
const payload={protocol:'VSR1_VTPS_SOURCE_REFLECTION_V1',status:'PASS_VSR1_VTPS_SOURCE_REFLECTION',source_parent:'0c98adf69682f6fee48be85dcbcb7667d77dea2e',source_pins_sha256:crypto.createHash('sha256').update(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))).digest('hex'),checks,negative_controls,fixtures,summary:{exact_checks:checks.length,symbolic_certificates:0,negative_controls:negative_controls.length,source_deformations:sources.length,off_cut_mirror_pairs:mirrors.length},scope:{actual_constitutive_source:true,VTPS_order:['V','T','P','S'],exact_reflection:'(mu,v,theta) -> (-mu,-v,-theta)',all_finite_grade_similarity:'written theorem; depths 1 and 2 executed',diagonal_metric_has_nonzero_curvature:true,Maxwell_is_response_integrability:true,mu_is_thermodynamic_lambda_component:false,negative_derivative_orders:false,arbitrary_source_cut_symmetry:false,physical_interpretation:'open',local_scientific_replays:0}};
const bytes=JSON.stringify(payload,(_,v)=>typeof v==='bigint'?v.toString():v,2)+'\n',hash=crypto.createHash('sha256').update(bytes).digest('hex'),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(args[0]==='--write'){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
