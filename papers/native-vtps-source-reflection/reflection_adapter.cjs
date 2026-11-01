'use strict';
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const J=require('../native-tvsp-lambda-matrices/tvsp_adapter.cjs');
const {F,Cut}=N;
const sum=xs=>xs.reduce((a,b)=>J.add(a,b));
function inverse2(g){
 const [[a,b],[c,d]]=g.rows,det=J.sub(J.mul(a,d),J.mul(b,c));
 return J.matrix([[J.div(d,det),J.scale(J.div(b,det),-1)],[J.scale(J.div(c,det),-1),J.div(a,det)]]);
}
function determinant(g){const [[a,b],[c,d]]=g.rows;return J.sub(J.mul(a,d),J.mul(b,c));}
function source(mu,order=4){
 mu=F.of(mu);if(!Number.isSafeInteger(order)||order<2)throw new RangeError('At least two retained layers required');
 const context=J.createContext(order,'VTPS source: adapter s=v, adapter v=theta; mu='+mu);
 const c=x=>J.scalar(context,x),v=J.coordinate(context,'s'),theta=J.coordinate(context,'v');
 const a=J.add(c(2),J.scale(v,mu)),s=J.div(J.sub(theta,v),a),b=J.add(c(1),J.scale(s,mu));
 const Delta=J.sub(J.scale(a,5),J.power(b,2)),S=J.add(c(1),s),V=J.add(c(1),v),T=J.add(c(10),theta);
 const P=J.sub(J.sub(J.sub(c(10),s),J.scale(v,5)),J.scale(J.power(s,2),mu.div(2)));
 const potential=J.sub(J.add(J.sub(J.sub(c(-10),theta),J.scale(v,10)),J.scale(J.power(v,2),'5/2')),J.div(J.power(J.sub(theta,v),2),J.scale(a,2)));
 const H=J.matrix([[a,b],[b,c(5)]]);
 const Q=J.matrix([[J.scale(J.div(b,a),-1),J.inverse(a)],[c(1),c(0)]]);
 const metric=J.matrix([[J.div(Delta,a),c(0)],[c(0),J.inverse(a)]]);
 const response=J.matrix([[J.scale(J.div(Delta,a),-1),J.scale(J.div(b,a),-1)],[J.scale(J.div(b,a),-1),J.inverse(a)]]);
 return {mu,context,v,theta,s,a,b,Delta,S,V,T,P,potential,H,Q,metric,response,reading:J.scale(V,5)};
}
function leviCivita(g){
 const depth=g.context.order-1;if(depth<0)throw new RangeError('One metric layer required');
 const inv=J.matrixProject(inverse2(g),depth),dg=['s','v'].map(x=>J.matrixDerivative(g,x));
 return [0,1].map(i=>J.matrix([0,1].map(k=>[0,1].map(j=>J.scale(sum([0,1].map(r=>J.mul(inv.rows[k][r],J.sub(J.add(dg[i].rows[r][j],dg[j].rows[r][i]),dg[r].rows[i][j])))),'1/2')))));
}
function curvature(gamma){
 const depth=gamma[0].context.order-1;
 const comm=J.matrixSub(J.matrixMul(gamma[0],gamma[1]),J.matrixMul(gamma[1],gamma[0]));
 return J.matrixAdd(J.matrixSub(J.matrixDerivative(gamma[1],'s'),J.matrixDerivative(gamma[0],'v')),J.matrixProject(comm,depth));
}
function transportedConnection(z,includeDerivative=true){
 const C=z.context,c=x=>J.scalar(C,x),m=z.mu,invH=inverse2(z.H),depth=C.order-1;
 const old=[J.matrix([[c(0),c(m)],[c(m),c(0)]]),J.matrix([[c(m),c(0)],[c(0),c(0)]])].map(h=>J.matrixScale(J.matrixMul(invH,h),'1/2'));
 const Q=J.matrixProject(z.Q,depth),Qi=J.matrixProject(inverse2(z.Q),depth);
 return [0,1].map(i=>{
  const contracted=J.matrix(old[0].rows.map((row,k)=>row.map((_,j)=>sum([0,1].map(r=>J.mul(old[r].rows[k][j],z.Q.rows[r][i]))))));
  let value=J.matrixMul(J.matrixMul(Qi,J.matrixProject(contracted,depth)),Q);
  if(includeDerivative)value=J.matrixAdd(value,J.matrixMul(Qi,J.matrixDerivative(z.Q,i===0?'s':'v')));
  return value;
 });
}
function transformedCurvature(z){
 const c=x=>J.scalar(z.context,x),f=J.div(c(z.mu.pow(2).div(4)),J.power(z.Delta,2));
 const old=J.matrix([[J.mul(f,J.scale(z.b,5)),J.scale(f,25)],[J.mul(f,J.scale(z.a,-5)),J.mul(f,J.scale(z.b,-5))]]);
 const result=J.matrixMul(J.matrixMul(inverse2(z.Q),old),z.Q),jac=J.scale(J.inverse(z.a),-1);
 return J.matrix(result.rows.map(row=>row.map(x=>J.mul(jac,x))));
}
function scalarCurvature(z,Fcurv){
 const depth=Fcurv.context.order,g=J.matrixProject(z.metric,depth);
 return J.div(J.matrixMul(g,Fcurv).rows[0][1],J.project(determinant(z.metric),depth));
}
function reflectionEqual(a,b,offset=0){
 if(a.context.order!==b.context.order)return false;
 return J.basisIndices(a.context).every(([x,y])=>J.coefficient(a,x,y).eq(J.coefficient(b,x,y).mul((x+y+offset)%2?-1:1)));
}
function matrixReflectionEqual(a,b){return a.height===b.height&&a.width===b.width&&a.rows.every((row,i)=>row.every((x,j)=>reflectionEqual(x,b.rows[i][j])));}
function grade(context,width=2){const signs=J.basisIndices(context).flatMap(([x,y])=>Array(width).fill((x+y)%2?-1:1));return N.matrix(signs.map((s,i)=>signs.map((_,j)=>i===j?s:0)));}
function exactPoint(mu,v,theta){
 mu=F.of(mu);v=F.of(v);theta=F.of(theta);const a=new F(2).add(mu.mul(v));if(a.le(0))throw new RangeError('Positive coordinate denominator required');
 const s=theta.sub(v).div(a),b=new F(1).add(mu.mul(s)),Delta=a.mul(5).sub(b.pow(2));
 if(Delta.le(0)||new F(1).add(s).le(0)||new F(1).add(v).le(0))throw new RangeError('Regular positive source chart required');
 return {mu,v,theta,s,a,b,Delta,S:new F(1).add(s),V:new F(1).add(v),P:new F(10).sub(s).sub(v.mul(5)).sub(mu.mul(s.pow(2)).div(2)),metric:N.matrix([[Delta.div(a),0],[0,new F(1).div(a)]]),response:N.matrix([[Delta.div(a).neg(),b.div(a).neg()],[b.div(a).neg(),new F(1).div(a)]]),reading:new F(1).add(v).mul(5),kappa:mu.pow(2).mul(5).div(Delta.pow(2).mul(4))};
}
module.exports=Object.freeze({inverse2,determinant,source,leviCivita,curvature,transportedConnection,transformedCurvature,scalarCurvature,reflectionEqual,matrixReflectionEqual,grade,exactPoint});
