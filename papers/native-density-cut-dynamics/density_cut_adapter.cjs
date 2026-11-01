'use strict';
/** DCD-1: a source-specific adapter over frozen QRT, TET and native matrices.
 * Source derivatives consume jet depth. Pulse series are exact finite Taylor
 * coefficients in a separate parameter, not numerical exponentials or limits.
 */
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const T=require('../native-tvsp-exponential-transport/transport_adapter.cjs');
const QRT=require('../native-quantum-response-transport/response_adapter.cjs');
const comm=(a,b)=>T.sub(T.mul(a,b),T.mul(b,a));
function densityCutPath(source,path={a:1,b:1,p:1,q:2}){
  const qrt=QRT.responsePath(source,path),C=qrt.Q.context,order=C.order,p=C.path,m=source.atoms.length;
  const psi=T.restrictMatrix(source.psi,path),Pi=T.mul(psi,T.dagger(psi)),W=T.add(Pi,T.scale(qrt.Q,'1/2'));
  const d=Array.from({length:order+1},()=>N.zeros(m));
  if(p.p-1<=order)for(let k=0;k<m;k++)d[p.p-1][k][k]=d[p.p-1][k][k].add(source.atoms[k][0].mul(p.a).mul(p.p));
  if(p.q-1<=order)for(let k=0;k<m;k++)d[p.q-1][k][k]=d[p.q-1][k][k].add(source.atoms[k][1].mul(p.b).mul(p.q));
  const D=T.series(C,d),remainder=T.sub(T.sub(T.identity(C,m),Pi),qrt.Q),depth=order-1;
  const B=T.project(comm(W,D),depth),Bprojection=qrt.Kgen,K=T.transport(T.scale(B,-1));
  return Object.freeze({...qrt,psi,Pi,W,D,remainder,B,Bprojection,K,Kprojection:qrt.K,
    contract:'D=dot(gamma) dot X; W=Pi+Q/2; all derivative identities at depth N-1'});
}
function weightedGenerator(pathData,alpha,beta){
  const W=T.add(T.scale(pathData.Pi,N.Cut.of(alpha)),T.scale(pathData.Q,N.Cut.of(beta)));
  return T.project(comm(W,pathData.D),pathData.B.context.order);
}
function pulseJet(W,D,order=3,orientation=1){
  W=N.matrix(W);D=N.matrix(D);
  if(W.length!==W[0].length||D.length!==W.length||D[0].length!==W.length||!N.equal(W,N.dagger(W))||!N.equal(D,N.dagger(D)))throw new TypeError('Matched Hermitian matrices required for coherent pulses');
  if(!Number.isSafeInteger(order)||order<2)throw new RangeError('Pulse order at least two required');
  if(orientation!==1&&orientation!==-1)throw new RangeError('Pulse orientation must be +1 or -1');
  const C=T.createContext(order,'DCD independent pulse parameter'),m=W.length;
  // Rightmost acts first: G=exp(-iW)exp(+iD)exp(+iW)exp(-iD).
  const pair=orientation===1?[W,D]:[D,W],generators=[N.scale(pair[1],N.IOTA.neg()),N.scale(pair[0],N.IOTA),N.scale(pair[1],N.IOTA),N.scale(pair[0],N.IOTA.neg())];
  const pulses=generators.map(A=>T.ordinaryExponential(T.series(C,Array.from({length:order+1},(_,k)=>k===1?A:N.zeros(m)))));
  let product=T.identity(C,m);for(const pulse of pulses)product=T.mul(pulse,product);
  const target=N.scale(N.commutator(W,D),orientation);
  return Object.freeze({context:C,orientation,pulses:Object.freeze(pulses),product,target,
    chronological:orientation===1?'-iD,+iW,+iD,-iW':'-iW,+iD,+iW,-iD'});
}
function pulseResourceBound({delta,w=1,d,b,steps=1}){
  delta=N.F.of(delta);w=N.F.of(w);d=N.F.of(d);b=b===undefined?w.mul(d).mul(2):N.F.of(b);
  if(delta.le(0)||w.le(0)||d.le(0)||b.le(0)||!Number.isSafeInteger(steps)||steps<1)throw new RangeError('Positive rational resources and positive integer steps required');
  const cubic=delta.pow(3).mul(w.add(d)).mul(b),quartic=delta.pow(4).mul(b.pow(2)).div(2),local=cubic.add(quartic);
  return Object.freeze({delta,w,d,b,steps,cubic,quartic,local,total:local.mul(steps),effectiveDuration:delta.pow(2).mul(steps),
    contract:'conditional norm bound: ||G_delta-exp(delta^2 B)|| <= delta^3(w+d)b+delta^4 b^2/2; fixed generators'});
}
module.exports=Object.freeze({densityCutPath,weightedGenerator,pulseJet,pulseResourceBound});
