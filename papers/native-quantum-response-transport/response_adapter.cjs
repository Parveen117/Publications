'use strict';
/** QRT-1 exact finite spectral comparison adapter. The analytic continuum
 * theorem is separate. Arithmetic is exclusively the frozen native Cut,
 * Paninian-backed TLM jets and TET path transport; no numerical integration.
 */
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const L=require('../native-tvsp-lambda-matrices/tvsp_adapter.cjs');
const T=require('../native-tvsp-exponential-transport/transport_adapter.cjs');
const {F,Cut,ZERO,ONE}=N;
const real=x=>{const c=Cut.of(x);if(!c.turn.zero())throw new TypeError('Real exact rational data required');return c;};
function positive(x,name){const c=real(x);if(c.rad.le(0))throw new RangeError(name+' must be strictly positive');return c;}
function pointsOf(points){if(!Array.isArray(points)||points.length<2)throw new TypeError('At least two spectral atoms required');return Object.freeze(points.map(p=>{if(!Array.isArray(p)||p.length!==2)throw new TypeError('Two observable values per atom required');return Object.freeze(p.map(x=>positive(x,'Spectral value')));}));}
function expNilpotent(x){if(!L.outer(x).zero())throw new RangeError('Zero-constant exponential argument required');let out=L.scalar(x.context,1),term=out;for(let k=1;k<=x.context.order;k++){term=L.scale(L.mul(term,x),new F(1,k));out=L.add(out,term);}return out;}
function logOne(x){if(!L.outer(x).eq(1))throw new RangeError('Logarithm branch requires constant one');const z=L.sub(x,L.scalar(x.context,1));let out=L.scalar(x.context),term=L.scalar(x.context,1);for(let k=1;k<=x.context.order;k++){term=L.mul(term,z);out=L.add(out,L.scale(term,new F(k%2?1:-1,k)));}return out;}
function powerOne(x,exponent){if(!L.outer(x).eq(1))throw new RangeError('Binomial branch requires constant one');const alpha=F.of(exponent),z=L.sub(x,L.scalar(x.context,1));let out=L.scalar(x.context,1),term=out,binomial=new F(1);for(let k=1;k<=x.context.order;k++){binomial=binomial.mul(alpha.sub(k-1)).div(k);term=L.mul(term,z);out=L.add(out,L.scale(term,binomial));}return out;}
function finiteSource(amplitudes,points,order=4,label='finite positive spectral comparison'){
  if(!Number.isSafeInteger(order)||order<1)throw new RangeError('Response source depth at least one required');
  const atoms=pointsOf(points);if(!Array.isArray(amplitudes)||amplitudes.length!==atoms.length)throw new TypeError('One amplitude per atom required');
  const amps=Object.freeze(amplitudes.map(x=>positive(x,'Amplitude'))),weights=Object.freeze(amps.map(a=>a.mul(a)));
  if(!weights.reduce((s,p)=>s.add(p),ZERO).eq(1))throw new RangeError('Squared amplitudes must sum exactly to one');
  const high=L.createContext(order+2,label),s=L.coordinate(high,'s'),v=L.coordinate(high,'v');
  let Z=L.scalar(high);const phiHigh=[];
  for(let k=0;k<atoms.length;k++){
    const exponent=L.scale(L.add(L.scale(s,atoms[k][0]),L.scale(v,atoms[k][1])),-1);
    Z=L.add(Z,L.scale(expNilpotent(exponent),weights[k]));
    phiHigh.push(L.scale(expNilpotent(L.scale(exponent,'1/2')),amps[k]));
  }
  const u=logOne(Z),duS=L.derivative(u,'s'),duV=L.derivative(u,'v');
  const a=L.derivative(duS,'s'),b=L.derivative(duS,'v'),c=L.derivative(duV,'v'),context=a.context;
  const means=Object.freeze([L.scale(L.project(duS,order),-1),L.scale(L.project(duV,order),-1)]);
  const normalizer=powerOne(Z,new F(-1,2));
  const psi=L.matrix(phiHigh.map(phi=>[L.project(L.mul(normalizer,phi),order)]));
  const E=L.matrix(atoms.map((atom,k)=>atom.map((value,i)=>L.mul(L.sub(L.scalar(context,value),means[i]),psi.rows[k][0]))));
  const H=L.matrix([[a,b],[b,c]]),Delta=L.sub(L.mul(a,c),L.mul(b,b));
  return Object.freeze({kind:'QRT_FINITE_SPECTRAL_COMPARISON',context,atoms,amplitudes:amps,weights,means,
    partition:L.project(Z,order),potential:L.project(u,order),psi,E,H,a,b,c,Delta,
    unnormalized:L.matrix(phiHigh.map(phi=>[L.project(phi,order)])),
    precursor:Object.freeze({context:high,partition:Z,potential:u})});
}
function responseGeometry(source){
  const z=source,inv=L.matrix([[L.div(z.c,z.Delta),L.scale(L.div(z.b,z.Delta),-1)],[L.scale(L.div(z.b,z.Delta),-1),L.div(z.a,z.Delta)]]);
  const Q=L.matrixMul(z.E,L.matrixMul(inv,L.matrixDagger(z.E)));
  const selected=T.sourceGeometry(z);
  return Object.freeze({Hinv:inv,Q,...selected});
}
function responsePath(source,path={a:1,b:1,p:1,q:2}){
  const geometry=responseGeometry(source),H=T.restrictMatrix(source.H,path),E=T.restrictMatrix(source.E,path),Q=T.restrictMatrix(geometry.Q,path);
  const A=T.connectionFromMetric(H),U=T.transport(A),depth=Q.context.order-1,Qprime=T.derivative(Q),q=T.project(Q,depth);
  const Kgen=T.sub(T.mul(Qprime,q),T.mul(q,Qprime)),K=T.transport(T.scale(Kgen,-1));
  return Object.freeze({path,H,E,Q,A,U,Kgen,K,Qprime,geometry});
}
function determinant3(M){const [[a,b,c],[d,e,f],[g,h,i]]=N.matrix(M);return a.mul(e.mul(i).sub(f.mul(h))).sub(b.mul(d.mul(i).sub(f.mul(g)))).add(c.mul(d.mul(h).sub(e.mul(g))));}
function moments(probabilities,points){
  const atoms=pointsOf(points);if(!Array.isArray(probabilities)||probabilities.length!==atoms.length)throw new TypeError('One probability per atom required');
  const weights=Object.freeze(probabilities.map(x=>positive(x,'Probability')));if(!weights.reduce((s,p)=>s.add(p),ZERO).eq(1))throw new RangeError('Probabilities must sum exactly to one');
  const means=[0,1].map(i=>atoms.reduce((s,x,k)=>s.add(weights[k].mul(x[i])),ZERO));
  const centered=atoms.map(x=>x.map((value,i)=>value.sub(means[i])));
  const raw=N.matrix([0,1].map(i=>[0,1].map(j=>atoms.reduce((s,x,k)=>s.add(weights[k].mul(x[i]).mul(x[j])),ZERO))));
  const H=N.matrix([0,1].map(i=>[0,1].map(j=>centered.reduce((s,x,k)=>s.add(weights[k].mul(x[i]).mul(x[j])),ZERO))));
  const third=[[0,0,0],[0,0,1],[0,1,1],[1,1,1]].map(indices=>centered.reduce((s,x,k)=>s.add(indices.reduce((v,i)=>v.mul(x[i]),weights[k])),ZERO));
  const u3=third.map(x=>x.neg()),[[a,b],[,c]]=H,Delta=a.mul(c).sub(b.mul(b));
  const D=determinant3([[a,b,c],[u3[0],u3[1],u3[2]],[u3[1],u3[2],u3[3]]]);
  const curvature=Delta.zero()?null:D.neg().div(Delta.mul(Delta).mul(4));
  return Object.freeze({atoms,weights,means:Object.freeze(means),H,raw,third:Object.freeze(third),u3:Object.freeze(u3),Delta,D,curvature});
}
function continuumBounds(delta=new F(1,10000)){
  delta=F.of(delta);if(delta.le(0))throw new RangeError('Positive window half-width required');
  const covarianceError=delta.mul(4),thirdError=delta.mul(6),delta0=new F(1,27),determinantMagnitude=new F(1,729),determinantError=delta.mul(new F(5,3));
  const deltaLower=delta0.sub(covarianceError),deltaUpper=delta0.add(covarianceError),determinantLower=determinantMagnitude.sub(determinantError),determinantUpper=determinantMagnitude.add(determinantError);
  if(deltaLower.le(0)||determinantLower.le(0))throw new RangeError('Requested bound does not certify positive determinant margins');
  return Object.freeze({delta,covarianceError,thirdError,delta0,deltaLower,deltaUpper,determinantMagnitude,determinantError,determinantLower,determinantUpper,
    curvatureLower:determinantLower.div(deltaUpper.pow(2).mul(4)),curvatureUpper:determinantUpper.div(deltaLower.pow(2).mul(4))});
}
module.exports=Object.freeze({finiteSource,responseGeometry,responsePath,moments,determinant3,continuumBounds,expNilpotent,logOne,powerOne});
