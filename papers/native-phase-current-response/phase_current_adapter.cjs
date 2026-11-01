'use strict';
/** PCR-1 adapter. Finite matrices represent commuting densities and states.
 * Current operators are NEVER represented by finite matrices; their affine
 * commutators and mean-readout formulas belong to the separate symbolic lane.
 */
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const T=require('../native-tvsp-exponential-transport/transport_adapter.cjs');
const DCD=require('../native-density-cut-dynamics/density_cut_adapter.cjs');
const {F,Cut,ZERO}=N;
function real(x){const c=Cut.of(x);if(!c.turn.zero())throw new TypeError('Exact real rational required');return c;}
function hValue(h){h=F.of(h);if(h.le(0)||new F(1).le(h))throw new RangeError('Phase normalization h must satisfy 0<h<1');return h;}
function freezeMatrix(a){return Object.freeze(N.matrix(a).map(r=>Object.freeze(r)));}
function blocks(a,b,c,d){return freezeMatrix(a.map((r,i)=>r.concat(b[i])).concat(c.map((r,i)=>r.concat(d[i]))));}
function covariance(H){H=N.matrix(H);if(H.length!==2||H[0].length!==2||H.flat().some(x=>!x.turn.zero())||!N.equal(H,N.dagger(H)))throw new TypeError('Real symmetric two-dimensional covariance required');const Delta=H[0][0].mul(H[1][1]).sub(H[0][1].mul(H[1][0]));if(H[0][0].rad.le(0)||Delta.rad.le(0))throw new RangeError('Positive definite covariance required');return H;}
function quantumGeometry(H,h='1/2'){
  H=covariance(H);h=hValue(h);const zero=N.zeros(2),I=N.identity(2),ih=N.IOTA.mul(new F(1).div(h.mul(2)));
  const QGT=blocks(N.scale(H,'1/4'),N.scale(H,ih),N.scale(H,ih.neg()),N.scale(H,new F(1).div(h.pow(2))));
  const g=blocks(N.scale(H,'1/4'),zero,zero,N.scale(H,new F(1).div(h.pow(2))));
  const J=blocks(zero,N.scale(I,new F(-2).div(h)),N.scale(I,h.div(2)),zero);
  const BerryF=blocks(zero,N.scale(H,new F(-1).div(h)),N.scale(H,new F(1).div(h)),zero);
  const omega=freezeMatrix(N.mul(N.dagger(J),g));
  const factor=freezeMatrix(N.scale(I,'-1/2').map((row,i)=>row.concat(N.scale(I,N.IOTA.mul(new F(-1).div(h)))[i])));
  return Object.freeze({H:freezeMatrix(H),h,QGT,g,J,BerryF,omega,factor});
}
function readoutGeometry(means,H,beta,h='1/2'){
  H=covariance(H);h=hValue(h);if(!Array.isArray(means)||means.length!==2||!Array.isArray(beta)||beta.length!==2)throw new TypeError('Two means and two phase parameters required');
  const m=means.map(real),b=beta.map(real);if(m.some(x=>x.rad.le(0)))throw new RangeError('Positive means are required for logarithmic coordinates');
  const diag=values=>N.matrix([[values[0],0],[0,values[1]]]),zero=N.zeros(2),I=N.identity(2),dm=diag(m),db=diag(b),dminv=diag(m.map(x=>x.inv()));
  const p=Object.freeze(m.map((x,i)=>x.mul(b[i]))),qTheta=N.scale(N.mul(dminv,H),-1),pTheta=N.scale(N.mul(db,H),-1);
  const qpJacobian=blocks(qTheta,zero,pTheta,dm),mpJacobian=blocks(N.scale(H,-1),zero,pTheta,dm);
  const canonical=blocks(zero,N.scale(I,new F(1).div(h)),N.scale(I,new F(-1).div(h)),zero);
  const pullback=freezeMatrix(N.mul(N.dagger(qpJacobian),N.mul(canonical,qpJacobian)));
  const Delta=H[0][0].mul(H[1][1]).sub(H[0][1].mul(H[1][0]));
  return Object.freeze({means:Object.freeze(m),beta:Object.freeze(b),p,betaRecovered:Object.freeze(p.map((x,i)=>x.div(m[i]))),
    qContract:'q_i=log(m_i); no finite rational logarithm value is asserted',qpJacobian,mpJacobian,canonical,pullback,
    mpDeterminant:Delta.mul(m[0]).mul(m[1]),qpDeterminant:Delta});
}
function scalarIdentity(x,n){if(x.height!==1||x.width!==1)throw new TypeError('Scalar path record required');return T.series(x.context,T.coefficients(x).map(a=>N.scale(N.identity(n),a[0][0])));}
function phaseCurrentPath(source,{h='1/2',c='3/5',d='-2/7',thetaPath={a:1,b:1,p:1,q:2}}={}){
  h=hValue(h);c=real(c);d=real(d);const base=DCD.densityCutPath(source,thetaPath),C=base.K.context,order=C.order,m=source.atoms.length;
  const exponent=Array.from({length:order+1},()=>N.zeros(m)),velocity=Array.from({length:order+1},()=>N.zeros(m));
  for(let k=0;k<m;k++){
    exponent[1][k][k]=source.atoms[k][0].mul(c).mul(N.IOTA).mul(new F(-1).div(h));
    if(order>=2)exponent[2][k][k]=source.atoms[k][1].mul(d).mul(N.IOTA).mul(new F(-1).div(h));
    velocity[0][k][k]=source.atoms[k][0].mul(c);
    if(order>=1)velocity[1][k][k]=source.atoms[k][1].mul(d).mul(2);
  }
  const V=T.ordinaryExponential(T.series(C,exponent)),dotBetaX=T.series(C,velocity),psi=T.mul(V,base.psi),lift=T.mul(V,base.K),depth=order-1;
  const v=T.project(V,depth),liftGenerator=T.add(T.project(dotBetaX,depth),T.scale(T.mul(v,T.mul(base.B,T.dagger(v))),N.IOTA.mul(h)));
  const psiPrime=T.derivative(psi),psiLow=T.project(psi,depth),inner=T.mul(T.dagger(psiLow),psiPrime),BerryA=T.scale(inner,N.IOTA);
  const horizontal=T.sub(psiPrime,T.mul(psiLow,inner));
  const betaVelocity=T.series(C,Array.from({length:order+1},(_,k)=>N.matrix([[k===0?c:ZERO],[k===1?d.mul(2):ZERO]])));
  const path=C.path,thetaVelocity=T.series(C,Array.from({length:order+1},(_,k)=>N.matrix([[k===path.p-1?path.a.mul(path.p):ZERO],[k===path.q-1?path.b.mul(path.q):ZERO]])));
  const H=T.project(base.H,depth),thetaLow=T.project(thetaVelocity,depth),betaLow=T.project(betaVelocity,depth);
  const horizontalNormExpected=T.add(T.scale(T.mul(T.dagger(thetaLow),T.mul(H,thetaLow)),'1/4'),T.scale(T.mul(T.dagger(betaLow),T.mul(H,betaLow)),new F(1).div(h.pow(2))));
  const meanPaths=source.means.map(x=>T.restrictScalar(x,thetaPath));
  const betaPaths=[T.series(C,Array.from({length:order+1},(_,k)=>N.matrix([[k===1?c:ZERO]]))),T.series(C,Array.from({length:order+1},(_,k)=>N.matrix([[k===2?d:ZERO]])))];
  const pReadouts=meanPaths.map((x,i)=>T.mul(x,betaPaths[i]));
  const betaRecovered=pReadouts.map((x,i)=>T.mul(x,T.inverse(meanPaths[i])));
  const expectedBerry=T.scale(T.add(T.mul(T.project(meanPaths[0],depth),T.constant(BerryA.context,[[c]])),
    T.mul(T.project(meanPaths[1],depth),T.series(BerryA.context,Array.from({length:depth+1},(_,k)=>[[k===1?d.mul(2):ZERO]])))),new F(1).div(h));
  const scalarPhase=T.ordinaryExponential(T.scale(T.integrateZero(BerryA),N.IOTA)),horizontalState=T.mul(scalarIdentity(scalarPhase,m),psi);
  return Object.freeze({source,base,h,c,d,V,dotBetaX,psi,lift,liftGenerator,BerryA,expectedBerry,horizontal,horizontalNormExpected,
    thetaVelocity,betaVelocity,meanPaths:Object.freeze(meanPaths),betaPaths:Object.freeze(betaPaths),pReadouts:Object.freeze(pReadouts),betaRecovered:Object.freeze(betaRecovered),scalarPhase,horizontalState,
    currentContract:'p_i=beta_i m_i is the analytically specified current mean; no finite current matrix is supplied'});
}
module.exports=Object.freeze({quantumGeometry,readoutGeometry,phaseCurrentPath,scalarIdentity});
