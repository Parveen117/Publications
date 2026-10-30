'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const f=N.F.of,C=N.Cut.of,trace=a=>a.reduce((s,row,i)=>s.add(row[i]),C(0));
const I=N.identity(2),X=N.matrix([[0,1],[1,0]]),Y=N.matrix([[0,[0,-1]],[[0,1],0]]),Z=N.matrix([[1,0],[0,-1]]);
const words=[N.kron(I,X),N.kron(Z,X),N.kron(I,Y),N.kron(Z,Y),N.kron(X,I),N.kron(X,Z)];
function state(name){
 if(name==='vacuum'){const p=N.zeros(4);p[0][0]=C(1);return p;}
 if(name==='bell')return N.scale(N.matrix([[1,0,0,1],[0,0,0,0],[0,0,0,0],[1,0,0,1]]),'1/2');
 assert(['controlled_s','controlled_s_dagger'].includes(name));
 const a=N.matrix([[1],[1],[1],[[0,name==='controlled_s'?1:-1]]]);return N.scale(N.mul(a,N.dagger(a)),'1/4');
}
const means=(p,hs)=>hs.map(a=>trace(N.mul(p,a)));
function covariance(p,hs){const m=means(p,hs);return hs.map((a,i)=>hs.map((b,j)=>trace(N.mul(p,N.mul(a,b))).sub(m[i].mul(m[j]))));}
function controls(row){
 const a=f(row.alpha),b=f(row.beta),g=f(row.gamma),D=N.zeros(4),J=N.identity(4);
 D[1][0]=D[0][1]=C(a);J[1][0]=C([0,b]);J[0][1]=C([0,b.neg()]);J[2][0]=J[0][2]=C(g);return[D,J];
}
function weights(row){const a=f(row.alpha).div(2),b=f(row.beta).div(2),g=f(row.gamma).div(2);return N.matrix([[a,a,0,0,0,0],[0,0,b,b,g,g]]);}
function scalar(entry){const[a,b,c,d]=entry.coefficients;assert(f(b).zero()&&f(d).zero());return C([a,c]);}
function parse(report){assert.equal(report.status,'completed');const a=report.answers[0];assert.equal(a.algebraically_exact,true);return{G:N.matrix(a.geometry.tensor.map(r=>r.map(scalar))),m:N.matrix(a.geometry.means.map(x=>[scalar(x.exact)]))};}
const dir=process.argv[2],source=JSON.parse(fs.readFileSync(path.join(dir,'source_intervals.json')));
const checks=[],negatives=[],check=(name,v)=>{assert(v,name);checks.push(name);},reject=(name,v)=>{assert(v,name);negatives.push(name);};
const results={};
for(const name of ['vacuum','bell','controlled_s','controlled_s_dagger'])for(const backend of ['observer','clifford_t']){
 const key=name+'_'+backend,p=state(name),r=parse(JSON.parse(fs.readFileSync(path.join(dir,key+'.result.json'))));results[key]=r;
 check(key+':independent_Pauli_covariance',N.equal(r.G,covariance(p,words)));
 check(key+':independent_Pauli_means',N.equal(r.m,N.matrix(means(p,words).map(x=>[x]))));
 for(const row of source.cases){
  const W=weights(row),Q=N.mul(N.mul(W,r.G),N.dagger(W)),Hs=controls(row);
  check(key+':weighted_covariance:'+row.h,N.equal(Q,covariance(p,Hs)));
  const mm=N.mul(W,r.m);mm[1][0]=mm[1][0].add(1);
  check(key+':weighted_means:'+row.h,N.equal(mm,N.matrix(means(p,Hs).map(x=>[x]))));
  for(const ts of source.times){
   const t=f(ts),S=N.matrix([[1,t],[0,1]]),pulseHs=[N.add(Hs[0],N.scale(Hs[1],t)),Hs[1]];
   check(key+':shear_covariance:h='+row.h+':t='+ts,N.equal(N.mul(N.mul(S,Q),N.dagger(S)),covariance(p,pulseHs)));
   check(key+':shear_means:h='+row.h+':t='+ts,N.equal(N.mul(S,mm),N.matrix(means(p,pulseHs).map(x=>[x]))));
  }
 }
}
const comparisons=[];
for(const row of source.cases){
 const h=f(row.h),W=weights(row),Q=N.mul(N.mul(W,results.vacuum_observer.G),N.dagger(W));
 const A=f(row.A),c=f(row.c),lo=f(row.B.lower),hi=f(row.B.upper),tol=f('1/1000000');
 check('source_A:'+row.h,Q[0][0].rad.sub(A).abs().le(tol));
 check('source_B_enclosure:'+row.h,Q[1][1].rad.sub(lo).abs().le(tol)&&Q[1][1].rad.sub(hi).abs().le(tol));
 check('source_cross_sign:'+row.h,Q[0][1].turn.sub(c).abs().le(tol));
 check('source_positive_Gram:'+row.h,N.rank(Q)===2&&!f(row.determinant_lower).le(0));
 for(const ts of source.times){
  const t=f(ts),S=N.matrix([[1,t],[0,1]]),qt=N.mul(N.mul(S,Q),N.dagger(S));
  const measured=qt[0][0].rad.add(h.div(4)),exactLo=h.div(2).add(t.pow(2).mul(lo)),exactHi=h.div(2).add(t.pow(2).mul(hi));
  const err=f(row.absolute_error_bounds.A).add(t.pow(2).mul(row.absolute_error_bounds.B));
  check('actual_variance_error_bound:h='+row.h+':t='+ts,measured.sub(exactLo).abs().le(err)&&measured.sub(exactHi).abs().le(err));
  const naiveLo=h.div(2).add(t.pow(2).mul(f(row.I2.lower).sub(1))),naiveHi=h.div(2).add(t.pow(2).mul(f(row.I2.upper).sub(1)));
  if(!t.zero())reject('naive_density_shift_variance:h='+row.h+':t='+ts,!naiveLo.sub(exactHi).le(0)&&!naiveLo.sub(measured).le(err));
  else check('zero_pulse_common_variance:'+row.h,exactLo.eq(naiveLo)&&exactHi.eq(naiveHi));
  comparisons.push({h:row.h,t:ts,encoded_measured_variance:measured,certified_error_bound:err,actual:{lower:exactLo,upper:exactHi},naive:{lower:naiveLo,upper:naiveHi}});
 }
 const hs=controls(row),J=hs[1],D=hs[0],comm=N.commutator(D,J),affine=N.scale(J,C([0,h]));
 reject('finite_encoding_not_affine_representation:'+row.h,!N.equal(comm,affine)&&trace(comm).zero()&&!trace(affine).zero());
 reject('missing_measurement_noise:'+row.h,!h.div(4).le(tol));
 const bad=controls({...row,beta:f(row.beta).neg().toString()});
 reject('reversed_imaginary_cross:'+row.h,!N.equal(covariance(state('vacuum'),bad),Q));
}
reject('phase_dagger_state_is_observable',!N.equal(results.controlled_s_observer.G,results.controlled_s_dagger_observer.G));
for(const row of source.bounded_response)reject('naive_bounded_sine_response:'+row.h,!f(row.sine_response_slope_gap.lower).le(0));
for(const block of source.cutoff_moments){
 const h=f(block.h),rows=block.rows;
 for(let i=0;i<rows.length;i++){
  const row=rows[i];
  check('cutoff_positive_error:h='+block.h+':N='+row.N,!f(row.mass_L2_error_squared.lower).le(0));
  if(i)check('cutoff_strict_progress:h='+block.h+':N='+row.N,!f(row.I_N.lower).sub(rows[i-1].I_N.upper).le(0));
 }
 // Independent N=1 check by the positive Bessel I_0(h) scalar series,
 // rather than the Laurent-polynomial convolution used by the producer.
 let term=f(1),sum=f(1);
 for(let j=1;j<=40;j++){term=term.mul(h.pow(2)).div(4*j*j);sum=sum.add(term);}
 const tail=term.mul(h.pow(2)).div(4*41*41).div(f(1).sub(h.pow(2).div(4*42*42)));
 const low=f(rows[0].I_N.lower).mul(sum),high=f(rows[0].I_N.upper).mul(sum.add(tail));
 check('cutoff_N1_independent_Bessel_overlap:'+block.h,low.le(rows[1].I_N.upper)&&f(rows[1].I_N.lower).le(high));
}
// Exact truncated formal series check, independently expanding BCH translation
// and Gaussian overlap: coefficient of i*k^3 is t*h^2/24-t*h/4-t^3*m3/6.
const P=(...xs)=>xs.map(C),pad=a=>Array.from({length:4},(_,i)=>a[i]||C(0));
const pmul=(a,b)=>Array.from({length:4},(_,n)=>pad(a).reduce((s,x,j)=>j<=n?s.add(x.mul(pad(b)[n-j])):s,C(0)));
for(const hs of ['1/8','1/4','3/8'])for(const ts of ['-1','1/2','1']){
 const h=f(hs),t=f(ts),m2=f('7/3'),m3=f(9),a=P(0,t,0,t.mul(h.pow(2)).div(24));
 const a2=pmul(a,a),a3=pmul(a2,a),phi=P(1);
 for(let j=0;j<4;j++)phi[j]=(phi[j]||C(0)).add(a[j].mul([0,1])).sub(a2[j].mul(m2.div(2))).add(a3[j].mul([0,m3.neg().div(6)]));
 const out=pmul(P(1,0,h.neg().div(4),0),phi);
 check('bounded_characteristic_jet:h='+hs+':t='+ts,out[1].eq([0,t])&&out[2].eq(h.neg().div(4).sub(t.pow(2).mul(m2).div(2)))&&out[3].eq([0,t.mul(h.pow(2)).div(24).sub(t.mul(h).div(4)).sub(t.pow(3).mul(m3).div(6))]));
}
process.stdout.write(JSON.stringify({status:'PASS',checks:checks.length,negative_controls:negatives.length,
 check_names:checks,negative_control_names:negatives,source_values:source.cases.map(r=>({h:r.h,...r.display})),
 bounded_response:source.bounded_response,cutoff_moments:source.cutoff_moments,comparisons,
 scope:'Eight fresh NRQC finite executions; three conditional source points and five pulse times; four states, two engines. Only vacuum encodes source moments.',
 continuum_pulse_simulated:false,hardware_called:false,advantage_established:false,
 analytic_advancement:'Exact characteristic function and bounded sine linear-response witness, proved separately; formal finite jet check is not a proof of the continuum result.'},null,2)+'\n');
