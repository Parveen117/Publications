'use strict';
const assert=require('node:assert/strict'),R=require('./source.cjs');
const {f,iv,point,add,sub,neg,mul,scale,div,exp,sqrt,PI,trig,widen,absmax,contains,overlap,
 gamma,gammaCDF,bessel0,bessel1,halfGammaTails,radialWeighted,angular0}=R;

function angularQuadrature(x,M=64){
 x=f(x);let sum=point(0);const step=scale(PI,f(2).div(M));
 for(let j=0;j<M;j++){
  const theta=scale(PI,f(2*j+1).div(M));
  sum=add(sum,trig(scale(trig(theta).cos,x)).cos);
 }
 return widen(scale(sum,f(1).div(M)),absmax(mul(step,step)).mul(x.abs().add(x.pow(2))).div(24));
}

function verify(){
 const checks=[],negatives=[],check=(n,ok)=>{assert(ok,n);checks.push(n);},reject=(n,ok)=>{assert(ok,n);negatives.push(n);};
 check('half_Gamma_independent_normalizer',overlap(gamma('1/2'),sqrt(PI)));
 check('three_half_Gamma_recurrence',overlap(gamma('3/2'),scale(gamma('1/2'),'1/2')));
 for(const x of ['1/2','2','8']){
  const tails=halfGammaTails(x),q=sub(point(1),gammaCDF('3/2',x));
  check('half_Gamma_recurrence:'+x,overlap(tails[0],q));
  check('half_Gamma_shape_monotone:'+x,tails.every((t,j)=>j===0||tails[j-1].lo.le(t.hi)));
  check('ordinary_Bessel_vs_angular_integral:'+x,overlap(angular0(point(x)).value,angularQuadrature(x)));
 }
 check('ordinary_Bessel_zero',contains(angular0(point(0)).value,1));
 check('ordinary_Bessel_even',overlap(angular0(point(2)).value,angular0(point(-2)).value));
 check('large_argument_fallback_kept',angular0(point(33)).fallback&&contains(angular0(point(33)).value,-1)&&contains(angular0(point(33)).value,1));
 check('modified_Bessel_zero',contains(bessel1(0),0));
 check('half_Gamma_zero_endpoint',halfGammaTails(0).every(t=>contains(t,1)));
 check('half_Gamma_infinite_endpoint',halfGammaTails(Infinity).every(t=>contains(t,0)));
 const h=f('1/4'),coarse=radialWeighted(h,16),fine=radialWeighted(h,32);
 check('positive_radial_mass',fine.rows.every(r=>r.first1.lo.n>0n&&r.first1.hi.le(r.first.lo)));
 check('radial_probability_and_uniform_mean',contains(fine.totalP,1)&&contains(fine.totalA,1));
 check('integrated_B1_mixture_normalization',overlap(fine.first1Total,fine.coefficientTotal));
 check('nested_B1_second_increases',coarse.second1.hi.le(fine.second1.lo));
 check('B1_second_below_continuous',fine.second1.hi.le(bessel1(h).lo));
 check('unbounded_radial_tail_retained',fine.rows.at(-1).right===Infinity&&fine.rows.at(-1).first1.lo.n>0n);
 reject('discarding_radial_tail_loses_probability',!contains(sub(fine.totalP,fine.rows.at(-1).p),1));
 reject('fixed_aligned_angle_changes_mean',!contains(add(point(1),scale(fine.first1Total,'1/2')),1));
 reject('dropping_relative_angle_loses_second_moment',bessel1(h).lo.n>0n);
 reject('ordinary_is_not_modified_Bessel',!overlap(angular0(point('1/2')).value,bessel0('1/2')));
 reject('angular_factor_is_not_one',!contains(angular0(point('1/2')).value,1));
 reject('half_shape_is_not_integer_shape',!overlap(halfGammaTails(2)[0],exp(point(-2))));
 const factor=f(1).add(f('1/4').mul('3/8').div(f(2).mul(f(2).sub('3/8'))));
 check('weighted_second_exact_witness',factor.eq('107/104'));
 reject('uniform_mass_law_changes_weighted_second',!factor.eq(1));
 assert.throws(()=>radialWeighted(h,24));negatives.push('unfrozen_radial_partition_rejected');
 assert.throws(()=>halfGammaTails(17));negatives.push('unsupported_radial_endpoint_rejected');
 return{checks,negatives};
}
if(require.main===module){const r=verify();console.log(JSON.stringify({status:'PASS',checks:r.checks.length,negative_controls:r.negatives.length}));}
module.exports={verify,angularQuadrature};
