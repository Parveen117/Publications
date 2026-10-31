'use strict';
const assert=require('node:assert/strict');
const R=require('./source.cjs');
const {f,point,add,sub,mul,scale,div,exp,overlap,contains,gamma,gammaLower,gammaCDF,dyadicPower,exponentialPartition}=R;
function verify(){
 const checks=[],negatives=[],check=(name,ok)=>{assert(ok,name);checks.push(name);};
 const reject=(name,ok)=>{assert(ok,name);negatives.push(name);};
 check('Gamma_one',contains(gamma(1),1));
 check('Gamma_two',contains(gamma(2),1));
 for(const s of ['1/2','7/8','13/16','15/16']){
  check('Gamma_recurrence:'+s,overlap(gamma(f(s).add(1)),scale(gamma(s),s)));
  const x=f('3/2');
  check('incomplete_Gamma_recurrence:'+s,overlap(gammaLower(f(s).add(1),x),
    sub(scale(gammaLower(s,x),s),mul(dyadicPower(x,s),exp(point(x.neg()))))));
 }
 for(const x of ['1/65536','1/2','2','16']){
  check('exponential_CDF:'+x,overlap(gammaCDF(1,x),sub(point(1),exp(point(f(x).neg())))));
 }
 check('Gamma_CDF_zero',contains(gammaCDF('7/8',0),0));
 check('Gamma_CDF_infinity',contains(gammaCDF('7/8',Infinity),1));
 const g=f('1/4'),coarse=exponentialPartition(g,2),fine=exponentialPartition(g,4);
 check('q_independent_mass_normalizer',coarse.normalizer.lo.n>0n);
 check('positive_conditional_mass',fine.rows.every(r=>r.mass.lo.n>0n));
 check('all_tails_normalized',contains(fine.totalP,1)&&contains(fine.totalA,1));
 check('nested_mass_second_increases',coarse.second.hi.le(fine.second.lo));
 check('first_cell_keeps_singularity',fine.rows[0].left.zero()&&fine.rows[0].right.eq('1/65536'));
 check('unbounded_tail_retained',fine.rows.at(-1).right===Infinity);
 check('fine_bin_count',fine.rows.length===82);
 check('coarse_bin_count',coarse.rows.length===42);
 const left=coarse.rows[0],right=coarse.rows.at(-1);
 reject('discarding_singular_end_loses_probability',!contains(sub(coarse.totalP,left.p),1));
 reject('discarding_large_E_tail_loses_probability',!contains(sub(coarse.totalP,right.p),1));
 reject('omitting_Gamma_normalizer_changes_mean',!contains(mul(coarse.totalA,coarse.normalizer),1));
 reject('omitting_zero_mode_Wick_subtraction_changes_mean',!contains(exp(point(g.div(2))),1));
 const wrong=div(gamma(f(1).sub(g.mul(2))),mul(gamma(f(1).sub(g)),gamma(f(1).sub(g))));
 reject('wrong_coupling_exponent_changes_second_moment',!overlap(wrong,coarse.fullSecond));
 assert.throws(()=>gamma(0));negatives.push('nonpositive_Gamma_shape_rejected');
 assert.throws(()=>dyadicPower(2,'2/3'));negatives.push('unsupported_fractional_power_rejected');
 assert.throws(()=>exponentialPartition(g,8));negatives.push('unfrozen_partition_rejected');
 return{checks,negatives};
}
if(require.main===module){const r=verify();console.log(JSON.stringify({status:'PASS',checks:r.checks.length,negative_controls:r.negatives.length}));}
module.exports={verify};
