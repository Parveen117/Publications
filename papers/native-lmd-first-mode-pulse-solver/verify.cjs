'use strict';
const assert=require('node:assert/strict');
const R=require('./source.cjs');
const {f,point,add,sub,mul,exp,contains,overlap,bessel0,radialPartition,gaussianPartition}=R;
function verify(){
 const checks=[],negatives=[],check=(name,ok)=>{assert(ok,name);checks.push(name);};
 check('Bessel_zero_is_one',contains(bessel0(0),1));
 const g=f('1/8'),b=bessel0(g),partial=f(1).add(g.pow(2).div(4));
 check('Bessel_positive_correction',partial.le(b.lo));
 const coarse=radialPartition(g,2),fine=radialPartition(g,4);
 for(const [name,p] of [['coarse',coarse],['fine',fine]]){
  check(name+'_normalized',contains(p.totalP,1)&&contains(p.totalA,1));
  check(name+'_positive_cells',p.rows.every(r=>r.p.lo.n>0n&&r.first.lo.n>0n));
  check(name+'_mass_second_between_one_and_Bessel',f(1).le(p.second.lo)&&p.second.hi.le(b.lo));
  check(name+'_tail_retained',p.rows.at(-1).right===Infinity&&p.rows.at(-1).p.lo.n>0n);
 }
 check('radial_refinement_progress',coarse.second.hi.le(fine.second.lo));
 for(let j=0;j<coarse.rows.length-1;j++){
  check('nested_radial_probability:'+j,overlap(coarse.rows[j].p,add(fine.rows[2*j].p,fine.rows[2*j+1].p)));
  check('nested_radial_mass:'+j,overlap(coarse.rows[j].first,add(fine.rows[2*j].first,fine.rows[2*j+1].first)));
 }
 const q=gaussianPartition(g,4);
 check('q_probability_and_mass',contains(q.totalP,1)&&contains(q.totalA,1));
 check('q_quantization_second_moment',q.second.hi.le(exp(point(g)).lo));
 function reject(name,ok){assert(ok,name);negatives.push(name);}
 reject('dropping_radial_tail_loses_probability',!contains(sub(fine.totalP,fine.rows.at(-1).p),1));
 reject('omitting_Wick_factor_changes_mean',!contains(mul(fine.totalA,exp(point(g.div(2)))),1));
 reject('substituting_unit_radial_mass_loses_variance',!overlap(fine.second,point(1)));
 assert.throws(()=>radialPartition('-1/8',4));negatives.push('negative_source_variance_rejected');
 assert.throws(()=>radialPartition(g,64));negatives.push('over_budget_partition_rejected');
 return{checks,negatives};
}
if(require.main===module){const v=verify();console.log(JSON.stringify({status:'PASS',checks:v.checks.length,negative_controls:v.negatives.length}));}
module.exports={verify};
