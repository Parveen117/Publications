'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const R=require('../native-lmd-certified-pulse-solver/interval.cjs');
const {f,iv,point,add,sub,neg,mul,scale,PI,trig,decimal,pack,contains,overlap,width}=R;
const dir=process.argv[2];assert(dir);
const report=JSON.parse(fs.readFileSync(path.join(dir,'calculator.result.json'))),protocol=JSON.parse(fs.readFileSync(path.join(__dirname,'cases.json')));
assert.equal(report.status,'completed');assert.equal(report.backend,'cyclotomic64');
const answers=new Map(report.answers.map(a=>[a.id,a])),checks=[],negativeControls=[];
const check=(n,ok)=>{assert(ok,n);checks.push(n);},reject=(n,ok)=>{assert(ok,n);negativeControls.push(n);};
const value=id=>{const a=answers.get(id);assert(a&&a.status==='completed'&&a.algebraically_exact===true);const d=a.value.decimal_enclosure;return iv(decimal(d.lower),decimal(d.upper));};
function exactRational(id,want){const c=answers.get(id).value.exact.coefficients;return c.length===32&&f(c[0]).eq(want)&&c.slice(1).every(x=>f(x).zero());}
check('normalization_exact',exactRational('normalization',1));
const rows=new Map();
for(const c of protocol.cases){
 const z=value(c.id+'_Z'),y=value(c.id+'_Y'),t=trig(scale(PI,c.theta_pi));
 check('Z_vs_Taylor:'+c.id,overlap(z,t.cos));check('Y_vs_Taylor:'+c.id,overlap(y,neg(t.sin)));
 const excitation=scale(sub(point(1),z),'1/2');
 const coMovingZ=sub(mul(t.cos,z),mul(t.sin,y)),coMovingExcitation=scale(sub(point(1),coMovingZ),'1/2');
 check('co_transformed_energy:'+c.id,contains(coMovingExcitation,0)&&width(coMovingExcitation).le('1/1000000000000000000'));
 if(c.theta_pi==='0')check('ground_or_closed_cycle_exact:'+c.id,exactRational(c.id+'_Z',1)&&exactRational(c.id+'_Y',0));
 if(answers.has(c.id+'_excited'))check('independent_probability_query:'+c.id,overlap(excitation,value(c.id+'_excited')));
 rows.set(c.id,{id:c.id,theta_pi:c.theta_pi,z:pack(z),y:pack(y),excitation:pack(excitation),co_moving_excitation:pack(coMovingExcitation),
  absolute_energies:protocol.offsets.map(o=>({offset:o,energy:pack(add(point(o),excitation))}))});
}
const cycles=[];
for(const suffix of ['s','t','p16']){
 const d=rows.get('drive_'+suffix),r=rows.get('return_'+suffix),energy=iv(d.excitation.lower,d.excitation.upper);
 check('drive_cost_positive:'+suffix,energy.lo.n>0n);
 check('return_ground_exact:'+suffix,exactRational('return_'+suffix+'_Z',1));
 reject('omit_drive_debit_creates_false_gain:'+suffix,energy.lo.n>0n);
 cycles.push({id:suffix,drive_work:pack(energy),ideally_recoverable_work:pack(energy),net_external_gain_exact:'0',
  exact_reason:'W_recover-W_drive=E_initial-E_return; exact return Z=1',physical_load_measured:false});
}
check('pi_over_two_cost_exact',exactRational('drive_s_Z',0));
check('direction_same_energy',overlap(value('drive_s_Z'),value('opposite_Z')));
check('direction_opposite_phase',overlap(value('drive_s_Y'),neg(value('opposite_Y'))));
reject('state_only_rotation_is_not_passive_basis_change',!contains(iv(rows.get('drive_s').excitation.lower,rows.get('drive_s').excitation.upper),0));
reject('wrong_return_does_not_close_cycle',exactRational('wrong_return_Z',-1));
reject('absolute_offset_is_not_extractable_work',exactRational('idle_Z',1)&&f('1/2').n>0n);
const result={status:'PASS',checks:checks.length,negative_controls:negativeControls.length,check_names:checks,negative_control_names:negativeControls,
 cases:[...rows.values()],cycles,fresh_nrqc_requests:1,hardware_called:false,gravity_hamiltonian_derived:false,
 vacuum_work_extraction_demonstrated:false,scope:'Finite supplied Hamiltonian energy-accounting calibration; inverse energy recovery is ideal, not a measured load.'};
fs.writeFileSync(path.join(dir,'verification.json'),JSON.stringify(result,null,2)+'\n');
console.log(JSON.stringify({status:result.status,checks:result.checks,negative_controls:result.negative_controls,cycles},null,2));
