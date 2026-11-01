'use strict';
// Post-execution planning diagnostic, not another probability simulation.
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict'),crypto=require('node:crypto');
const R=require('./source.cjs'),{f,iv,point,add,sub,mul,scale,exp,sinh,PI,absmax,pack}=R;
const filename=path.join(__dirname,'../native-lmd-pulse-characteristic/results/source_intervals.json');
const bytes=fs.readFileSync(filename),source=JSON.parse(bytes),h=f('3/8'),t=f('1/2'),k=f(2),target=f('1/200');
const damping=exp(point(h.mul(k.pow(2)).div(4).neg()));
const a=scale(sinh(point(h.mul(k).div(2))),t.mul(2).div(h));
const transfer=scale(mul(damping,mul(a,a)),'1/4');
const step=scale(PI,'1/32'),phaseError=absmax(mul(mul(step,step),mul(step,step))).div(24);
const phaseBudget=scale(damping,phaseError.div(2));
const rows=source.cutoff_moments.find(c=>c.h==='3/8').rows.map(r=>{
 const residual=iv(r.mass_L2_error_squared.lower,r.mass_L2_error_squared.upper);
 const budget=mul(transfer,residual),margin=sub(point(target),add(budget,phaseBudget));
 return{N:r.N,ideal_mode_probability_budget:pack(budget),
  remaining_budget_after_phase_interpolation:pack(margin),
  leaves_room_for_finite_source_compression:margin.lo.n>0n};
});
assert(rows[0].ideal_mode_probability_budget&&rows.length===6);
// If all Fourier orientations are discarded independently, conditioning on
// every radius captures at most exp(h + h^2*pi^2/24) as its second moment.
// I_0(x)<=exp(x^2/4), directly from the positive power-series coefficients.
const amplitudeSecondUpper=exp(add(point(h),scale(mul(PI,PI),h.pow(2).div(24))));
const fullMoment=iv(source.integrals['3/8'].lower,source.integrals['3/8'].upper);
const amplitudeVarianceLowerExpression=sub(fullMoment,amplitudeSecondUpper);
const amplitudeBudgetLowerExpression=mul(transfer,amplitudeVarianceLowerExpression);
assert(target.le(amplitudeBudgetLowerExpression.lo),'amplitude-only budget obstruction not certified');
const output={purpose:'Post-run source-resolution planning; no new characteristic expectation evaluated',
 h:h.toString(),t:t.toString(),k:k.toString(),target_halfwidth:target.toString(),
 phase_interpolation_probability_budget:pack(phaseBudget),rows,
 first_tested_cutoff_leaving_room:rows.find(r=>r.leaves_room_for_finite_source_compression)?.N??null,
 all_amplitudes_only:{second_moment_upper_bound_expression:pack(amplitudeSecondUpper),
   residual_variance_lower_bound_expression:pack(amplitudeVarianceLowerExpression),
   probability_certificate_lower_bound_expression:pack(amplitudeBudgetLowerExpression),
   exceeds_requested_halfwidth:true,
   meaning:'A bound on this quadratic conditional-variance certificate even with all mode radii; not a lower bound on actual probability approximation error or a no-go for other certificates.'},
 caveat:'This compares previously certified cutoffs 0,1,2,4,8,16 only. It does not identify the smallest possible integer cutoff or prove an executable finite partition meets tolerance.',
 source_sha256:crypto.createHash('sha256').update(bytes).digest('hex'),fresh_nrqc_requests:0};
const text=JSON.stringify(output,null,2)+'\n';
if(process.argv[2])fs.writeFileSync(process.argv[2],text);else process.stdout.write(text);
