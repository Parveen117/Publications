'use strict';
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const N=require('../native-common-cone-limit/native_engine/native_operator.cjs');
const f=N.F.of,C=N.Cut.of;
const trace=a=>a.reduce((s,row,i)=>s.add(row[i]),C(0));
const I=N.identity(2),X=N.matrix([[0,1],[1,0]]),Y=N.matrix([[0,[0,-1]],[[0,1],0]]),Z=N.matrix([[1,0],[0,-1]]);
const words=[N.kron(I,X),N.kron(Z,X),N.kron(I,Y),N.kron(Z,Y),N.kron(X,I),N.kron(X,Z)];
function controls(row){
 const [a,b,g]=[row.alpha,row.beta,row.gamma].map(f);
 const R=N.zeros(4),J=N.zeros(4);
 R[1][0]=R[0][1]=C(a);J[1][0]=C([0,b]);J[0][1]=C([0,b.neg()]);
 J[2][0]=J[0][2]=C(g);
 return [R,J];
}
function state(name){
 if(name==='vacuum'||name==='flipped'){const P=N.zeros(4),i=name==='vacuum'?0:1;P[i][i]=C(1);return P;}
 if(name==='bell')return N.scale(N.matrix([[1,0,0,1],[0,0,0,0],[0,0,0,0],[1,0,0,1]]),'1/2');
 const a=N.matrix([[1],[[0,1]],[1],[[0,1]]]);return N.scale(N.mul(a,N.dagger(a)),'1/4');
}
function covariance(P,H){
 const means=H.map(a=>trace(N.mul(P,a)));
 return H.map((a,i)=>H.map((b,j)=>trace(N.mul(P,N.mul(a,b))).sub(means[i].mul(means[j]))));
}
function weights(row){const a=f(row.alpha).div(2),b=f(row.beta).div(2),g=f(row.gamma).div(2);return N.matrix([[a,a,0,0,0,0],[0,0,b,b,g,g]]);}
function parseTensor(report){
 assert.equal(report.answers.length,1);const answer=report.answers[0];
 assert.equal(answer.algebraically_exact,true);
 return N.matrix(answer.geometry.tensor.map(row=>row.map(entry=>{
  const [a,b,c,d]=entry.coefficients;assert(f(b).zero()&&f(d).zero());return C([a,c]);
 })));
}
function median(xs){return [...xs].sort((a,b)=>a-b)[Math.floor(xs.length/2)];}
if(process.argv[2]==='--baseline'){
 const source=JSON.parse(fs.readFileSync(process.argv[3]));
 const start=process.hrtime.bigint();const result=source.cases.map(row=>covariance(state('vacuum'),controls(row)));
 process.stdout.write(JSON.stringify({result,compute_seconds:Number(process.hrtime.bigint()-start)/1e9})+'\n');
}else{
 const dir=process.argv[2],source=JSON.parse(fs.readFileSync(path.join(dir,'source_intervals.json')));
 const checks=[],negatives=[],timings=[];
 const check=(name,value)=>{assert(value,name);checks.push(name);};
 const reject=(name,value)=>{assert(value,name);negatives.push(name);};
 const reports={};
 for(const name of ['vacuum','flipped','bell','complex_product'])for(const backend of ['stabilizer','observer']){
  const report=JSON.parse(fs.readFileSync(path.join(dir,name+'_'+backend+'.result.json')));
  const G=parseTensor(report);reports[name+'_'+backend]=G;
  check(name+':'+backend+':independent_Pauli_tensor',N.equal(G,covariance(state(name),words)));
  for(const row of source.cases){
   const W=weights(row),start=process.hrtime.bigint();
   const Q=N.mul(N.mul(W,G),N.dagger(W));
   timings.push(Number(process.hrtime.bigint()-start)/1e9);
   check(name+':'+backend+':weighted_tensor:h='+row.h,N.equal(Q,covariance(state(name),controls(row))));
  }
 }
 for(const row of source.cases){
  const W=weights(row),Q=N.mul(N.mul(W,reports.vacuum_stabilizer),N.dagger(W));
  const A=f(row.A),lo=f(row.B.lower),hi=f(row.B.upper),c=f(row.c),tol=f('1/100000');
  check('vacuum:source_A_error:'+row.h,Q[0][0].rad.sub(A).abs().le(tol));
  check('vacuum:source_B_error_both_endpoints:'+row.h,Q[1][1].rad.sub(lo).abs().le(tol)&&Q[1][1].rad.sub(hi).abs().le(tol));
  check('vacuum:source_signed_cross_error:'+row.h,Q[0][1].turn.sub(c).abs().le(tol));
  check('vacuum:complex_rank_two:'+row.h,N.rank(Q)===2);
  const gram=covariance(state('vacuum'),controls({...row,gamma:'0'}));
  reject('drop_second_tangent_changes_finite_tensor:'+row.h,!N.equal(gram,Q));
  reject('one_pure_qubit_rank_one:'+row.h,N.rank(gram)===1&&N.rank(Q)===2);
  reject('leading_only_A_outside_tolerance:'+row.h,!A.sub(f(row.h).div(4)).abs().le(tol));
 }
 const corrupted=N.matrix(reports.vacuum_stabilizer);corrupted[0][2]=corrupted[0][2].neg();
 reject('wrong_XY_phase_is_detected',!N.equal(corrupted,covariance(state('vacuum'),words)));
 const row=source.cases[1],Q=covariance(state('vacuum'),controls(row));
 const wrong=controls(row);wrong[1]=N.scale(wrong[1],-1);
 reject('reversed_current_density_cross_sign',!N.equal(covariance(state('vacuum'),wrong),Q));
 const W=weights(row),bad=N.matrix(W);bad[0][1]=C(0);
 reject('omitted_projector_factor',!N.equal(N.mul(N.mul(bad,reports.vacuum_stabilizer),N.dagger(bad)),Q));
 const execution=JSON.parse(fs.readFileSync(path.join(dir,'execution.json')));
 const wall={};for(const engine of ['native','nrqc'])wall[engine]=median(execution.timings.filter(r=>r.engine===engine).map(r=>r.wall_seconds));
 const result={status:'PASS',scope:'three source covariance points, four finite states, two existing NRQC backends; moment encoding only',
   checks:checks.length,negative_controls:negatives.length,check_names:checks,negative_control_names:negatives,
   tolerance_per_entry:'1/100000',source_preparation_seconds:source.source_preparation_seconds,
   source_values:source.cases.map(r=>({h:r.h,...r.display})),
   median_cold_command_seconds:wall,median_weighted_assembly_seconds:median(timings),
   benchmark_limit:'native command computes three weighted tensors; NRQC command computes one shared six-Pauli tensor, reused for all three; NRQC weighted assembly is additional; shared source coefficients supplied to both',
   advantage_established:false,gravity_dynamics_simulated:false,
   inherited_continuum_claims:'conditional on LMD; finite executions do not prove them'};
 process.stdout.write(JSON.stringify(result,null,2)+'\n');
}
