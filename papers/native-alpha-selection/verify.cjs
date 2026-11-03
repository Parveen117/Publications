'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const here=__dirname,repo=path.resolve(here,'../..'),args=process.argv.slice(2),pos=args.indexOf('--rkf-root');
if(pos<0||!args[pos+1]||(args.includes('--write')&&args.includes('--check')))throw Error('Supply --rkf-root PATH and --check or --write');
const rkf=path.resolve(args[pos+1]),sha=b=>crypto.createHash('sha256').update(b).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(here,'SOURCE_PINS.json'),'utf8'));
for(const [root,files] of [[repo,pins.publications_files],[rkf,pins.rkf_files]])for(const f of files)
 assert.equal(sha(fs.readFileSync(path.join(root,f.path))),f.sha256,'Source drift: '+f.path);
const o=require(path.join(rkf,'operator_foundation/core/native_operator.cjs'));
const p=require(path.join(rkf,'operator_foundation/core/paninian_operator.cjs'));
const v=require(path.join(rkf,'research/recognition_return/r2/variable_return.cjs'));
const r1=require(path.join(rkf,'research/recognition_return/return_solver.cjs'));
const {WeightedCompletion}=require(path.join(rkf,'operator_foundation/core/weighted_completion.cjs'));
const {F,Cut}=o,f=x=>new F(x),eq=(a,b)=>assert(F.of(a).eq(b));
const counts={native_pair_proof_replays:0,cut_families:0,jet_coefficients:0,slope_enclosures:0,
 full_source_inversions:0,basis_transports:0,scale_controls:0};
const results={};
function prove(s,a,b){const e=p.proveEquality(a,b);assert.equal(e.status,'EQUAL_IN_DECLARED_QUOTIENT');assert(s.replay(e.certificate));return e.certificate;}
const s=r1.emk(),I=s.one(),R=s.word(['R']),K=s.word(['K']),L=K.times(R),Q=I.plus(L).scale(new Cut('1/2'));
results.cut_proof=prove(s,Q.times(Q),Q);
const families=[];
for(const bs of ['1/4','1/2','1','2','3','5','10']){
 const b=f(bs),a=b.add(1),slope=f(1).div(b.mul(2).add(1)),second=b.mul(2).mul(b.mul(4).add(3)).neg().div(b.mul(2).add(1).pow(3));
 const design=v.synthesize(1,b);eq(design.period[0],a);eq(design.normalizedCutDefect,0);
 eq(f(1).sub(design.otherBoundaryReturn).div(f(1).add(design.otherBoundaryReturn)),slope);
 const w=v.pairWitness(a,a,{u:1,v:design.otherBoundaryReturn});
 for(const proof of Object.values(w.proofs)){assert(s.replay(proof));counts.native_pair_proof_replays++;}
 const c=[f(1)];
 for(let n=1;n<=6;n++){
  let sum=f(0);for(let i=1;i<n;i++)sum=sum.add(c[i].mul(c[n-i]));
  for(let i=0;i<n;i++)sum=sum.add(c[i].mul(c[n-1-i]));
  c.push((n===1?a:f(0)).sub(b.mul(sum)).div(b.mul(2).add(1)));
 }
 // Independent polynomial convolution of the full implicit equation through degree six.
 const sq=Array.from({length:7},()=>f(0));
 for(let i=0;i<7;i++)for(let j=0;j<7-i;j++)sq[i+j]=sq[i+j].add(c[i].mul(c[j]));
 for(let n=0;n<7;n++){
  const residual=b.mul(sq[n].add(n?sq[n-1]:0)).add(c[n]).sub(n<2?a:0);
  eq(residual,0);counts.jet_coefficients++;
 }
 eq(c[1],slope);eq(c[2].mul(2),second);
 const slopeBounds=[];
 for(const ds of ['1/16','1/32']){
  const delta=f(ds),zminus=f(1).sub(delta),zplus=f(1).add(delta);
  const minus=v.solvePeriodic({pattern:[a.mul(zminus),b.mul(zminus)],tolerance:'1/10000000000',maxCells:1024});
  const plus=v.solvePeriodic({pattern:[a.mul(zplus),b.mul(zplus)],tolerance:'1/10000000000',maxCells:1024});
  assert.equal(minus.status,'CERTIFIED_RECOGNITION_RESPONSE');assert.equal(plus.status,'CERTIFIED_RECOGNITION_RESPONSE');
  const lower=plus.interval.lower.sub(1).div(delta),upper=f(1).sub(minus.interval.lower).div(delta);
  assert(lower.le(slope)&&slope.le(upper));
  slopeBounds.push({delta,lower,upper,minusCells:minus.cellsUsed,plusCells:plus.cellsUsed});counts.slope_enclosures++;
 }
 families.push({b,a,otherPhase:design.otherBoundaryReturn,slope,second,jets:c,slopeBounds});counts.cut_families++;
}
assert(!families[2].slope.eq(families[3].slope));
assert(families[3].slopeBounds[1].upper.le(families[2].slopeBounds[1].lower)); // aperture-certified distinction
results.native_return_families=families;
// Genuine noncommutative source dressing through the unchanged native WC5 API.
const W=new WeightedCompletion(s.spec(),{R:1,K:1}),ws=W.presentation,wi=W.one(),wr=ws.word(['R']),wk=ws.word(['K']);
const data={a:wi.scale(5),b:wr,c:wr.scale(-1),d:wi.scale(2),f:wi.scale(2),g:wk,
 hiddenApprox:wi.scale(new Cut('1/2')),effectiveApprox:wi.scale(new Cut('2/9'))};
const dressed=W.schurResponse(data);assert.equal(dressed.status,'NATIVE_SCHUR_RESPONSE_CERTIFIED');
eq(dressed.visible_error,0);eq(dressed.hidden_error,0);
prove(ws,W.add(W.mul(data.a,dressed.visible),W.mul(data.b,dressed.hidden)),data.f);
prove(ws,W.add(W.mul(data.c,dressed.visible),W.mul(data.d,dressed.hidden)),data.g);
const erased=W.schurResponse({...data,g:W.zero()});
assert(!W.mass(W.sub(dressed.visible,erased.visible)).zero());
results.native_source_dressing={visible:dressed.visible.toJSON(),hidden:dressed.hidden.toJSON(),
 effective_source:dressed.effective_source,error:dressed.visible_error,erased_source_detected:true};
// Exact positive quadratic adapter: source residue and contact term against full inversion.
const scalar=a=>{assert(a.turn.zero());return a.rad;};
const matEq=(a,b)=>assert(o.equal(a,b));
for(const zeta of ['1/2','1','2','3'])for(const B of ['-1','0','1'])for(const C of ['1','2','3'])for(const kk of ['1/2','1','2']){
 const z=f(zeta),b=f(B),c=f(C),k=f(kk),qv=f(2),qh=f(1),Z=z.mul(k.pow(2));
 const H=o.matrix([[Z.add(b.pow(2).div(c)),b],[b,c]]),source=o.matrix([[qv],[qh]]);
 const inv=o.inverse(H);matEq(o.mul(H,inv),o.identity(2));
 const full=scalar(o.mul(o.dagger(source),o.mul(inv,source))[0][0]);
 const qe=qv.sub(b.mul(qh).div(c)),contact=qh.pow(2).div(c),residue=qe.pow(2).div(z);
 eq(full,residue.div(k.pow(2)).add(contact));
 const x=scalar(o.mul(inv,source)[0][0]);eq(x,qe.div(Z));
 // Hidden and visible changes of coordinates transport both source and cost.
 const T=o.matrix([['1/2',0],[0,3]]),Ht=o.mul(o.dagger(T),o.mul(H,T)),qt=o.mul(o.dagger(T),source);
 eq(scalar(o.mul(o.dagger(qt),o.mul(o.inverse(Ht),qt))[0][0]),full);counts.basis_transports++;
 eq(qe.div(2).pow(2).div(z.div(4)),residue);
 if(!qe.zero())assert(!qe.pow(2).div(z.mul(2)).eq(residue));
 counts.full_source_inversions++;
}
results.residue_example={q_effective:'3/2',contact:'1/2',residue_zeta_1:'9/4',residue_zeta_2:'9/8',hidden_source_erased_residue_zeta_1:'4'};
// CT path and response-connection scale invariance; energy scale itself is not invariant.
const H=o.matrix([[2,1],[1,5]]),dH=o.matrix([[0,1],[1,0]]),gradient=o.matrix([[2],[3]]);
const X=o.scale(o.mul(o.inverse(H),gradient),-1),A=o.scale(o.mul(o.inverse(H),dH),'1/2');
for(const lambda of ['1/3','1','2','5']){
 const Hl=o.scale(H,lambda);matEq(o.scale(o.mul(o.inverse(Hl),o.scale(gradient,lambda)),-1),X);
 matEq(o.scale(o.mul(o.inverse(Hl),o.scale(dH,lambda)),'1/2'),A);
 counts.scale_controls++;
}
const sources=fs.readdirSync(here).filter(x=>/\.(md|cjs|json)$/.test(x)&&x!=='CERTIFICATE.json').sort();
const report=p.stable(JSON.parse(JSON.stringify({protocol:'NATIVE_ALPHA_SELECTION_R1',status:'PASS_EXACT_NATIVE_SENSITIVITY_AND_DRESSED_SOURCE_CONTROLS',
 alpha_status:'NOT_DERIVED',experimental_alpha_input:false,counts,results,rkf_commit:pins.rkf_commit,
 source_sha256:Object.fromEntries(sources.map(x=>[x,sha(fs.readFileSync(path.join(here,x)))])),
 scope:{selected_native_vacuum:false,physical_electron_source:false,photon_matching:false,
 unique_action_normalization:false,upstream_master_recertified:false,proof_assistant:false}})));
const text=JSON.stringify(report,null,2)+'\n',digest=sha(text);
if(args.includes('--write')){fs.writeFileSync(path.join(here,'CERTIFICATE.json'),text);fs.writeFileSync(path.join(here,'EXPECTED.sha256'),digest+'\n');}
else{assert.equal(fs.readFileSync(path.join(here,'CERTIFICATE.json'),'utf8'),text,'Certificate mismatch; no rewrite');assert.equal(fs.readFileSync(path.join(here,'EXPECTED.sha256'),'utf8').trim(),digest);}
console.log(JSON.stringify({status:report.status,alpha_status:report.alpha_status,sha256:digest,counts,runtime:process.version,mode:args.includes('--write')?'write':'check'}));
