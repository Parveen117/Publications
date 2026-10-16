'use strict';
// RDB-1 finite exact witnesses for BRIDGE.md. Regression evidence only: no
// continuum claim is certified; neither D's density nor its failure is decided
// by any finite check. Arithmetic: unchanged native engine.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const checks=[],negative_controls=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const Z=f(0),ONE=f(1);
// ---- W7: a conjugation-closed nontrivial set generates a simple group ------
const permMat=p=>N.matrix(p.map((_,r)=>p.map((_,c)=>p[c]===r?1:0))); // column c has 1 at row p[c]
const key=M=>JSON.stringify(M.map(r=>r.map(x=>x.toJSON())));
const closure=gens=>{const seen=new Map(gens.map(g=>[key(g),g]));let frontier=[...seen.values()];
 while(frontier.length){const next=[];
  for(const a of frontier)for(const b of seen.values())for(const m of [N.mul(a,b),N.mul(b,a)]){const k=key(m);if(!seen.has(k)){seen.set(k,m);next.push(m);}}
  frontier=next;}
 return seen.size;};
const cycles3=[];// all 3-cycles on 5 points: conjugation-closed in A5
for(let i=0;i<5;i++)for(let j=0;j<5;j++)for(let k=0;k<5;k++){if(i===j||j===k||i===k)continue;const p=[0,1,2,3,4];p[i]=j;p[j]=k;p[k]=i;cycles3.push(permMat(p));}
ck('conjugation_closed_class_generates_simple_A5',closure(cycles3)===60);
const dbl=[[1,0,3,2],[2,3,0,1],[3,2,1,0]].map(p=>permMat(p)); // double transpositions: a conjugacy class of A4
const v4=closure(dbl);
ck('nonsimple_A4_normal_closure_is_proper_V4',v4===4);
no('normal_closure_is_whole_group_without_simplicity',v4!==12);
// ---- W8: K_inf structure in a finite model ---------------------------------
// V=Q^4, core W=span(e1,e2,e3), ell=(1,0,0,0)^T dual, constraints map W into ker(ell|_W)=span(e2,e3).
const Q1=N.matrix([[0,0,0,0],[1,1,0,0],[0,0,0,0],[0,0,0,0]]); // Q1 e1=e2, Q1 e2=e2
const Q2=N.matrix([[0,0,0,0],[0,0,0,0],[1,0,1,0],[0,0,0,0]]); // Q2 e1=e3, Q2 e3=e3
const W=[[1,0,0,0],[0,1,0,0],[0,0,1,0]].map(v=>N.dagger(N.matrix([v]))); // column vectors
const images=[];for(const Q of [Q1,Q2])for(const w of W){const col=N.mul(Q,w);if(!col.every(r=>r[0].zero()))images.push(col.map(r=>r[0]));}
const rankImg=N.rank(N.matrix(images));
ck('solution_space_dimension_is_one_plus_codim',4-rankImg===1+(4-3));
// every solution restricted to W is proportional to ell|_W: solutions kill e2,e3 (images span them)
ck('constraint_images_span_kernel_of_ell_on_W',rankImg===2&&N.rank(N.matrix([...images,[0,1,0,0]]))===2&&N.rank(N.matrix([...images,[0,0,1,0]]))===2);
no('functional_with_F_e2_nonzero_solves_constraints',!(N.rank(N.matrix([...images,[0,1,0,0]]))===3));
// ---- W9: the two-point cylinder gap ----------------------------------------
// values psi(0,0)=0, psi(0,r)=1; any c depending on x1 only takes one value t at both points.
const grid=['-2','-1/2','0','1/4','1/2','3/4','1','2'].map(f);
ck('cylinder_two_point_error_at_least_half',grid.every(t=>{const a=t.abs(),b=ONE.sub(t).abs();const mx=a.le(b)?b:a;return f('1/2').le(mx);}));
ck('triangle_identity_abs_t_plus_abs_1_minus_t_at_least_1',grid.every(t=>ONE.le(t.abs().add(ONE.sub(t).abs()))));
{const c00=Z,c0r=ONE; // a function of both coordinates fits the two prescribed values exactly
 no('two_coordinate_function_also_blocked',!(!c00.sub(Z).zero()||!c0r.sub(ONE).zero()));}
// ---- W10: infinitesimal annihilation integrates to flow invariance ----------
const Nn=N.matrix([[0,1,0,0],[0,0,1,0],[0,0,0,1],[0,0,0,0]]);
const expN=s=>{let T=N.identity(4),P=N.identity(4),c=C(1);for(let k=1;k<4;k++){P=N.mul(P,Nn);c=c.mul(C(f(1).div(k)));const sc=C(f(s).pow(k)).mul(c);T=N.add(T,P.map(r=>r.map(x=>x.mul(sc))));}return T;};
const v=N.dagger(N.matrix([[1,1,1,1]]));
const Fv=M=>N.mul(N.matrix([[0,0,0,1]]),N.mul(M,v))[0][0]; // F = e4 dual; images N v=(1,1,1,0), N^2 v=(1,1,0,0), N^3 v=(1,0,0,0)
ck('F_annihilates_all_generator_images',[1,2,3].every(k=>{let M=N.identity(4);for(let i=0;i<k;i++)M=N.mul(M,Nn);return N.mul(N.matrix([[0,0,0,1]]),N.mul(M,v))[0][0].zero();}));
ck('flow_invariance_from_infinitesimal_annihilation',['-3','-1','1/2','2','7/3'].every(s=>Fv(expN(s)).eq(Fv(N.identity(4)))));
const F2=N.matrix([[1,-1,0,0]]); // kills N v and N^2 v but not N^3 v
no('partial_annihilation_gives_flow_invariance',!['1','2'].every(s=>N.mul(F2,N.mul(expN(s),v))[0][0].eq(N.mul(F2,v)[0][0])));
// ---- W11: slotwise constraint-annihilating forms factor ---------------------
// Model: gauge space Q^2 with core W=span(e1,e2)... take V2=Q^2, ell=(1,0), constraint image=span(e2).
// Forms on V2 (x) C^2: 4x4 M with M(Q w (x) s, .)=0 and (., Q w (x) s)=0 => rows/cols indexed by e2(x)C^2 vanish.
const E=(i,j)=>{const M=N.zeros(4,4);M[i][j]=N.ONE;return M;};
// basis of solutions after slot annihilation: M = e1e1^T (x) B, B 2x2 => 4 complex params at entries (0,0),(0,1),(1,0),(1,1) in the (gauge,sector) ordering index=2*g+s
const slotKill=M=>[2,3].every(i=>M[i].every(x=>x.zero())&&M.every(r=>r[i].zero())); // e2-slots are indices 2,3
const B4=(a,b,c,d)=>N.matrix([[a,b,0,0],[c,d,0,0],[0,0,0,0],[0,0,0,0]]);
ck('slot_annihilated_forms_live_on_ell_block',slotKill(B4(1,2,3,4)));
// sector condition: commute with sz on the sector factor: [M, I2 (x) sz]=0 kills off-diagonal b,c
const Isz=N.matrix([[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,-1]]);
const comm0=M=>N.sub(N.mul(M,Isz),N.mul(Isz,M)).every(r=>r.every(x=>x.zero()));
ck('sector_symmetry_forces_diagonal_weights',comm0(B4(1,0,0,4))&&!comm0(B4(1,2,3,4)));
// dimension count: {B4(a,0,0,d)} is 2-dimensional; adding sign mixing sx forces a=d (dimension 1)
const Isx=N.matrix([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]]);
const commx=M=>N.sub(N.mul(M,Isx),N.mul(Isx,M)).every(r=>r.every(x=>x.zero()));
ck('sign_mixing_collapses_weights_to_one',commx(B4(1,0,0,1))&&!commx(B4(1,0,0,4)));
// zero-weight degenerate members of the trichotomy at finite level
ck('trichotomy_ranks_2_1_0',N.rank(B4(1,0,0,4))===2&&N.rank(B4(1,0,0,0))===1&&N.rank(B4(0,0,0,0))===0);
no('zero_weight_member_has_full_block_rank',N.rank(B4(1,0,0,0))!==2);
// ---- certificate -------------------------------------------------------------
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'RDB1_CONSTRAINT_BRIDGE_V1',status:'PASS_FINITE_WITNESSES_FOR_CONSTRAINT_BRIDGE',
 source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,
 results:{continuum_proofs_machine_certified:false,
  bridge:'written proof in BRIDGE.md: K_inf = C ell (+) annihilator of closure(D); yes on the core',
  density_of_D_in_19:'not decided; equivalent open approximation problem recorded in BRIDGE.md 5(e)'},
 summary:{exact_checks:checks.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
