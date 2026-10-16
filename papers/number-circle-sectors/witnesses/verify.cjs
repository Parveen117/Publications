'use strict';
// NCS-1 finite exact witnesses for NUMBER_CIRCLE.md. Regression evidence only:
// no continuum statement (Theorem C, Proposition C3, Conjecture C4) is
// certified by any finite check. Arithmetic: the unchanged native engine over
// Q + iota Q. Finite groups below are models, not Diff(S^1).
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../../native-common-cone-limit/native_engine/native_operator.cjs'),f=N.F.of,C=N.Cut.of;
const checks=[],negative_controls=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const col=v=>N.matrix(v.map(x=>[x]));
const dot=(u,v)=>N.mul(N.dagger(u),v)[0][0];
const pairs=[['3','5'],['-7/2','1/3'],['0','4'],['5/6','5/6'],['-2','-9'],['11/7','-3/4']];

// ---- W1: the common/disagreement transform ---------------------------------
// S = sqrt(2) R(-pi/4) acts on (o,x): S(o,x) = (o+x, x-o) = sqrt(2)(c,d). All checks are
// exact over Q because the sqrt(2) is carried as the scale factor 2 = (sqrt 2)^2.
const S=N.matrix([[1,1],[-1,1]]), I2=N.identity(2), A=N.matrix([[0,-1],[1,0]]);
ck('W1_S_transpose_S_is_2I',N.equal(N.mul(N.dagger(S),S),N.scale(I2,2)));
ck('W1_det_S_positive_rotation_not_reflection',N.equal(N.mul(S,N.matrix([[1,-1],[1,1]])),N.scale(I2,2)));
ck('W1_S_squared_is_2_times_quarter_turn',N.equal(N.mul(S,S),N.scale(N.matrix([[0,1],[-1,0]]),2)));
ck('W1_S_fourth_is_minus_4I_eighth_is_16I',N.equal(N.power(S,4),N.scale(I2,-4))&&N.equal(N.power(S,8),N.scale(I2,16)));
ck('W1_S_is_sqrt2_times_I_minus_A_over_sqrt2',N.equal(S,N.sub(I2,A)));
ck('W1_invariant_c2_plus_d2_equals_o2_plus_x2',pairs.every(([o,x])=>{const v=col([o,x]),w=N.mul(S,v);return dot(w,w).eq(dot(v,v).mul(C(2)));}));
ck('W1_inverse_recovers_o_and_x',pairs.every(([o,x])=>{const v=col([o,x]);return N.equal(N.scale(N.mul(N.dagger(S),N.mul(S,v)),C('1/2')),v);}));
ck('W1_d_zero_iff_o_equals_x',pairs.every(([o,x])=>N.mul(S,col([o,x]))[1][0].zero()===(f(o).eq(f(x)))));
// complex form on the number circle: c + iota d = ((1 - iota)/sqrt2)(o + iota x)
ck('W1_complex_form_multiplication_by_unit_circle_point',pairs.every(([o,x])=>{const q=new N.Cut(o,x),w=N.mul(S,col([o,x]));return C([1,-1]).mul(q).eq(new N.Cut(w[0][0].rad,w[1][0].rad));})&&C([1,-1]).norm2().eq(2));
// common shift (o+b,x+b) moves only c; opposite shift (o-b,x+b) moves only d
ck('W1_common_shift_fixes_d_opposite_shift_fixes_c',pairs.every(([o,x])=>['1/2','-3','7'].every(b=>{
 const w0=N.mul(S,col([o,x])),wc=N.mul(S,col([f(o).add(b),f(x).add(b)])),wo=N.mul(S,col([f(o).sub(b),f(x).add(b)]));
 return wc[1][0].eq(w0[1][0])&&!wc[0][0].eq(w0[0][0])&&wo[0][0].eq(w0[0][0])&&!wo[1][0].eq(w0[1][0]);})));
no('W1_disagreement_alone_recovers_inputs',!(()=>{const d1=N.mul(S,col(['1','3']))[1][0],d2=N.mul(S,col(['5','7']))[1][0];return !d1.eq(d2);})());
no('W1_unscaled_S_preserves_the_quadratic_reading',!pairs.every(([o,x])=>{const v=col([o,x]),w=N.mul(S,v);return dot(w,w).eq(dot(v,v));}));

// ---- W2: frame-labelled operators ------------------------------------------
// U(t) = cayley(A,t) = exp(Phi(t) A), rational points of the circle (NAC-1 (2)).
const U=t=>N.cayley(A,C(t)), Ui=t=>N.inverse(U(t));
const lab=(M,b,a)=>N.mul(N.mul(Ui(b),M),U(a)); // M_{b<-a}
const M=N.matrix([['2','-1/3'],['5','7/2']]), Nn=N.matrix([['-1','4'],['1/2','3']]);
const ts=['0','1/2','-1/3','2','5/7'];
ck('W2_cayley_points_are_rotations',ts.every(t=>N.equal(N.mul(N.dagger(U(t)),U(t)),I2)&&N.equal(N.mul(U(t),A),N.mul(A,U(t)))));
ck('W2_frame_composition_middle_frame_cancels',ts.every(a=>ts.every(b=>['3','-1/4'].every(g=>N.equal(N.mul(lab(Nn,g,b),lab(M,b,a)),N.mul(N.mul(Ui(g),N.mul(Nn,M)),U(a)))))));
no('W2_mismatched_middle_frames_compose',!ts.every(a=>N.equal(N.mul(lab(Nn,'3','1/2'),lab(M,'2',a)),N.mul(N.mul(Ui('3'),N.mul(Nn,M)),U(a)))));
// tangent-addition: U(t1)U(t2) = U((t1+t2)/(1-t1 t2)); frame dependence only through
// the angle difference holds iff M commutes with A (M = a I + b A, the number-circle case).
const tadd=(t1,t2)=>f(t1).add(f(t2)).div(f(1).sub(f(t1).mul(f(t2))));
ck('W2_cayley_tangent_addition',[['1/2','1/3'],['-1/3','2'],['5/7','-1/4']].every(([p,q])=>N.equal(N.mul(U(p),U(q)),U(tadd(p,q)))));
const Mq=N.add(N.scale(I2,'3/2'),N.scale(A,'-4')); // q-type operator r(cos th + A sin th)
const shiftPairs=[['1/2','1/3'],['2','-1/3'],['1/4','5/7']]; // (a,b) vs (tadd(a,s),tadd(b,s)) share the difference
ck('W2_q_type_operator_depends_only_on_angle_difference',shiftPairs.every(([a,b])=>['1/5','-2'].every(s=>N.equal(lab(Mq,b,a),lab(Mq,tadd(b,s),tadd(a,s))))));
ck('W2_q_type_operator_commutes_with_generator',N.isZero(N.commutator(Mq,A)));
no('W2_generic_operator_depends_only_on_angle_difference',!shiftPairs.every(([a,b])=>N.equal(lab(M,b,a),lab(M,tadd(b,'1/5'),tadd(a,'1/5')))));
// two-angle content: for M = diag(3,1) rotated, some (alpha,beta) make M_{beta<-alpha} diagonal
const Md=N.mul(N.mul(U('1/2'),N.matrix([[3,0],[0,1]])),Ui('1/3'));
const diagLab=lab(Md,'1/2','1/3');
ck('W2_two_frame_angles_diagonalize_a_non_q_operator',diagLab[0][1].zero()&&diagLab[1][0].zero()&&diagLab[0][0].eq(C(3))&&diagLab[1][1].eq(C(1)));
no('W2_one_common_angle_diagonalizes_it',!ts.some(t=>{const L=lab(Md,t,t);return L[0][1].zero()&&L[1][0].zero();}));

// ---- W3: sectored block operators (mixed-block gate, MZC-2 section 3) -------
const B=N.matrix([[1,2],[2,1]]),v=col(['1','-1']);
ck('W3_sector_compressions_positive',B[0][0].rad.n>0n&&B[1][1].rad.n>0n);
no('W3_positive_compressions_imply_positive_operator',!(dot(v,N.mul(B,v)).rad.n>=0n));
// conjugating the two-sector operator by S turns it block-diagonal: its c/d channels reduce it
const BS=N.mul(N.mul(S,B),N.dagger(S));
ck('W3_cd_frame_reduces_the_mixed_block',BS[0][1].zero()&&BS[1][0].zero()&&BS[0][0].eq(C(6))&&BS[1][1].eq(C(-2)));
// four-sector convention M_ij : sector j -> sector i; M_12 reads sector 2 into sector 1
const E=(i,j)=>{const Z=N.zeros(4,4);Z[i][j]=N.ONE;return Z;};
const e=k=>col([0,1,2,3].map(i=>i===k?1:0));
ck('W3_block_convention_M12_maps_sector2_to_sector1',N.equal(N.mul(E(0,1),e(1)),e(0))&&N.isZero(N.mul(E(0,1),e(0))));
no('W3_block_order_is_irrelevant',!N.equal(N.mul(E(0,1),E(1,2)),N.mul(E(1,2),E(0,1))));

// ---- W4: the chiral constraint pair is a c/d transform ----------------------
// IGQ-T4: T_eps = -eps Q_eps, H = T_+ + T_-, D = T_+ - T_-. Labels as coordinate vectors.
const K=N.matrix([[-1,0],[0,1]]);              // (Q_+,Q_-) -> (T_+,T_-)
const HD_T=N.matrix([[1,1],[1,-1]]);           // (T_+,T_-) -> (H,D)
const HD_Q=N.mul(HD_T,K);                       // (Q_+,Q_-) -> (H,D)
const swap=N.matrix([[0,1],[1,0]]);
// o = T_-, x = T_+ : S(o,x) = (o+x, x-o) = (H, D)
ck('W4_T_labels_H_D_equal_S_on_(T_minus,T_plus)',N.equal(N.mul(S,swap),HD_T));
// geometric labels: c-direction Q_+ + Q_- is -D, d-direction Q_+ - Q_- is -H (o = Q_-, x = Q_+)
const cdQ=N.mul(S,swap); // (Q_+,Q_-) -> (Q_- + Q_+, Q_+ - Q_-)
ck('W4_Q_labels_common_is_minus_D_disagreement_is_minus_H',N.equal(N.mul(N.matrix([[0,-1],[-1,0]]),HD_Q),cdQ));
ck('W4_diagonal_flow_direction_is_D',N.equal(N.mul(HD_Q,col([1,1])),col([0,-2])));
no('W4_chiral_sign_is_irrelevant_to_which_is_common',!N.equal(N.mul(HD_T,col([1,1])),N.mul(HD_Q,col([1,1]))));

// ---- W5: relational comparison in finite group models -----------------------
// X = G acting on itself by left multiplication on each factor: x_eps = chi_eps.
// Delta = {(g,g)}, H-type elements (h, h^{-1}), rho(x_+,x_-) = x_+^{-1} x_-.
const permMat=p=>N.matrix(p.map((_,r)=>p.map((_,c)=>p[c]===r?1:0)));
const key=Mx=>JSON.stringify(Mx.map(r=>r.map(x=>x.toJSON())));
const groupFrom=gens=>{const seen=new Map(gens.map(g=>[key(g),g]));let fr=[...seen.values()];
 while(fr.length){const nx=[];for(const a of fr)for(const b of gens){const m=N.mul(a,b);const k=key(m);if(!seen.has(k)){seen.set(k,m);nx.push(m);}}fr=nx;}
 return [...seen.values()];};
const model=(G)=>{
 const idx=new Map(G.map((g,i)=>[key(g),i])),inv=G.map(g=>idx.get(key(N.inverse(g)))),mul=(a,b)=>idx.get(key(N.mul(G[a],G[b])));
 const n=G.length,T=[];for(let a=0;a<n;a++){T.push([]);for(let b=0;b<n;b++)T[a].push(mul(a,b));}
 const rho=(p,m)=>T[inv[p]][m];
 const orbits=(acts)=>{const lab=new Int32Array(n*n).fill(-1);let k=0;
  for(let s=0;s<n*n;s++){if(lab[s]>=0)continue;const st=[s];lab[s]=k;while(st.length){const z=st.pop(),p=Math.floor(z/n),m=z%n;for(const [gp,gm] of acts){const w=T[gp][p]*n+T[gm][m];if(lab[w]<0){lab[w]=k;st.push(w);}}}k++;}
  return {count:k,lab};};
 const diag=[...Array(n).keys()].map(g=>[g,g]),htype=[...Array(n).keys()].map(h=>[h,inv[h]]);
 const full=[];for(let a=0;a<n;a++)full.push([a,0],[0,a]);
 const genSub=acts=>{const seen=new Set(acts.map(([a,b])=>a*n+b));let fr=[...seen];while(fr.length){const nx=[];for(const z of fr)for(const [a,b] of acts){const w=T[Math.floor(z/n)][a]*n+T[z%n][b];if(!seen.has(w)){seen.add(w);nx.push(w);}}fr=nx;}return seen.size;};
 return {n,T,inv,rho,orbits,diag,htype,full,genSub};};
// F21 = {y -> a y + b mod 7, a in {1,2,4}}: odd order, every element a square
const aff=(a,b)=>permMat([...Array(7).keys()].map(y=>(a*y+b)%7));
const F21=groupFrom([aff(1,1),aff(2,0)]);
const m21=model(F21);
ck('W5_F21_has_order_21',m21.n===21);
ck('W5_F21_every_element_is_a_square',(()=>{const sq=new Set();for(let a=0;a<21;a++)sq.add(m21.T[a][a]);return sq.size===21;})());
const oFull=m21.orbits(m21.full),oDiag=m21.orbits(m21.diag),oBoth=m21.orbits([...m21.diag,...m21.htype]);
ck('W5_full_product_action_single_orbit_invariant_dim_1',oFull.count===1);
ck('W5_diagonal_action_has_n_orbits_invariant_dim_21',oDiag.count===21);
ck('W5_rho_is_diagonal_invariant',[...Array(21).keys()].every(g=>{for(let p=0;p<21;p++)for(let m=0;m<21;m++)if(m21.rho(m21.T[g][p],m21.T[g][m])!==m21.rho(p,m))return false;return true;}));
ck('W5_rho_is_a_complete_diagonal_invariant',(()=>{const map=new Map();for(let p=0;p<21;p++)for(let m=0;m<21;m++){const o=oDiag.lab[p*21+m],r=m21.rho(p,m);if(map.has(o)&&map.get(o)!==r)return false;map.set(o,r);}return new Set(map.values()).size===21;})());
ck('W5_H_type_acts_on_rho_by_h_minus_2_conjugated',[...Array(21).keys()].every(h=>{const hi=m21.inv[h];for(let p=0;p<21;p++)for(let m=0;m<21;m++){
 const r=m21.rho(m21.T[h][p],m21.T[hi][m]),want=m21.T[m21.T[m21.inv[p]][m21.T[hi][hi]]][m];if(r!==want)return false;}return true;}));
ck('W5_diagonal_plus_H_type_generate_full_product',m21.genSub([...m21.diag,...m21.htype])===441);
ck('W5_adding_H_type_collapses_comparison_to_one_channel',oBoth.count===1);
no('W5_H_type_preserves_rho',![...Array(21).keys()].every(h=>{for(let p=0;p<21;p++)for(let m=0;m<21;m++)if(m21.rho(m21.T[h][p],m21.T[m21.inv[h]][m])!==m21.rho(p,m))return false;return true;}));
// S3: squares generate only A3, so Delta + H-type is index 2 and sign(rho) survives
const S3=groupFrom([permMat([1,0,2]),permMat([1,2,0])]),m6=model(S3);
const sgn=g=>(()=>{const P=S3[g];let s=0;for(let i=0;i<3;i++)for(let j=i+1;j<3;j++){const pi=P.findIndex(r=>r[i].rad.n===1n),pj=P.findIndex(r=>r[j].rad.n===1n);if(pi>pj)s^=1;}return s;})();
ck('W5_S3_control_generated_subgroup_index_2',m6.genSub([...m6.diag,...m6.htype])===18);
const o6=m6.orbits([...m6.diag,...m6.htype]);
ck('W5_S3_control_two_orbits_labelled_by_sign_of_rho',o6.count===2&&(()=>{for(let p=0;p<6;p++)for(let m=0;m<6;m++)if(o6.lab[p*6+m]!==o6.lab[0]&&sgn(m6.rho(p,m))===sgn(m6.rho(0,0)))return false;return true;})());
no('W5_H_type_generation_holds_without_square_generation',!(m6.genSub([...m6.diag,...m6.htype])===36));
// abelian toy: Z_7 rotations; Delta-invariant is the disagreement x - o
const Z7=groupFrom([aff(1,1)]),m7=model(Z7);
ck('W5_abelian_toy_diagonal_invariant_is_difference',m7.orbits(m7.diag).count===7);
// positive forms: rank of the Gram matrix of orbit indicators = #orbits (counting measure)
const gramRank=o=>{const n=o.count,len=o.lab.length,rows=[];for(let k=0;k<n;k++)rows.push([...Array(n).keys()].map(j=>j===k?[...o.lab].filter(x=>x===k).length:0));return N.rank(N.matrix(rows));};
ck('W5_invariant_positive_form_rank_21_diagonal_vs_1_full',gramRank(oDiag)===21&&gramRank(oFull)===1);

// ---- certificate -------------------------------------------------------------
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'NCS1_NUMBER_CIRCLE_SECTORS_V1',status:'PASS_FINITE_WITNESSES_FOR_NUMBER_CIRCLE_SECTORS',
 source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,
 results:{continuum_proofs_machine_certified:false,
  comparison_theorem:'written proof in NUMBER_CIRCLE.md section 5: momentum-constraint (diagonal) invariant dual is infinite-dimensional; full G_c dual is C ell',
  diagonal_classification:'Conjecture C4 (Delta-invariants are exactly functions of rho) is not decided'},
 summary:{exact_checks:checks.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
