'use strict';
// JED-1 finite exact witnesses for JET_ENTRY_DUALS.md. Regression evidence only:
// no continuum statement (Theorems 1-5 for arbitrary N, the flat-function
// obstruction, the analytic completeness statement) is certified by any finite
// check. Arithmetic: the unchanged native engine over Q + iota Q through the
// unchanged TLM-1 jet adapter and the TET-1 path adapter. Nothing is forked.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const N=require('../../native-common-cone-limit/native_engine/native_operator.cjs');
const L=require('../../native-tvsp-lambda-matrices/tvsp_adapter.cjs');
const T=require('../../native-tvsp-exponential-transport/transport_adapter.cjs');
const f=N.F.of,C=N.Cut.of,I=N.IOTA;
const checks=[],negative_controls=[];
const ck=(name,ok)=>{assert.ok(ok,name);checks.push({name,passed:true});};
const no=(name,ok)=>{assert.ok(ok,name);negative_controls.push({name,rejected:true});};
const realNonneg=z=>z.turn.zero()&&f(0).le(z.rad);
const fact=n=>{let r=1n;for(let k=2n;k<=BigInt(n);k++)r*=k;return r;};
const binom=(n,k)=>Number(fact(n)/(fact(k)*fact(n-k)));

// ---- helpers on top of the adapter ---------------------------------------------
const sc=(Cx,c)=>L.scalar(Cx,c);
const mscal=(A,x)=>L.matrix(A.rows.map(r=>r.map(e=>L.mul(e,x))));           // matrix jet times scalar jet
const constM=(Cx,M)=>L.matrix(N.matrix(M).map(r=>r.map(c=>sc(Cx,c))));      // constant matrix as a jet matrix
const trJ=M=>M.rows.reduce((s,r,i)=>L.add(s,r[i]),sc(M.context,0));         // J_N-valued trace
const trN=M=>M.reduce((s,r,i)=>s.add(r[i]),N.ZERO);                           // ordinary trace of a native matrix
const vecJ=x=>L.basisIndices(x.context).map(([a,b])=>L.coefficient(x,a,b));  // coefficient column of a scalar jet
const vecA=M=>M.rows.flatMap(r=>r.flatMap(vecJ));                             // coefficient column of a matrix jet
const dim=Cx=>L.basisIndices(Cx).length;
function expJet(Cx,X,t){ // exp(t X) for a constant matrix X and a scalar jet t with t_0 = 0 (finite sum)
  const m=X.length;let out=L.identity(Cx,m),term=L.identity(Cx,m),Xj=constM(Cx,X);
  for(let k=1;k<=Cx.order;k++){term=L.matrixScale(mscal(L.matrixMul(term,Xj),t),new N.F(1,k));out=L.matrixAdd(out,term);}
  return out;
}
function invJet(G){ // (G_0 + G_+)^{-1} = sum_k (-G_0^{-1} G_+)^k G_0^{-1}, exact in the nilpotent quotient
  const Cx=G.context,m=G.height,G0=L.matrixOuter(G),G0i=constM(Cx,N.inverse(G0));
  const Gp=L.matrixSub(G,constM(Cx,G0)),Q=L.matrixScale(L.matrixMul(G0i,Gp),-1);
  let out=L.identity(Cx,m),term=L.identity(Cx,m);
  for(let k=1;k<=Cx.order;k++){term=L.matrixMul(term,Q);out=L.matrixAdd(out,term);}
  return L.matrixMul(out,G0i);
}
const conj=(G,M)=>L.matrixMul(L.matrixMul(invJet(G),M),G);                  // G^{-1} * M * G
function subst(x,p,q){ // chart substitution: x(eps_s,eps_v) -> x(p,q), p and q scalar jets in the same context
  let out=sc(x.context,0);
  for(const [a,b,c] of L.entries(x))out=L.add(out,L.scale(L.mul(L.power(p,a),L.power(q,b)),c));
  return out;
}
const substM=(M,p,q)=>L.matrix(M.rows.map(r=>r.map(x=>subst(x,p,q))));
function spanRank(vectors){ // rank of the span of coefficient columns (exact rref over Q + iota Q)
  const rows=vectors.filter(v=>v.some(c=>!c.zero()));return rows.length?N.rank(N.matrix(rows)):0;
}
function sampleScalars(Cx){ // a few scalar jets with every layer populated
  const s=L.coordinate(Cx,'s'),v=L.coordinate(Cx,'v');
  const a=L.add(sc(Cx,2),L.add(s,L.scale(v,I)));
  const b=L.add(sc(Cx,'3/2'),L.sub(L.scale(s,-2),L.mul(v,v)));
  const c=L.add(sc(Cx,'-1/3'),L.add(L.mul(s,v),L.scale(L.mul(s,s),'1/5')));
  const d=L.add(sc(Cx,[1,1]),L.add(L.scale(s,[0,-1]),L.mul(L.mul(s,v),v)));
  return {s,v,a,b,c,d,one:sc(Cx,1),zero:sc(Cx,0)};
}
function sampleMatrices(Cx){
  const {s,v,a,b,c,d,one,zero}=sampleScalars(Cx);
  return {
    M:L.matrix([[a,b],[s,c]]),P:L.matrix([[v,d],[one,L.mul(a,b)]]),
    R:L.matrix([[a,s,v],[b,c,one]]),Sm:L.matrix([[d,zero],[s,c],[one,v]])
  };
}
const monomials=Cx=>L.basisIndices(Cx).map(([a,b])=>L.monomial(Cx,a,b));
const unitBasis=(Cx,m)=>{const out=[];for(let i=0;i<m;i++)for(let j=0;j<m;j++)for(const e of monomials(Cx)){
  const rows=Array.from({length:m},(_,r)=>Array.from({length:m},(_,cc)=>r===i&&cc===j?e:sc(Cx,0)));out.push(L.matrix(rows));}return out;};
// a functional given by layer weights {'a,b': Cut}
const evalW=(w,x)=>Object.entries(w).reduce((z,[k,c])=>{const [a,b]=k.split(',').map(Number);return z.add(C(c).mul(L.coefficient(x,a,b)));},N.ZERO);

const depths=[1,2,3];
for(const Nd of depths){
  const Cx=L.createContext(Nd,'JED witness chart'),d=dim(Cx),S=sampleScalars(Cx),Mx=sampleMatrices(Cx);
  const {s,v,a,b,c,one,zero}=S;
  ck('dimension:d_N:'+Nd,d===(Nd+1)*(Nd+2)/2);

  // ---- Theorem 1: characters and irreducible representations --------------------
  ck('T1_formal_coordinates_nilpotent_nonzero:'+Nd,L.equal(L.power(s,Nd+1),zero)&&L.equal(L.power(v,Nd+1),zero)&&!L.equal(s,zero));
  no('T1_a_nonzero_value_of_eps_defines_a_character:'+Nd,['1','-2','1/3'].every(t=>!C(t).mul(C(t)).zero()&&!(()=>{let p=C(1);for(let k=0;k<=Nd;k++)p=p.mul(C(t));return p.zero();})()));
  no('T1_a_layer_coefficient_is_multiplicative:'+Nd,(()=>{const chi=x=>L.outer(x).add(L.coefficient(x,1,0));return !chi(L.mul(a,b)).eq(chi(a).mul(chi(b)));})());
  { // regular module J_N: the positive-degree ideal acts nilpotently; its common kernel is the top degree;
    // every one-dimensional subrepresentation carries the observer character.
    const Rs=L.regularMatrix(L.matrix([[s]])),Rv=L.regularMatrix(L.matrix([[v]]));
    ck('T1_ideal_acts_nilpotently_with_index_N_plus_1:'+Nd,N.isZero(N.power(Rs,Nd+1))&&!N.isZero(N.power(Rs,Nd)));
    const stacked=N.matrix(Rs.concat(Rv)),kernelDim=d-N.rank(stacked);
    ck('T1_common_kernel_of_layers_is_top_degree:'+Nd,kernelDim===Nd+1);
    const basis=L.basisIndices(Cx),top=basis.map((ab,k)=>[ab,k]).filter(([[p,q]])=>p+q===Nd).map(([,k])=>k);
    const col=k=>N.matrix(basis.map((_,l)=>[l===k?1:0]));
    const Rx=L.regularMatrix(L.matrix([[a]])),Ry=L.regularMatrix(L.matrix([[c]]));
    ck('T1_top_degree_vectors_carry_observer_eigenvalue:'+Nd,top.every(k=>N.equal(N.mul(Rx,col(k)),N.scale(col(k),L.outer(a)))&&N.equal(N.mul(Ry,col(k)),N.scale(col(k),L.outer(c)))));
    no('T1_layers_act_trivially_on_the_regular_module:'+Nd,!N.isZero(Rs));
    no('T1_a_degree_zero_vector_is_a_common_eigenvector:'+Nd,!N.equal(N.mul(Rs,col(0)),N.zeros(d,1)));
  }

  // ---- Theorem 2: the positive dual collapses to the observer ------------------
  const pairs=[[a,b],[b,c],[a,S.d],[L.mul(a,s),c],[S.d,L.add(one,v)]];
  for(const rho of [[[2,1],[1,1]],[[1,0],[0,3]],[[1,[0,1]],[[0,-1],1]]]){ // PSD: det 1, trace>0; diag; [[1,i],[-i,1]] (det 0)
    const phi=M=>trN(N.mul(N.matrix(rho),L.matrixOuter(M)));
    const Ms=[Mx.M,Mx.P,L.matrixMul(Mx.M,Mx.P),L.matrixAdd(Mx.M,L.matrixScale(Mx.P,I))];
    ck('T2_point_state_is_positive_on_squares:'+Nd+':'+JSON.stringify(rho),Ms.every(M=>realNonneg(phi(L.matrixMul(L.matrixDagger(M),M)))));
    ck('T2_point_state_is_hermitian:'+Nd+':'+JSON.stringify(rho),Ms.every(A=>Ms.every(B=>phi(L.matrixMul(L.matrixDagger(B),A)).eq(phi(L.matrixMul(L.matrixDagger(A),B)).dagger()))));
  }
  { // every functional that reads a layer fails positivity on an explicit element 1 + t eps^alpha + r eps^(2 alpha)
    const ts=['1','-1','2','-2','1/2','-1/2','3','-3','100','-100','1000','-1000',[0,1],[0,-1]],rs=['0','-1','-4','-9','-16','1'];
    const violates=w=>L.basisIndices(Cx).some(([p,q])=>p+q>0&&ts.some(t=>rs.some(r=>{
      let x=L.add(one,L.scale(L.monomial(Cx,p,q),t));if(2*(p+q)<=Nd)x=L.add(x,L.scale(L.monomial(Cx,2*p,2*q),r));
      return !realNonneg(evalW(w,L.mul(L.dagger(x),x)));})));
    const readers=[];for(const [p,q] of L.basisIndices(Cx))if(p+q>0)for(const cw of ['1','-1/2',[0,1],'3'])readers.push({'0,0':'1',[p+','+q]:cw});
    readers.push({'0,0':'1','1,0':'1','0,1':'-2'},{'0,0':'2','1,0':'1','2,0':'1'},{'0,0':'1','1,1':'1','0,2':'1'});
    ck('T2_every_layer_reading_functional_violates_positivity:'+Nd,readers.filter(w=>Object.keys(w).every(k=>k.split(',').reduce((u,z)=>u+Number(z),0)<=Nd)).every(violates));
    ck('T2_the_observer_itself_is_positive_on_the_same_elements:'+Nd,!violates({'0,0':'1'}));
    // sampling cannot certify positivity: a reading functional nonnegative on every fixture square still fails
    const z=L.buildFixture('1',Nd),fixtureRecords=[z.a,z.b,z.c,z.Delta,z.T,z.P,z.r,...Object.values(z.responses)];
    const w={'0,0':'1','0,1':'1/100'};
    const batchOK=fixtureRecords.every(x=>realNonneg(evalW(w,L.mul(L.dagger(x),x))));
    no('T2_nonnegativity_on_fixture_squares_certifies_positivity:'+Nd,batchOK&&violates(w));
  }

  // ---- Theorem 3: frame-invariant linear functionals = layer traces ------------
  ck('T3_jet_trace_is_cyclic:'+Nd,L.equal(trJ(L.matrixMul(Mx.M,Mx.P)),trJ(L.matrixMul(Mx.P,Mx.M)))&&L.equal(trJ(L.matrixMul(Mx.R,Mx.Sm)),trJ(L.matrixMul(Mx.Sm,Mx.R))));
  const X=[[I,2],[-2,I.neg()]].map(r=>r.map(x=>C(x).mul(C('1/3')))),Y=[[0,'3/7'],['-3/7',0]];                  // skew-Hermitian
  const U0=[['3/5','-4/5'],['4/5','3/5']],V0=[[I,0],[0,1]];                                                        // constant unitaries
  const unitaryJets=[L.matrixMul(L.matrixMul(expJet(Cx,X,s),expJet(Cx,Y,v)),constM(Cx,U0)),L.matrixMul(constM(Cx,V0),expJet(Cx,Y,L.add(s,L.mul(s,v)))),expJet(Cx,X,L.sub(v,L.mul(s,s)))];
  ck('T3_sampled_jets_are_unitary:'+Nd,unitaryJets.every(G=>L.matrixEqual(L.matrixMul(L.matrixDagger(G),G),L.identity(Cx,2))));
  ck('T3_unitary_jet_dagger_is_inverse:'+Nd,unitaryJets.every(G=>L.matrixEqual(L.matrixDagger(G),invJet(G))));
  const G0=constM(Cx,[[2,1],[1,1]]),invertibleJets=[...unitaryJets,L.matrixMul(G0,L.matrixAdd(L.identity(Cx,2),L.matrixAdd(mscal(constM(Cx,[[1,'1/2'],[0,-1]]),s),mscal(constM(Cx,[[0,1],[1,0]]),v))))];
  ck('T3_sampled_jets_are_invertible:'+Nd,invertibleJets.every(G=>L.matrixEqual(L.matrixMul(invJet(G),G),L.identity(Cx,2))));
  ck('T3_jet_trace_invariant_under_all_sampled_frame_jets:'+Nd,invertibleJets.every(G=>[Mx.M,Mx.P].every(M=>L.equal(trJ(conj(G,M)),trJ(M)))));
  ck('T3_conjugation_actually_moves_the_record:'+Nd,invertibleJets.every(G=>!L.matrixEqual(conj(G,Mx.M),Mx.M)));
  no('T3_offdiagonal_layer_coefficient_is_frame_invariant:'+Nd,!unitaryJets.every(G=>L.coefficient(conj(G,Mx.M).rows[0][1],1,0).eq(L.coefficient(Mx.M.rows[0][1],1,0))));
  no('T3_jet_trace_invariant_under_independent_in_out_frames:'+Nd,!L.equal(trJ(L.matrixMul(L.matrixMul(L.matrixDagger(unitaryJets[0]),Mx.M),unitaryJets[1])),trJ(Mx.M)));
  { // necessity: constant signed permutations already force the trace on M_2(C); then the invariant dual of M_2(J_N) has dimension exactly d
    const signed=[[[1,0],[0,-1]],[[0,1],[1,0]]].map(U=>constM(Cx,U)),units=unitBasis(Cx,2);
    const diffs=Gs=>Gs.flatMap(G=>units.map(M=>vecA(L.matrixSub(conj(G,M),M))));
    const constDiffs=[[[1,0],[0,-1]],[[0,1],[1,0]]].flatMap(U=>[[[1,0],[0,0]],[[0,1],[0,0]],[[0,0],[1,0]],[[0,0],[0,1]]].map(E=>N.flatten(N.sub(N.mul(N.mul(N.inverse(N.matrix(U)),N.matrix(E)),N.matrix(U)),N.matrix(E)))));
    ck('T3_signed_permutations_force_the_trace_on_M2:'+Nd,spanRank(constDiffs)===3&&constDiffs.every(w=>w[0].add(w[3]).zero()));
    const r1=spanRank(diffs(signed)),r2=spanRank(diffs([...signed,...unitaryJets]));
    ck('T3_invariant_linear_dual_has_dimension_d:'+Nd,r1===3*d&&r2===3*d);
    const traceFunctionals=L.basisIndices(Cx).map(([p,q])=>M=>L.coefficient(trJ(M),p,q));
    ck('T3_layer_traces_vanish_on_all_differences:'+Nd,[...signed,...unitaryJets].every(G=>units.every(M=>traceFunctionals.every(t=>t(L.matrixSub(conj(G,M),M)).zero()))));
    no('T3_frame_jets_cut_the_invariant_dual_below_d:'+Nd,!(r2>r1));
  }
  { // the unique invariant positive functional is the trace state, which is layer blind; the jet trace reads the TLM deformation with its sign
    // each fixture carries its own chart context; comparisons across fixtures go through coefficient columns
    const readings=['-5','0','1/2','1','5'].map(mu=>{const z=L.buildFixture(mu,Nd);return {mu,tr:trJ(z.H),H:z.H,zc:z.context};});
    const sameCoeffs=(x,y)=>{const p=vecJ(x),q=vecJ(y);return p.length===q.length&&p.every((c,k)=>c.eq(q[k]));};
    ck('T3_jet_trace_of_TLM_source_is_7_plus_mu_eps_v:'+Nd,readings.every(({mu,tr,zc})=>L.equal(tr,L.add(sc(zc,7),L.scale(L.coordinate(zc,'v'),mu)))));
    ck('T3_trace_state_is_layer_blind_on_TLM_source:'+Nd,readings.every(({tr})=>L.outer(tr).eq(7)));
    ck('T3_jet_trace_separates_mu_from_minus_mu:'+Nd,!sameCoeffs(readings[0].tr,readings[4].tr)&&readings[0].H.height===2);
    no('T3_ordinary_trace_of_regular_encoding_reads_layers:'+Nd,readings.every(({H})=>trN(L.regularMatrix(H)).eq(C(7).mul(C(d)))));
    const Rs1=L.regularMatrix(L.matrix([[s]])),Rs2=L.regularMatrix(L.matrix([[L.scale(s,2)]]));
    ck('T3_jordan_type_detects_a_nonzero_layer:'+Nd,N.rank(Rs1)>0&&N.rank(L.regularMatrix(L.matrix([[zero]])))===0);
    no('T3_jordan_type_reads_the_layer_value:'+Nd,[1,2,3].every(k=>N.rank(N.power(Rs1,k))===N.rank(N.power(Rs2,k)))&&!L.equal(s,L.scale(s,2)));
    const {H:H1,zc:zc1}=readings[3],two=constM(zc1,[[2,0],[0,2]]);
    no('T3_a_linear_congruence_invariant_of_the_metric_record_exists:'+Nd,!L.equal(trJ(L.matrixScale(H1,4)),trJ(H1))&&L.matrixEqual(L.matrixScale(H1,4),L.matrixMul(L.matrixMul(two,H1),two)));
  }

  // ---- Theorem 4: chart jets, the observer and the Laplacian tower ------------------
  const theta=(p,q)=>x=>subst(x,p,q);
  const shearS=theta(L.add(s,v),v),shearV=theta(s,L.add(v,s)),scaleD=theta(L.scale(s,2),L.scale(v,'1/2'));
  const rot=(cs,sn)=>theta(L.sub(L.scale(s,cs),L.scale(v,sn)),L.add(L.scale(s,sn),L.scale(v,cs)));
  const rot1=rot('3/5','4/5'),rot2=rot('5/13','12/13');
  const curved=theta(L.add(s,L.add(v,L.mul(s,s))),L.sub(L.scale(v,2),L.mul(s,v)));
  const autos=[shearS,shearV,scaleD,rot1,rot2,curved];
  ck('T4_substitutions_are_unital_star_homomorphisms_fixing_the_observer:'+Nd,autos.every(th=>pairs.every(([x,y])=>L.equal(th(L.mul(x,y)),L.mul(th(x),th(y)))&&L.equal(th(L.add(x,y)),L.add(th(x),th(y)))&&L.outer(th(x)).eq(L.outer(x))&&L.outer(th(y)).eq(L.outer(y)))&&L.equal(th(one),one)&&L.equal(th(L.dagger(a)),L.dagger(th(a)))));
  ck('T4_substitutions_are_automorphisms_full_rank:'+Nd,autos.every(th=>spanRank(monomials(Cx).map(e=>vecJ(th(e))))===d));
  { const inv=theta(L.sub(s,v),v);ck('T4_shear_inverse_composes_to_identity:'+Nd,pairs.every(([x])=>L.equal(inv(shearS(x)),x)&&L.equal(shearS(inv(x)),x)));
    const ps=L.add(s,L.mul(s,s)),series=[s,L.sub(s,L.mul(s,s)),L.add(L.sub(s,L.mul(s,s)),L.scale(L.power(s,3),2))][Math.min(Nd,3)-1];
    const th=theta(ps,v),thInv=theta(series,v);ck('T4_curved_chart_inverse_through_depth_N:'+Nd,pairs.every(([x])=>L.equal(thInv(th(x)),x)));}
  no('T4_singular_linear_part_is_an_automorphism:'+Nd,spanRank(monomials(Cx).map(e=>vecJ(subst(e,s,s))))<d);
  no('T4_a_substitution_moving_the_point_respects_the_quotient:'+Nd,!L.equal(L.power(L.add(one,s),Nd+1),zero)&&L.equal(L.power(s,Nd+1),zero)); // eps_s -> 1 + eps_s would send the relation eps_s^(N+1) = 0 to a nonzero element
  { const diffs=ths=>ths.flatMap(th=>monomials(Cx).map(e=>vecJ(L.sub(th(e),e))));
    ck('T4_SL2_chart_jets_leave_only_the_observer:'+Nd,spanRank(diffs([shearS,shearV,scaleD]))===d-1&&[shearS,shearV,scaleD].every(th=>monomials(Cx).every(e=>L.outer(L.sub(th(e),e)).zero())));
    const lap=k=>x=>{let z=N.ZERO;for(let j=0;j<=k;j++)z=z.add(C(binom(k,j)).mul(C(fact(2*j)*fact(2*k-2*j))).mul(L.coefficient(x,2*j,2*k-2*j)));return z;};
    const laps=Array.from({length:Math.floor(Nd/2)+1},(_,k)=>lap(k));
    ck('T4_rotations_leave_exactly_the_Laplacian_tower:'+Nd,spanRank(diffs([rot1,rot2]))===d-laps.length&&[rot1,rot2].every(th=>monomials(Cx).every(e=>laps.every(l=>l(L.sub(th(e),e)).zero()))));
    ck('T4_Laplacian_readings_are_independent_of_the_observer:'+Nd,spanRank(laps.map(l=>monomials(Cx).map(l)))===laps.length);
    if(Nd>=2)no('T4_rotation_invariant_Laplacian_reading_survives_a_shear:'+Nd,!L.equal(scaleD(L.mul(s,s)),L.mul(s,s))&&!laps[1](L.sub(scaleD(L.mul(s,s)),L.mul(s,s))).zero());
    if(Nd>=2)no('T4_a_compact_frame_group_collapses_the_dual_to_the_observer:'+Nd,spanRank(diffs([rot1,rot2]))<d-1);
  }
  // the coefficient functionals are the derivative distributions of the observer (finite part of the rigging)
  ck('T4_layer_coefficients_are_point_derivatives:'+Nd,L.basisIndices(Cx).every(([p,q])=>[a,b,c,S.d].every(x=>{let y=x;for(let k=0;k<p;k++)y=L.derivative(y,'s');for(let k=0;k<q;k++)y=L.derivative(y,'v');return L.outer(y).eq(L.coefficient(x,p,q).mul(C(fact(p)*fact(q))));})));
}

// ---- Theorem 5: complete finite towers versus finite transport --------------------
for(const Nd of [1,2,3]){
  // TET-T6 family H_eta = diag(1 + eta t^(N+1), 1) on the s-axis: equal transport jets through N, exact rational endpoint differs.
  const eta=f('5/4').mul(f(2).pow(Nd+1));                 // chosen so that 1 + eta/2^(N+1) = 9/4
  const Cp=T.createContext(Nd+1,'JED TET-T6 path');
  const H=e=>T.series(Cp,Array.from({length:Nd+2},(_,k)=>k===0?N.identity(2):k===Nd+1?N.matrix([[e,0],[0,0]]):N.zeros(2)));
  const U=e=>T.transport(T.connectionFromMetric(H(e)));
  const Ue=U(eta),U0=U(0);
  ck('T5_dependency_TET_T6_equal_transport_jets_through_N:'+Nd,Array.from({length:Nd+1},(_,k)=>k).every(k=>N.equal(T.coefficient(Ue,k),T.coefficient(U0,k))));
  ck('T5_dependency_TET_T6_first_difference_is_minus_eta_half:'+Nd,N.equal(T.coefficient(Ue,Nd+1),N.matrix([[eta.neg().div(2),0],[0,0]]))&&N.isZero(T.coefficient(U0,Nd+1)));
  const endpoint=f(1).add(eta.div(f(2).pow(Nd+1)));
  ck('T5_exact_endpoint_is_diag_two_thirds_one:'+Nd,endpoint.eq('9/4')&&f('2/3').pow(2).mul(endpoint).eq(1));
  no('T5_truncated_tower_evaluated_at_the_endpoint_gives_the_transport:'+Nd,(()=>{let acc=N.zeros(2);for(let k=0;k<=Nd;k++)acc=N.add(acc,N.scale(T.coefficient(Ue,k),f(1,1).div(f(2).pow(k))));return N.equal(acc,N.identity(2))&&!N.equal(acc,N.matrix([['2/3',0],[0,1]]));})());
}
{ // analytic positive control: on TET-T4's affine path the retained tower is the Taylor tower of the closed form (binomial series), order by order
  const z=L.buildFixture('5',4),path=T.sourcePath(z,{a:0,b:1});
  ck('T5_analytic_affine_path_tower_matches_closed_form',T.equal(T.project(path.U,path.H.context.order),T.affineBinomialTransport(path.H)));
  // the flat-function obstruction of Theorem 5(a) has no rational data; nothing here evaluates it, and no finite check could.
}

// ---- certificate -------------------------------------------------------------------
const sha=x=>crypto.createHash('sha256').update(x).digest('hex');
const payload={protocol:'JED1_JET_ENTRY_DUALS_V1',status:'PASS_FINITE_WITNESSES_FOR_JET_ENTRY_DUALS',
 source_pins_sha256:sha(fs.readFileSync(path.join(__dirname,'SOURCE_PINS.json'))),checks,negative_controls,
 results:{continuum_proofs_machine_certified:false,
  theorem_1:'written proof in JET_ENTRY_DUALS.md: every irreducible representation of M_m(J_N) factors through the degree-zero observer',
  theorem_2:'written proof: every positive functional on M_m(J_N) factors through the observer; no state reads a layer',
  theorem_3:'written proof: the frame-invariant linear dual is exactly the layer traces (dimension d_N)',
  theorem_4:'written proof: under chart jets containing SL(2) the invariant linear dual is the observer alone; under rotations alone the Laplacian tower survives',
  theorem_5:'written proof plus TET-T6 dependency: the complete smooth tower does not determine finite transport; the analytic tower does'},
 summary:{exact_checks:checks.length,negative_controls:negative_controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',hash=sha(bytes),cp=path.join(__dirname,'CERTIFICATE.json'),hp=path.join(__dirname,'EXPECTED.sha256');
if(process.argv.includes('--write')){fs.writeFileSync(cp,bytes);fs.writeFileSync(hp,hash+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(cp,'utf8'),bytes);assert.equal(fs.readFileSync(hp,'utf8'),hash+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,...payload.summary,sha256:hash}));
