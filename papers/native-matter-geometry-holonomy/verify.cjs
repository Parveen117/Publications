'use strict';
// MGH-1: exact metric fixtures and quantum operations on the unchanged native engines.
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto'),assert=require('node:assert/strict');
const root=__dirname,args=process.argv.slice(2),nativeRoot=args[args.indexOf('--native-root')+1],write=args.includes('--write');
assert.ok(args.includes('--native-root')&&nativeRoot);assert.ok(write!==args.includes('--check'));
const hash=x=>crypto.createHash('sha256').update(x).digest('hex');
const pins=JSON.parse(fs.readFileSync(path.join(root,'SOURCE_PINS.json'),'utf8'));
for(const[n,h]of Object.entries(pins.native_engine))assert.equal(hash(fs.readFileSync(path.join(nativeRoot,n))),h,n);
for(const[n,h]of Object.entries(pins.packet_files))assert.equal(hash(fs.readFileSync(path.join(root,n))),h,n);
for(const[n,h]of Object.entries(pins.predecessor_files))assert.equal(hash(fs.readFileSync(path.resolve(root,'../..',n))),h,n);
const N=require(path.join(nativeRoot,'native_operator.cjs')),M=require(path.join(nativeRoot,'paninian_operator.cjs'));
const f=N.F.of,C=N.Cut.of,mm=(...a)=>a.reduce(N.mul),tr=a=>a.reduce((s,r,i)=>s.add(r[i]),N.ZERO),real=z=>{assert.ok(z.turn.zero());return z.rad;};
const col=a=>N.matrix(a.map(x=>[x])),cols=a=>N.dagger(a.map(v=>v.map(r=>r[0]))),outer=(a,b=a)=>mm(a,N.dagger(b));
const dot=(a,b)=>real(mm(N.dagger(a),b)[0][0]),conj=(u,r)=>mm(u,r,N.dagger(u)),positive=x=>!f(x).le(0),cpow=(z,n)=>{let a=N.ONE;for(let j=0;j<n;j++)a=a.mul(z);return a;};
const checks=[],controls=[],proofs=[],presentations=[],systems=new Map(),results={};
const ck=(name,value,detail={})=>{assert.ok(value,name);checks.push({name,passed:true,...detail});};
const no=(name,value,detail={})=>{assert.ok(value,name);controls.push({name,rejected:true,...detail});};
const sum=a=>a.reduce(N.add,N.zeros(a[0].length,a[0][0].length));
function det(a){a=N.matrix(a);if(a.length===1)return a[0][0];return a[0].reduce((s,z,j)=>s.add(z.mul(det(a.slice(1).map(r=>r.filter((_,k)=>k!==j)))).mul(j%2?-1:1)),N.ZERO);}
const positiveGram=g=>g.every((_,k)=>positive(real(det(g.slice(0,k+1).map(r=>r.slice(0,k+1))))));

// 1. Reuse N4H-1's whole-tetrahedron unfolding construction, with a new
// curved metric on precisely its certified six-cell triangulation.
const I4=N.identity(4),e=i=>col([0,1,2,3].map(j=>i===j?1:0)),Ph=N.matrix([[1,0,0,0],[0,1,0,0],[0,0,0,0],[0,0,0,0]]),Qh=N.sub(I4,Ph),Rh=N.matrix([[0,0,0,0],[0,0,0,0],[0,0,0,-1],[0,0,1,0]]);
const pair=t=>{t=f(t);const d=f(1).add(t.pow(2));return{c:f(1).sub(t.pow(2)).div(d),s:t.mul(2).div(d)};};
const normal=(c,s)=>N.add(Ph,N.add(N.scale(Qh,c),N.scale(Rh,s)));
function star(name,parameters,scale=1){
 const L=f(scale),count=parameters.length,cells=[],edges=new Map(),transports=[];
 for(let i=0;i<count;i++){
  const{c,s}=pair(parameters[i]),coords={A:col([0,0,0,0]),B:N.scale(e(0),L),C:N.scale(e(1),L)};
  coords['D'+i]=N.scale(e(2),L);coords['D'+((i+1)%count)]=N.scale(col([0,0,c,s]),L);
  const names=['A','B','C','D'+i,'D'+((i+1)%count)],basis=cols(names.slice(1).map(v=>coords[v])),gram=mm(N.dagger(basis),basis),volume=real(det(basis)).div(24);
  ck(name+':cell_'+i+':positive_geometry',positiveGram(gram)&&volume.eq(L.pow(4).mul(s).div(24))&&positive(volume));
  for(let a=0;a<5;a++)for(let b=0;b<a;b++){
   const key=[names[a],names[b]].sort().join('|'),length2=N.energy(N.sub(coords[names[a]],coords[names[b]]));
   if(edges.has(key))ck(name+':cell_'+i+':shared_edge_'+key,edges.get(key).eq(length2));else edges.set(key,length2);
  }
  cells.push({names,coords,c,s,gram,volume});
 }
 for(let i=0;i<count;i++){
  const source=cells[i],target=cells[(i+1)%count],shared=['B','C','D'+((i+1)%count)],E=cols(shared.map(v=>source.coords[v])),F=cols(shared.map(v=>target.coords[v])),gram=mm(N.dagger(E),E);
  ck(name+':face_'+i+':matching_Gram',N.equal(gram,mm(N.dagger(F),F))&&N.rank(E)===3);
  const PE=mm(E,N.inverse(gram),N.dagger(E)),PF=mm(F,N.inverse(gram),N.dagger(F));
  const ns=N.scale(mm(N.sub(I4,PE),source.coords['D'+i]),f(1).div(L.mul(source.s))),nt=N.scale(mm(N.sub(I4,PF),target.coords['D'+((i+2)%count)]),f(1).div(L.mul(target.s)));
  const T=N.sub(mm(F,N.inverse(gram),N.dagger(E)),outer(nt,ns));
  ck(name+':face_'+i+':unit_inward_normals',dot(ns,ns).eq(1)&&dot(nt,nt).eq(1)&&N.isZero(mm(N.dagger(E),ns))&&N.isZero(mm(N.dagger(F),nt)));
  ck(name+':face_'+i+':complete_unfolding',N.equal(mm(T,E),F)&&N.equal(mm(T,ns),N.scale(nt,-1))&&N.rank(cols([e(0),e(1),source.coords[shared[2]],ns]))===4);
  ck(name+':face_'+i+':oriented_isometry',N.equal(mm(N.dagger(T),T),I4)&&det(T).eq(1));
  ck(name+':face_'+i+':transport_from_metric',N.equal(T,normal(source.c,source.s.neg())));
  if(i===0){const folded=N.add(mm(F,N.inverse(gram),N.dagger(E)),outer(nt,ns));no(name+':folded_face_rejected',N.equal(mm(folded,E),F)&&!N.equal(mm(folded,ns),N.scale(nt,-1)));}
  transports.push(T);
 }
 let H=I4;for(const t of transports)H=mm(t,H);
 const z=cells.reduce((u,c)=>u.mul(C([c.c,c.s.neg()])),N.ONE),defect=f(1).sub(z.rad),volume=cells.reduce((s,c)=>s.add(c.volume),f(0));
 ck(name+':metric_holonomy',N.equal(H,normal(z.rad,z.turn))&&N.energy(N.sub(H,I4)).div(4).eq(defect));
 const eigenCut=N.scale(N.sub(Qh,N.scale(Rh,N.IOTA)),'1/2');
 ck(name+':normal_character',N.equal(mm(eigenCut,eigenCut),eigenCut)&&tr(eigenCut).eq(1)&&N.equal(mm(H,eigenCut),N.scale(eigenCut,z))&&tr(mm(eigenCut,H)).eq(z));
 const frames=cells.map((_,i)=>{const a=N.zeros(4),u=i%4,v=(i+1)%4;a[u][v]=N.ONE;a[v][u]=N.ONE.neg();return N.cayley(a,f(i+1).div(i+2));});
 let transformed=I4;for(let i=0;i<count;i++)transformed=mm(frames[(i+1)%count],transports[i],N.dagger(frames[i]),transformed);
 ck(name+':independent_frame_covariance',N.equal(transformed,conj(frames[0],H)));
 ck(name+':character_is_frame_invariant',tr(mm(conj(frames[0],eigenCut),transformed)).eq(z));
 return{name,parameters:parameters.map(String),scale:L,cells,edge_lengths_squared:Object.fromEntries(edges),transports,holonomy:H,character:z,compact_defect:defect,coordinate_volume:volume,hinge_area_squared:L.pow(4).div(4)};
}
const flatParameters=['1/2','1/2','3/4','1/2','1/2','3/4'],curvedParameters=['1/2','1/2','2/11','1/2','1/2','2/11'];
const flat=star('flat',flatParameters),curved=star('curved',curvedParameters),scaled=star('curved_scaled',curvedParameters,3);
const predecessor=JSON.parse(fs.readFileSync(path.resolve(root,'../native-four-dimensional-hinge/CERTIFICATE.json'),'utf8'));
const oldFlat=predecessor.results.hinge_stars.find(s=>s.name==='six_wedges_flat');
ck('lineage:flat_geometry_reproduced',N.equal(flat.holonomy,N.matrix(oldFlat.holonomy))&&flat.parameters.every((t,i)=>t===oldFlat.parameters[i]));
ck('geometry:common_triangulation',JSON.stringify(flat.cells.map(c=>c.names))===JSON.stringify(curved.cells.map(c=>c.names)));
ck('geometry:common_hinge_area',flat.hinge_area_squared.eq(curved.hinge_area_squared)&&flat.hinge_area_squared.eq('1/4'));
ck('geometry:distinct_boundary_metric',!flat.edge_lengths_squared['D2|D3'].eq(curved.edge_lengths_squared['D2|D3']));
ck('geometry:scale_invariant_transport',N.equal(curved.holonomy,scaled.holonomy)&&scaled.coordinate_volume.eq(curved.coordinate_volume.mul(81)));
ck('geometry:exact_volume_values',flat.coordinate_volume.eq('16/75')&&curved.coordinate_volume.eq('61/375'));
const half=C(['3/5','4/5']),z=curved.character,kappa=curved.compact_defect;
ck('geometry:exact_hinge_phase',z.eq(['-7/25','24/25'])&&kappa.eq('32/25'));
ck('geometry:half_angle_phase',half.mul(half).eq(z));
const alpha=pair('1/2'),beta=pair('3/4'),betaMinusAlpha=pair('2/11');
ck('geometry:lift_triple_relations',C([alpha.c,alpha.s]).mul(C([alpha.c,alpha.s])).mul(C([beta.c,beta.s])).eq(-1)&&C([betaMinusAlpha.c,betaMinusAlpha.s]).mul(C([alpha.c,alpha.s])).eq(C([beta.c,beta.s])));
ck('geometry:positive_angle_branch',positive(alpha.c)&&positive(alpha.s)&&positive(beta.c)&&positive(beta.s)&&positive(betaMinusAlpha.c)&&positive(betaMinusAlpha.s));
results.geometry={flat,curved,scaled,flat_lift:'0',curved_lift:'2*alpha; cos(alpha)=3/5, sin(alpha)=4/5, 0<alpha<pi/2',interpolation:'Keep the four 1/2 Cayley parameters fixed and move both remaining parameters continuously from 3/4 to 2/11. Every parameter stays positive, so every simplex is nondegenerate. The entire finite region, including its differing boundary metric, belongs to the configuration register.',provenance:'The whole-tetrahedron unfolding adapter is reused from N4H-1. The two-state quantum configuration register and its mixing are explicit additional modelling choices.'};

// 2. Independent metric-configuration and matter-occupancy registers.
const I=N.identity(2),Z=N.matrix([[1,0],[0,-1]]),X=N.matrix([[0,1],[1,0]]),R=N.matrix([[0,-1],[1,0]]),Y=N.scale(R,N.IOTA),P0=N.scale(N.add(I,Z),'1/2'),P1=N.sub(I,P0),plus=N.scale(N.add(I,X),'1/2');
const B=N.scale(N.add(I,R),C(['1/2','1/2'])),curvature=N.scale(P1,kappa),volume=N.add(N.scale(P0,flat.coordinate_volume),N.scale(P1,curved.coordinate_volume));
ck('configuration:mixing_unitary',N.equal(mm(N.dagger(B),B),I));
ck('configuration:geometric_readout_cut',N.equal(N.sub(I,N.scale(curvature,f(2).div(kappa))),Z));
ck('configuration:active_mixing_reflection',N.equal(conj(B,Z),X));
const diagonal=values=>N.matrix(values.map((v,i)=>values.map((_,j)=>i===j?v:0)));
function interaction(phases){const d=phases.length,id=N.identity(d),blocks=phases.map(z=>N.add(P0,N.scale(P1,z))),D=N.zeros(2*d);for(let i=0;i<d;i++)for(let j=0;j<2;j++)for(let k=0;k<2;k++)D[2*i+j][2*i+k]=blocks[i][j][k];return{D,W:mm(N.kron(id,B),D,N.kron(id,N.dagger(B)))};}
function partial(r,da,db,keep){return keep==='a'?N.matrix(Array.from({length:da},(_,i)=>Array.from({length:da},(_,j)=>Array.from({length:db},(_,k)=>r[i*db+k][j*db+k]).reduce((s,z)=>s.add(z),N.ZERO)))):N.matrix(Array.from({length:db},(_,i)=>Array.from({length:db},(_,j)=>Array.from({length:da},(_,k)=>r[k*db+i][k*db+j]).reduce((s,z)=>s.add(z),N.ZERO))));}
const purity=r=>real(tr(mm(r,r))),pt=r=>N.matrix(Array.from({length:4},(_,i)=>Array.from({length:4},(_,j)=>r[2*Math.floor(i/2)+(j%2)][2*Math.floor(j/2)+(i%2)])));
const cut=(z,x,y=0)=>N.scale(N.add(I,N.add(N.scale(Z,z),N.add(N.scale(X,x),N.scale(Y,y)))),'1/2');
const sources=[{name:'coherent_half',state:plus},{name:'tilted',state:cut('3/5','4/5')},{name:'tilted_phase',state:cut('3/5',0,'4/5')},{name:'absent',state:P0},{name:'present',state:P1}];
const couplings=[{mass:f(0),phase:N.ONE},{mass:f('1/2'),phase:half},{mass:f(1),phase:z},{mass:f(2),phase:cpow(z,2)},{mass:f(3),phase:cpow(z,3)}],samples=[];
for(const coupling of couplings){
 const{mass,phase}=coupling,{D,W}=interaction([N.ONE,phase]),Mlocal=N.scale(P1,mass),Mjoint=N.kron(Mlocal,I),Gjoint=N.kron(Mlocal,conj(B,P1));
 ck('coupling:'+mass+':joint_unitarity',N.equal(mm(N.dagger(W),W),N.identity(4)));
 ck('coupling:'+mass+':conserved_mass_and_generator',N.isZero(N.commutator(W,Mjoint))&&N.isZero(N.commutator(W,Gjoint)));
 for(const source of sources){
  const input=N.kron(source.state,P0),out=conj(W,input),rm=partial(out,2,2,'a'),rg=partial(out,2,2,'b'),p=real(tr(mm(source.state,P1))),characteristic=N.ONE.sub(p).add(phase.mul(p)),entropy=f(1).sub(purity(rg)),pred=p.mul(f(1).sub(p)).mul(f(1).sub(phase.rad));
  const expectedGeometry=N.matrix([[C(f(1).add(characteristic.rad).div(2)),N.IOTA.mul(characteristic.turn.div(2))],[N.IOTA.mul(characteristic.turn.div(-2)),C(f(1).sub(characteristic.rad).div(2))]]);
  ck('sample:'+mass+':'+source.name+':pure_joint',N.equal(mm(out,out),out)&&tr(out).eq(1));
  ck('sample:'+mass+':'+source.name+':mass_moments',tr(mm(rm,Mlocal)).eq(tr(mm(source.state,Mlocal)))&&tr(mm(rm,N.power(Mlocal,2))).eq(tr(mm(source.state,N.power(Mlocal,2)))));
  ck('sample:'+mass+':'+source.name+':characteristic_readout',N.equal(rg,expectedGeometry));
  ck('sample:'+mass+':'+source.name+':information_curvature',entropy.eq(pred)&&f(1).sub(purity(rm)).eq(pred)&&entropy.eq(f(1).sub(characteristic.norm2()).div(2)));
  ck('sample:'+mass+':'+source.name+':metric_population',tr(mm(rg,P1)).eq(p.mul(f(1).sub(phase.rad)).div(2))&&tr(mm(rg,curvature)).eq(kappa.mul(p).mul(f(1).sub(phase.rad)).div(2)));
  samples.push({mass,source:source.name,initial_source:source.state,p,characteristic,out,source_reduced:rm,geometry_reduced:rg,linear_entropy:entropy,curved_probability:real(tr(mm(rg,P1))),mean_curvature:real(tr(mm(rg,curvature))),mean_volume:real(tr(mm(rg,volume)))});
 }
}
const{D,W}=interaction([N.ONE,z]),demo=samples.find(s=>s.mass.eq(1)&&s.source==='coherent_half');
const localCalibration=N.kron(N.add(P0,N.scale(P1,half.dagger())),N.add(P0,N.scale(P1,N.IOTA))),calibrated=conj(localCalibration,demo.out);
const expected=N.scale(N.matrix([[1,0,'3/5','4/5'],[0,0,0,0],['3/5',0,'9/25','12/25'],['4/5',0,'12/25','16/25']]),'1/2');
ck('demo:calibrated_state',N.equal(calibrated,expected));
ck('demo:exact_numbers',demo.linear_entropy.eq('8/25')&&demo.curved_probability.eq('8/25')&&demo.mean_curvature.eq('256/625')&&demo.mean_volume.eq('616/3125'));
ck('demo:local_calibration_preserves_geometric_and_mass_readings',N.isZero(N.commutator(localCalibration,N.kron(P1,I)))&&N.isZero(N.commutator(localCalibration,N.kron(I,curvature))));
const Uchat=N.add(N.kron(P0,I),N.kron(P1,N.add(N.scale(I,'3/5'),N.scale(X,N.IOTA.mul('-4/5'))))),Lmatter=N.kron(N.add(P0,N.scale(P1,half)),I);
ck('demo:chat_coupling_geometric_half_angle',N.equal(W,mm(Lmatter,Uchat)));
const Am=N.scale(N.kron(P1,X),N.IOTA.neg()),Ag=N.scale(N.kron(I,Z),N.IOTA.neg()),Omega=N.commutator(Am,Ag),mean=f('1/2'),centered=N.commutator(N.scale(N.kron(N.sub(P1,N.scale(I,mean)),X),N.IOTA.neg()),Ag);
ck('curvature:configuration_commutator_square',N.equal(mm(N.dagger(Omega),Omega),N.scale(N.kron(P1,I),4)));
ck('curvature:connected_fluctuation',tr(mm(N.kron(plus,P0),mm(N.dagger(centered),centered))).eq(1));
const sourceFrame=N.add(N.scale(I,'3/5'),N.scale(R,'4/5')),geometryFrame=N.add(N.scale(I,'5/13'),N.scale(R,'12/13')),frame=N.kron(sourceFrame,geometryFrame);
const transformedOut=conj(conj(frame,W),conj(frame,N.kron(plus,P0)));
ck('frame:joint_covariance',N.equal(transformedOut,conj(frame,demo.out)));
ck('frame:curvature_and_entanglement_invariance',tr(mm(transformedOut,conj(frame,N.kron(I,curvature)))).eq(demo.mean_curvature)&&purity(partial(transformedOut,2,2,'b')).eq(purity(demo.geometry_reduced)));
const common=N.matrix([['4/5',0],[0,'1/5']]),ea=N.add(N.scale(conj(W,N.kron(P0,P0)),'4/5'),N.scale(conj(W,N.kron(P1,P0)),'1/5')),eb=N.scale(N.add(conj(W,N.kron(cut('3/5','4/5'),P0)),conj(W,N.kron(cut('3/5','-4/5'),P0))),'1/2');
ck('ensemble:old_obstruction_pair_now_has_identical_output',N.equal(ea,eb)&&N.equal(ea,conj(W,N.kron(common,P0))));
no('mass:larger_mass_implies_global_monotonic_entanglement',samples.find(s=>s.mass.eq(3)&&s.source==='coherent_half').linear_entropy.le(samples.find(s=>s.mass.eq(2)&&s.source==='coherent_half').linear_entropy));
no('geometry:curvature_is_sufficient_without_source_uncertainty',samples.find(s=>s.mass.eq(1)&&s.source==='present').linear_entropy.zero()&&!samples.find(s=>s.mass.eq(1)&&s.source==='present').mean_curvature.zero());
results.quantum_model={configuration_labels:['full six-cell flat metric q0','full six-cell curved metric q1'],configuration_curvature:curvature,configuration_volume:volume,mixing:B,controlled_phase:D,joint_operation:W,phase_calibration:localCalibration,calibrated_demo:calibrated,samples,configuration_loop:Omega,connected_loop:centered,general_formula:'For initially pure product input rho_m tensor |q0><q0|, chi(delta)=Tr(rho_m exp(i delta M_o)); S_L=(1-|chi|^2)/2, Pr(q1)=(1-Re chi)/2. For M_o=mu Pi_o, S_L=p(1-p)[1-cos(mu delta)]. For unit mass ratio and the flat-reference hinge, S_L=p(1-p)*kappa_h.',additional_contract:'The metric configurations are represented by orthogonal native cuts, with coherent mixing B. A mass-weighted holonomy character selects the joint phase. These are explicit finite quantization and interaction choices; they are not derived solely from metric gluing.'};

// 3. Five-settings entanglement witness and a classical-mixture control.
const v=col([-1,-3,3,-1]),negativeProjector=N.scale(outer(v),'1/20'),witness=pt(negativeProjector),transposed=pt(calibrated);
ck('witness:negative_eigenvector',N.equal(mm(transposed,v),N.scale(v,'-2/5')));
ck('witness:projector_positive',N.equal(mm(negativeProjector,negativeProjector),negativeProjector)&&tr(negativeProjector).eq(1));
ck('witness:negative_expectation',tr(mm(witness,calibrated)).eq('-2/5'));
const paulis={I,X,Y,Z},terms=[{axes:['X','X'],coefficient:'-1/5',expectation:'4/5'},{axes:['X','Z'],coefficient:'-3/20',expectation:'3/5'},{axes:['Y','Y'],coefficient:'1/4',expectation:'-4/5'},{axes:['Z','X'],coefficient:'3/20',expectation:'-12/25'},{axes:['Z','Z'],coefficient:'-1/5',expectation:'16/25'}];
const expanded=terms.reduce((s,t)=>N.add(s,N.scale(N.kron(paulis[t.axes[0]],paulis[t.axes[1]]),t.coefficient)),N.scale(N.identity(4),'1/4'));
ck('witness:five_setting_expansion',N.equal(expanded,witness));
for(const term of terms)ck('witness:correlator_'+term.axes.join(''),tr(mm(calibrated,N.kron(paulis[term.axes[0]],paulis[term.axes[1]]))).eq(term.expectation));
const noisy=[];
for(const visibility of['0','1/3','5/13','1/2','1']){
 const r=f(visibility),state=N.add(N.scale(calibrated,r),N.scale(N.identity(4),f(1).sub(r).div(4))),value=real(tr(mm(witness,state))),expected=f(5).sub(r.mul(13)).div(20);
 ck('witness:visibility_'+visibility+':expectation',value.eq(expected));
 ck('witness:visibility_'+visibility+':eigenmode',N.equal(mm(pt(state),v),N.scale(v,expected)));
 noisy.push({visibility:r,witness_expectation:value,detected:positive(value.neg())});
}
no('witness:visibility_boundary_is_strictly_detected',real(tr(mm(witness,N.add(N.scale(calibrated,'5/13'),N.scale(N.identity(4),'2/13'))))).zero());
const classicalInput=N.kron(N.scale(I,'1/2'),P0),classicalOut=conj(W,classicalInput),branch0=conj(W,N.kron(P0,P0)),branch1=conj(W,N.kron(P1,P0)),classicalCalibrated=conj(localCalibration,classicalOut);
ck('classical:explicit_separable_decomposition',N.equal(classicalOut,N.scale(N.add(branch0,branch1),'1/2'))&&N.equal(branch0,N.kron(P0,partial(branch0,2,2,'b')))&&N.equal(branch1,N.kron(P1,partial(branch1,2,2,'b'))));
ck('classical:identical_geometry_marginal',N.equal(partial(classicalOut,2,2,'b'),demo.geometry_reduced));
no('classical:geometry_mixedness_proves_entanglement',!purity(classicalOut).eq(1)&&N.equal(partial(classicalOut,2,2,'b'),demo.geometry_reduced)&&positive(f(1).sub(purity(demo.geometry_reduced))));
ck('classical:witness_nonnegative',f(0).le(real(tr(mm(witness,classicalCalibrated)))));
const productSamples=[P0,P1,plus,cut(0,0,1),cut('3/5','4/5')];
for(let i=0;i<productSamples.length;i++)for(let j=0;j<productSamples.length;j++)ck('witness:product_control_'+i+'_'+j,f(0).le(real(tr(mm(witness,N.kron(productSamples[i],productSamples[j]))))));
results.witness={negative_vector:v,negative_projector:negativeProjector,operator:witness,exact_expectation:'-2/5',settings:terms,white_noise_results:noisy,detection_condition:'visibility > 5/13 for the specified white-noise family; absence of a negative witness alone does not establish separability.',classical_control:{source:N.scale(I,'1/2'),output:classicalOut,geometry_reduced:partial(classicalOut,2,2,'b'),witness_expectation:real(tr(mm(witness,classicalCalibrated))),decomposition:[{weight:'1/2',product:branch0},{weight:'1/2',product:branch1}]}};

// 4. General local energy: equal first two moments do not close a finite
// phase response. The projector source is a special case that does close.
const energy=diagonal([0,1,2,3]),phases=[0,1,2,3].map(n=>cpow(z,n)),energyInteraction=interaction(phases),momentCases=[];
for(const vector of[[1,5,0,2],[2,4,3,1]]){
 const v=col(vector),source=N.scale(outer(v),'1/30'),moments=[1,2,3,4].map(k=>real(tr(mm(source,N.power(energy,k))))),out=conj(energyInteraction.W,N.kron(source,P0)),rg=partial(out,4,2,'b'),rm=partial(out,4,2,'a'),characteristic=tr(mm(source,diagonal(phases))),entropy=f(1).sub(purity(rg));
 ck('moments:'+vector.join('_')+':pure_source',N.equal(mm(source,source),source)&&tr(source).eq(1));
 ck('moments:'+vector.join('_')+':general_characteristic_formula',entropy.eq(f(1).sub(characteristic.norm2()).div(2))&&entropy.eq(f(1).sub(purity(rm)))&&tr(mm(rg,P1)).eq(f(1).sub(characteristic.rad).div(2)));
 ck('moments:'+vector.join('_')+':conserved_energy_moments',moments.every((m,i)=>m.eq(real(tr(mm(rm,N.power(energy,i+1)))))));
 momentCases.push({vector,source,moments,characteristic,linear_entropy:entropy,curved_probability:real(tr(mm(rg,P1)))});
}
ck('moments:same_first_two_invariants',momentCases[0].moments[0].eq(momentCases[1].moments[0])&&momentCases[0].moments[1].eq(momentCases[1].moments[1]));
const variance=momentCases[0].moments[1].sub(momentCases[0].moments[0].pow(2)),fourthBracket=c=>c.moments[3].sub(c.moments[0].mul(c.moments[2]).mul(4)).add(c.moments[1].pow(2).mul(3));
ck('moments:common_variance_and_fourth_order_separation',variance.eq('461/900')&&fourthBracket(momentCases[0]).sub(fourthBracket(momentCases[1])).eq('16/25'));
no('moments:first_two_invariants_fix_finite_entanglement',!momentCases[0].linear_entropy.eq(momentCases[1].linear_entropy));
no('moments:first_two_invariants_fix_finite_metric_response',!momentCases[0].curved_probability.eq(momentCases[1].curved_probability));
const oldExcess=predecessor.results.hinge_stars.find(s=>s.name==='six_wedges_excess');
ck('lift:predecessor_compact_alias',N.equal(N.matrix(oldFlat.holonomy),N.matrix(oldExcess.holonomy)));
const lifted=interaction([N.ONE,C(-1)]),liftedOut=conj(lifted.W,N.kron(plus,P0));
ck('lift:half_mass_detects_full_turn_lift',f(1).sub(purity(partial(liftedOut,2,2,'b'))).eq('1/2'));
no('lift:compact_holonomy_fixes_fractional_mass_phase',N.equal(N.matrix(oldFlat.holonomy),N.matrix(oldExcess.holonomy))&&!N.ONE.eq(-1),{written_branch:'N4H-1 retains deficits 0 and -2*pi. At mass ratio 1/2 their selected phase characters are +1 and -1.'});
results.moments={energy,moment_cases:momentCases,general_small_deficit:'S_L=delta^2*(I2-I1^2)/2-delta^4*(I4-4*I1*I3+3*I2^2)/24+O(delta^6), for bounded finite-dimensional M_o and initially pure product input.',projector_closure:'M_o^2=mu*M_o implies all positive moments I_n=mu^n*p, and chi=1-p+p*exp(i*mu*delta).',lift_control:{compact_alias:true,half_mass_phase_pair:['+1','-1'],half_mass_linear_entropy:'1/2'}};

// 5. Serialized native word proofs; analytic and positivity arguments are
// written explicitly in THEOREM.md and are not relabelled machine proofs.
function system(name,tokens,rules=[],commute=false){if(commute)for(let j=1;j<tokens.length;j++)for(let i=0;i<j;i++)rules.push({id:'central_'+j+'_'+i,lhs:[tokens[j],tokens[i]],rhs:[[[tokens[i],tokens[j]],1]],source:'Declared commuting coefficients'});const s=new M.Presentation({tokens,rules}),w=Object.fromEntries(tokens.map(t=>[t,s.word([t])])),audit=s.audit();ck(name+':confluence',audit.status==='CONFLUENT_BY_CHECKED_DIAMONDS');presentations.push({name,spec:s.spec(),audit});systems.set(s.hash(),s);return{s,w};}
function prove(name,lhs,rhs){const s=lhs.system,out=s.reduce(lhs.minus(rhs),{witness:true});ck(name,out.normal.terms.size===0);const certificate=JSON.parse(JSON.stringify(out.certificate));ck(name+':replay',s.replay(certificate));proofs.push({name,lhs:lhs.toJSON(),rhs:rhs.toJSON(),certificate});}
{
 const{s,w:{p,u,v}}=system('controlled_phase',['p','u','v'],[{id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'Joint source/configuration cut'},{id:'inverse_uv',lhs:['u','v'],rhs:[[[],1]],source:'Conjugate unit phases'},{id:'inverse_vu',lhs:['v','u'],rhs:[[[],1]],source:'Conjugate unit phases'},{id:'u_central',lhs:['u','p'],rhs:[[['p','u'],1]],source:'Central phase'},{id:'v_central',lhs:['v','p'],rhs:[[['p','v'],1]],source:'Central phase'}]);
 const one=s.one(),d=one.plus(u.minus(one).times(p)),dd=one.plus(v.minus(one).times(p));
 prove('phase:unitarity',dd.times(d),one);
 prove('phase:preserves_source_cut',d.times(p),u.times(p));
 prove('phase:complement_unchanged',d.times(one.minus(p)),one.minus(p));
}
{
 const{s,w:{c,t,p,mu}}=system('geometric_information',['c','t','p','mu'],[{id:'circle',lhs:['t','t'],rhs:[[[],1],[['c','c'],-1]],source:'Holonomy cosine and sine'}],true),one=s.one(),a=one.minus(p).plus(p.times(c)),b=p.times(t),k=one.minus(c),variance=p.times(one.minus(p));
 prove('information:characteristic_norm',one.minus(a.times(a)).minus(b.times(b)),variance.times(k).scale(2));
 prove('information:curvature_entropy',one.minus(a.times(a)).minus(b.times(b)).scale('1/2'),variance.times(k));
 prove('geometry:population_response',one.minus(a).scale('1/2'),p.times(k).scale('1/2'));
 prove('geometry:mean_defect_response',k.times(one.minus(a)).scale('1/2'),p.times(k).times(k).scale('1/2'));
 prove('matter:local_mass_variance',mu.times(mu).times(p).minus(mu.times(p).times(mu.times(p))),mu.times(mu).times(variance));
}
{
 const{s,w:{a,b}}=system('geometric_half_angle',['a','b'],[{id:'circle',lhs:['b','b'],rhs:[[[],1],[['a','a'],-1]],source:'Half-angle phase has unit norm'}],true),one=s.one(),u=a.plus(b.scale(N.IOTA)),z=u.times(u);
 prove('half_angle:retained_amplitude',one.plus(z).scale('1/2'),u.times(a));
 prove('half_angle:transition_amplitude',one.minus(z).scale('1/2'),u.times(b).scale(N.IOTA.neg()));
 prove('half_angle:compact_defect',one.minus(a.times(a).minus(b.times(b))),b.times(b).scale(2));
}
{
 const{s,w:{m,a,x,z}}=system('configuration_curvature',['m','a','x','z'],[{id:'x_square',lhs:['x','x'],rhs:[[[],1]],source:'Configuration reflection'},{id:'z_square',lhs:['z','z'],rhs:[[[],1]],source:'Metric-label reflection'},{id:'anti',lhs:['z','x'],rhs:[[['x','z'],-1]],source:'Anticommuting configuration reflections'},{id:'am',lhs:['a','m'],rhs:[[['m','a'],1]],source:'Scalar mean'},{id:'xm',lhs:['x','m'],rhs:[[['m','x'],1]],source:'Independent matter register'},{id:'zm',lhs:['z','m'],rhs:[[['m','z'],1]],source:'Independent matter register'},{id:'xa',lhs:['x','a'],rhs:[[['a','x'],1]],source:'Scalar mean'},{id:'za',lhs:['z','a'],rhs:[[['a','z'],1]],source:'Scalar mean'}]);
 const omega=m.times(x.times(z).minus(z.times(x))).scale(-1),center=m.minus(a),connected=center.times(x.times(z).minus(z.times(x))).scale(-1);
 prove('curvature:source_second_moment_operator',omega.times(omega).scale(-1),m.times(m).scale(4));
 prove('curvature:centered_operator',connected.times(connected).scale(-1),center.times(center).scale(4));
}
{
 const{s,w:{p,mu}}=system('mass_cut_closure',['p','mu'],[{id:'cut',lhs:['p','p'],rhs:[[['p'],1]],source:'Local occupancy projector'}],true),m=mu.times(p);
 prove('moments:quadratic_closure',m.times(m),mu.times(m));
 prove('moments:cubic_closure',m.times(m).times(m),mu.times(mu).times(m));
 prove('moments:fourth_closure',m.times(m).times(m).times(m),mu.times(mu).times(mu).times(m));
}
{
 const{s,w:{v}}=system('witness_visibility',['v']);
 prove('witness:white_noise_threshold',v.scale('-2/5').plus(s.one().minus(v).scale('1/4')),s.one().scale('1/4').minus(v.scale('13/20')));
}
const chosen=proofs.find(p=>p.certificate.steps.length),sys=systems.get(chosen.certificate.presentation),bad=JSON.parse(JSON.stringify(chosen.certificate));bad.steps[0].coefficient=['997','0'];
no('integrity:corrupted_proof_rejected',(()=>{try{sys.replay(bad);return false;}catch{return true;}})());
const payload={protocol:'MGH1_NATIVE_MATTER_GEOMETRY_HOLONOMY_V1',status:'PASS_GEOMETRIC_REALIZATION_AND_FINITE_QUANTUM_DEMO',scope:'Actual four-simplex metrics and their uniquely unfolded transports; a selected coherent two-metric configuration model with a fixed mass-weighted phase interaction; exact entanglement-curvature relation, moment and lift controls, and a five-setting witness. Quantum configuration mixing and the phase-coupling prescription are stated assumptions. No physical gravitational coupling, continuum dynamics, experimental measurement or completed quantum-gravity theory is certified.',source_pins_sha256:hash(fs.readFileSync(path.join(root,'SOURCE_PINS.json'))),presentations,proofs,checks,negative_controls:controls,results,summary:{native_presentations:presentations.length,symbolic_certificates:proofs.length,serialized_replays:proofs.length,exact_checks:checks.length,negative_controls:controls.length}};
const bytes=JSON.stringify(payload,null,2)+'\n',h=hash(bytes),certificate=path.join(root,'CERTIFICATE.json'),expectedHash=path.join(root,'EXPECTED.sha256');
if(write){fs.writeFileSync(certificate,bytes);fs.writeFileSync(expectedHash,h+'  CERTIFICATE.json\n');}else{assert.equal(fs.readFileSync(certificate,'utf8'),bytes);assert.equal(fs.readFileSync(expectedHash,'utf8'),h+'  CERTIFICATE.json\n');}
console.log(JSON.stringify({status:payload.status,mode:write?'generation':'read-only replay',...payload.summary,sha256:h,geometric_phase:z,hinge_defect:kappa,linear_entropy:demo.linear_entropy,curved_probability:demo.curved_probability,mean_curvature:demo.mean_curvature,witness_expectation:'-2/5',noise_visibility_threshold:'5/13'},null,2));
