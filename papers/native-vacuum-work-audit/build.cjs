'use strict';
const fs=require('node:fs'),path=require('node:path');
const g=(name,q)=>[name,[q]],turn=(name,q)=>[g('h',q),g(name,q),g('h',q)];
const cases=[
 {id:'idle',q:0,theta_pi:'0',operations:[]},
 {id:'ground_phase',q:1,theta_pi:'0',operations:[g('p16',1)]},
 {id:'drive_s',q:2,theta_pi:'1/2',operations:turn('s',2)},
 {id:'return_s',q:3,theta_pi:'0',operations:[...turn('s',3),...turn('sdg',3)]},
 {id:'drive_t',q:4,theta_pi:'1/4',operations:turn('t',4)},
 {id:'return_t',q:5,theta_pi:'0',operations:[...turn('t',5),...turn('tdg',5)]},
 {id:'drive_p16',q:6,theta_pi:'1/8',operations:turn('p16',6)},
 {id:'return_p16',q:7,theta_pi:'0',operations:[...turn('p16',7),...turn('p16dg',7)]},
 {id:'opposite',q:8,theta_pi:'-1/2',operations:turn('sdg',8)},
 {id:'wrong_return',q:9,theta_pi:'1',operations:[...turn('s',9),...turn('s',9)]}
];
const request={schema:'nrqc-calculator-v1',n_qubits:10,backend:'cyclotomic64',decimal_places:24,
 operations:cases.flatMap(c=>c.operations),queries:cases.flatMap(c=>['Z','Y'].map(axis=>({id:c.id+'_'+axis,kind:'expectation',paulis:{[c.q]:axis}})))};
request.queries.push({id:'normalization',kind:'expectation',paulis:{}});
for(const c of cases.filter(c=>c.id.startsWith('drive_')||c.id==='opposite'||c.id==='wrong_return'))request.queries.push({id:c.id+'_excited',kind:'probability',outcomes:{[c.q]:1}});
const write=(n,v)=>fs.writeFileSync(path.join(__dirname,n),JSON.stringify(v,null,2)+'\n');
write('request.json',request);write('cases.json',{hamiltonian:'E_off I + Delta (I-Z)/2',energy_unit:'Delta; supplied calibration model',offsets:['0','1/2','100'],cases});
console.log(JSON.stringify({cases:cases.length,queries:request.queries.length}));
