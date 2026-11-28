import assert from 'node:assert/strict';
import {batchExample} from '../src/lib/multiomics/game-batch.js';
let n=0;
for(const crossed of [true,false])for(const batch of [0,.5,2,4])for(const bio of [0,.5,1,2]){const r=batchExample(crossed,batch,bio);assert.ok(Math.abs(r.difference-(crossed?bio:batch+bio))<1e-12);assert.equal(r.points.length,crossed?16:8);if(crossed)assert.ok(r.within.every(x=>Math.abs(x-bio)<1e-12));else assert.equal(r.within,null);n++;}console.log(n+' synthetic batch scenarios passed');
