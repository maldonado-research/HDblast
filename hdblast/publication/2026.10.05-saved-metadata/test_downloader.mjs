import fs from "node:fs";
import path from "node:path";
import {fileURLToPath} from "node:url";
import assert from "node:assert/strict";
import {createHash} from "node:crypto";
import {Sha256,Crc32,BUNDLE_NAME,COMMITS,validateSourceUrl,validateManifest,zipSize,
  buildStoredZip,prepareBundle,hashBlob} from "../../../site/zenodo-files/core.mjs";
import {mount} from "../../../site/zenodo-files/app.mjs";
const ROOT=path.dirname(fileURLToPath(import.meta.url));
const digest=bytes=>createHash("sha256").update(bytes).digest("hex");
const cases=[];
async function check(name,run) {await run();cases.push({name,status:"PASS"});}
const originalBytes=new Map();
const files=[];
for(let index=0;index<13;index++) {
  const filename="file-"+String(index).padStart(2,"0")+".txt";
  const bytes=new TextEncoder().encode("Exact bytes for "+filename+". "+"x".repeat(index*17));
  const count=index===2?3:(index===5 || index===9)?2:1,parts=[];
  for(let part=0;part<count;part++) {
    const chunk=bytes.slice(Math.floor(bytes.length*part/count),Math.floor(bytes.length*(part+1)/count));
    const url="https://raw.githubusercontent.com/maldonado-research/HDblast/"+COMMITS[index===12?1:0]+
      "/research/test-fixture/"+filename+".part"+part;
    parts.push({url,bytes:chunk.length,sha256:digest(chunk)});originalBytes.set(url,chunk);
  }
  files.push({filename,bytes:bytes.length,sha256:digest(bytes),parts});
}
const fixtureEntries=files.map(file=>{
  const bytes=new Uint8Array(file.bytes);let at=0;
  for(const part of file.parts) {bytes.set(originalBytes.get(part.url),at);at+=part.bytes;}
  return {filename:file.filename,blob:new Blob([bytes]),crc32:new Crc32().update(bytes).digest()};
});
const fixtureZip=buildStoredZip(fixtureEntries);
const fixture={record_id:23114217,total_bytes:files.reduce((n,f)=>n+f.bytes,0),files,
  zip:{filename:BUNDLE_NAME,bytes:zipSize(files),sha256:digest(new Uint8Array(await fixtureZip.arrayBuffer()))}};
function mock({alter,header,chunkSize=11}={}) {
  const state={requests:[],cancelled:0};
  state.fetchImpl=async(url,options)=>{
    state.requests.push({url,options});
    let bytes=originalBytes.get(url);
    assert.ok(bytes);assert.equal(options.method,"GET");assert.equal(options.mode,"cors");
    assert.equal(options.credentials,"omit");assert.equal(options.redirect,"error");
    assert.equal(options.referrerPolicy,"no-referrer");assert.equal(options.cache,"no-store");
    assert.ok(!options.headers && !options.body);
    if(alter)bytes=alter(bytes,url,state.requests.length);
    let at=0;
    const body=new ReadableStream({
      pull(controller) {
        if(at===bytes.length)controller.close();
        else {controller.enqueue(bytes.slice(at,at+chunkSize));at=Math.min(bytes.length,at+chunkSize);}
      },
      cancel() {state.cancelled++;}
    });
    return {status:200,url,redirected:false,type:"cors",body,...header};
  };
  return state;
}
await check("SHA256 standard and boundary vectors match native crypto",()=>{
  const vectors=[new Uint8Array(),new TextEncoder().encode("abc"),new Uint8Array(1000000).fill(97)];
  for(const length of [1,55,56,63,64,65,127,128,129,4097,131077])
    vectors.push(Uint8Array.from({length},(_,i)=>(i*29+11)&255));
  for(const bytes of vectors)for(const chunkSize of [1,7,64,4096]) {
    const hash=new Sha256();
    for(let at=0;at<bytes.length;at+=chunkSize)hash.update(bytes.subarray(at,at+chunkSize));
    assert.equal(hash.digestHex(),digest(bytes));assert.throws(()=>hash.update(new Uint8Array()));
  }
  assert.equal(new Sha256().update(new TextEncoder().encode("abc")).digestHex(),
    "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
});
await check("CRC32 known vector and chunk continuation",()=>{
  const crc=new Crc32().update(new TextEncoder().encode("1234")).update(new TextEncoder().encode("56789"));
  assert.equal(crc.digest(),0xcbf43926);assert.equal(new Crc32().digest(),0);
});
await check("ordinary stored ZIP headers and determinism",async()=>{
  const bytes=new Uint8Array(await fixtureZip.arrayBuffer()),view=new DataView(bytes.buffer);
  assert.equal(view.getUint32(0,true),0x04034b50);assert.equal(view.getUint16(6,true),0x0800);
  assert.equal(view.getUint16(8,true),0);assert.equal(view.getUint16(10,true),0);
  assert.equal(view.getUint16(12,true),0x21);assert.equal(view.getUint32(bytes.length-22,true),0x06054b50);
  assert.equal(view.getUint16(bytes.length-14,true),13);
  assert.deepEqual(new Uint8Array(await buildStoredZip(fixtureEntries).arrayBuffer()),bytes);
});
await check("raw-host, repository, immutable-commit and URL boundaries",()=>{
  const valid=files[0].parts[0].url;validateSourceUrl(valid);
  for(const bad of [valid.replace("https:","http:"),valid.replace("raw.githubusercontent.com","example.com"),
    valid.replace("/HDblast/","/OtherRepo/"),valid.replace(COMMITS[0],"main"),
    valid.replace(COMMITS[0],"a".repeat(40)),valid+"?access_token=fake",valid+"#fragment",
    valid.replace("https://","https://user@"),valid.replace("/research/","/research/../private/")])
    assert.throws(()=>validateSourceUrl(bad));
});
await check("all 17 mocked transports, split reconstruction and verified ZIP",async()=>{
  const state=mock(),events=[];
  const result=await prepareBundle(fixture,{fetchImpl:state.fetchImpl,onProgress:e=>events.push(e)});
  assert.equal(state.requests.length,17);assert.equal(result.sha256,fixture.zip.sha256);
  assert.deepEqual(new Uint8Array(await result.blob.arrayBuffer()),new Uint8Array(await fixtureZip.arrayBuffer()));
  assert.equal(events.filter(e=>e.phase==="file_verified").length,13);
  assert.ok(events.some(e=>e.phase==="download" && e.part_count===3));
  assert.ok(events.some(e=>e.phase==="verify_zip" && e.received===fixture.zip.bytes));
});
await check("oversized response stops early and cancels the reader",async()=>{
  const state=mock({chunkSize:10000,alter:bytes=>new Uint8Array(bytes.length+200000)});
  await assert.rejects(prepareBundle(fixture,{fetchImpl:state.fetchImpl}),/exceeded.*byte limit/);
  assert.equal(state.requests.length,1);assert.equal(state.cancelled,1);
});
await check("truncated response blocks preparation",async()=>{
  const state=mock({alter:bytes=>bytes.slice(0,-1)});
  await assert.rejects(prepareBundle(fixture,{fetchImpl:state.fetchImpl}),/shorter/);
  assert.equal(state.requests.length,1);
});
await check("same-size corrupted transport blocks preparation",async()=>{
  const state=mock({alter:bytes=>{const copy=bytes.slice();copy[0]^=1;return copy;}});
  await assert.rejects(prepareBundle(fixture,{fetchImpl:state.fetchImpl}),/Transport piece SHA256/);
  assert.equal(state.requests.length,1);
});
await check("valid pieces with incorrect complete-file SHA block preparation",async()=>{
  const changed=structuredClone(fixture);changed.files[2].sha256="f".repeat(64);
  await assert.rejects(prepareBundle(changed,{fetchImpl:mock().fetchImpl}),/Complete file SHA256/);
});
await check("complete ZIP hash mismatch offers no result",async()=>{
  const changed=structuredClone(fixture);changed.zip.sha256="f".repeat(64);
  await assert.rejects(prepareBundle(changed,{fetchImpl:mock().fetchImpl}),/Delivery ZIP SHA256/);
});
for(const [name,header] of [
  ["HTTP failure",{status:404}],["unexpected response URL",{url:"https://example.com/file"}],
  ["redirected response",{redirected:true}],["opaque response",{type:"opaque"}]]) {
  await check(name+" stops and cancels the response",async()=>{
    const state=mock({header});
    await assert.rejects(prepareBundle(fixture,{fetchImpl:state.fetchImpl}),/Source request/);
    assert.equal(state.requests.length,1);assert.equal(state.cancelled,1);
  });
}
await check("invalid manifest fails before requests",async()=>{
  for(const mutate of [
    m=>m.files[0].filename="../escape",m=>m.files[0].filename=undefined,
    m=>m.files[1].filename=m.files[0].filename,m=>m.files[0].parts[0].url="https://example.com/x",
    m=>m.files[0].bytes=2**32,m=>m.files[0].parts[0].bytes++,
    m=>m.files[0].parts[0].sha256="bad",m=>m.zip.bytes++,m=>m.zip.sha256=null,
    m=>m.files[1].parts[0].url=m.files[0].parts[0].url]) {
    const changed=structuredClone(fixture),state=mock();mutate(changed);
    await assert.rejects(prepareBundle(changed,{fetchImpl:state.fetchImpl}));assert.equal(state.requests.length,0);
  }
});
await check("ZIP builder rejects paths, duplicate names and oversized names",()=>{
  for(const name of ["../escape","absolute/name","back\\slash","a".repeat(65536)])
    assert.throws(()=>buildStoredZip([{...fixtureEntries[0],filename:name}]));
  assert.throws(()=>buildStoredZip([fixtureEntries[0],fixtureEntries[0]]));
});
await check("cancellation before requests is read-only",async()=>{
  const state=mock(),controller=new AbortController();controller.abort();
  await assert.rejects(prepareBundle(fixture,{fetchImpl:state.fetchImpl,signal:controller.signal}),{name:"AbortError"});
  assert.equal(state.requests.length,0);
});
await check("cancellation during collection cancels the reader",async()=>{
  const state=mock(),controller=new AbortController();
  await assert.rejects(prepareBundle(fixture,{fetchImpl:state.fetchImpl,signal:controller.signal,
    onProgress:e=>{if(e.phase==="download")controller.abort();}}),{name:"AbortError"});
  assert.equal(state.requests.length,1);assert.equal(state.cancelled,1);
});
await check("cancellation during final ZIP verification blocks its result",async()=>{
  const controller=new AbortController();
  await assert.rejects(prepareBundle(fixture,{fetchImpl:mock().fetchImpl,signal:controller.signal,
    onProgress:e=>{if(e.phase==="verify_zip")controller.abort();}}),{name:"AbortError"});
});
await check("streamed Blob SHA equals native crypto",async()=>{
  assert.equal(await hashBlob(fixtureZip),fixture.zip.sha256);
});
function dom() {
  const elements=new Map();
  for(const id of ["prepare","cancel","download","status","progress","detail"]) {
    elements.set(id,{hidden:["cancel","download","progress"].includes(id),disabled:false,textContent:"",
      listeners:{},addEventListener(event,handler){this.listeners[event]=handler;},
      removeAttribute(name){delete this[name];},emit(event){return this.listeners[event]();}});
  }
  return {getElementById:id=>elements.get(id)};
}
function urlMock() {
  return {created:0,revoked:0,createObjectURL(){this.created++;return "blob:verified-local-test";},
    revokeObjectURL(){this.revoked++;}};
}
const ready={blob:fixtureZip,sha256:fixture.zip.sha256,files:13};
await check("UI enables Download only after verified completion",async()=>{
  const document=dom(),url=urlMock();let finish;
  mount(document,{manifest:fixture,urlAPI:url,prepare:()=>new Promise(resolve=>{finish=resolve;})});
  const pending=document.getElementById("prepare").emit("click");
  assert.equal(document.getElementById("download").hidden,true);assert.equal(url.created,0);
  finish(ready);await pending;
  assert.equal(document.getElementById("download").hidden,false);assert.equal(url.created,1);
  assert.equal(document.getElementById("download").download,BUNDLE_NAME);
});
await check("UI integrity failure offers no Blob URL",async()=>{
  const document=dom(),url=urlMock();
  mount(document,{manifest:fixture,urlAPI:url,prepare:async()=>{throw new Error("Integrity failed.");}});
  await document.getElementById("prepare").emit("click");
  assert.equal(url.created,0);assert.equal(document.getElementById("download").hidden,true);
  assert.equal(document.getElementById("status").textContent,"Preparation stopped.");
});
await check("UI cancellation rejects a late resolved preparation",async()=>{
  const document=dom(),url=urlMock();let finish;
  mount(document,{manifest:fixture,urlAPI:url,prepare:()=>new Promise(resolve=>{finish=resolve;})});
  const pending=document.getElementById("prepare").emit("click");document.getElementById("cancel").emit("click");
  finish(ready);await pending;
  assert.equal(url.created,0);assert.equal(document.getElementById("download").hidden,true);
  assert.equal(document.getElementById("status").textContent,"Preparation cancelled.");
});
await check("new preparation revokes the preceding download",async()=>{
  const document=dom(),url=urlMock();let count=0;
  mount(document,{manifest:fixture,urlAPI:url,prepare:async()=>{if(++count===1)return ready;throw new Error("Integrity failed.");}});
  await document.getElementById("prepare").emit("click");await document.getElementById("prepare").emit("click");
  assert.equal(url.created,1);assert.equal(url.revoked,1);assert.equal(document.getElementById("download").hidden,true);
});
validateManifest(fixture);
const output=path.resolve(ROOT,process.argv[2] || "NODE_CHECKS.json");
assert.ok(output!==fileURLToPath(import.meta.url),"Receipt must not replace the test source.");
const receipt={status:"PASS_OFFLINE_DOWNLOADER_CHECKS",cases:cases.length,results:cases,
  actual_browser_execution:false,live_network_requests:0,
  source_hashes:Object.fromEntries(["core.mjs","app.mjs","manifest.mjs","index.html","style.css"].map(name=>[
    name,digest(fs.readFileSync(path.resolve(ROOT,"../../../site/zenodo-files",name)))]))};
fs.writeFileSync(output,JSON.stringify(receipt,null,2)+"\n",{flag:"wx"});
console.log(JSON.stringify({status:receipt.status,cases:receipt.cases,live_network_requests:0}));
