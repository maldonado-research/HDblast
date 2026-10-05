// Browser/Node byte-only downloader and deterministic ZIP_STORED builder.
// No fetched bytes are parsed, evaluated, uploaded, or published.
export const BUNDLE_NAME = "HDBLAST_ZENODO_UPLOAD_ADDITIONS_23114217.zip";
export const COMMITS = Object.freeze([
  "8aeb771e7f754b2950bfb0296c53cc3b4483b57b",
  "5a851f22c84f14b454d398cebc6a4cbfa21d19da"
]);
const LIMIT = 1024 * 1024 * 1024;
const SHA_PATTERN = /^[0-9a-f]{64}$/;
const NAME_PATTERN = /^[A-Za-z0-9][A-Za-z0-9_.-]*$/;
const encoder = new TextEncoder();
function requireValue(condition, message) { if (!condition) throw new Error(message); }
export function checkAbort(signal) {
  if (signal && signal.aborted) throw new DOMException("Preparation cancelled.", "AbortError");
}
const rotate = (x, n) => (x >>> n) | (x << (32 - n));
const K = new Uint32Array([
  0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
  0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
  0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
  0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
  0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
  0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
  0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
  0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
]);
// Incremental SHA-256 avoids making a second 425-MB ArrayBuffer for WebCrypto.
// Tested against Node crypto, standard vectors, chunk boundaries, and every full file.
export class Sha256 {
  constructor() {
    this.state = new Uint32Array([0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
      0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19]);
    this.buffer = new Uint8Array(64); this.words = new Uint32Array(64);
    this.used = 0; this.bytes = 0; this.finished = false;
  }
  block(data, offset) {
    const w = this.words;
    for (let i=0; i<16; i++) {
      const at = offset + i*4;
      w[i] = ((data[at]<<24)|(data[at+1]<<16)|(data[at+2]<<8)|data[at+3]) >>> 0;
    }
    for (let i=16; i<64; i++) {
      const a=w[i-15], b=w[i-2];
      w[i] = (w[i-16]+(rotate(a,7)^rotate(a,18)^(a>>>3))+w[i-7]+
        (rotate(b,17)^rotate(b,19)^(b>>>10))) >>> 0;
    }
    let [a,b,c,d,e,f,g,h] = this.state;
    for (let i=0; i<64; i++) {
      const first = (h+(rotate(e,6)^rotate(e,11)^rotate(e,25))+((e&f)^(~e&g))+K[i]+w[i])>>>0;
      const second = ((rotate(a,2)^rotate(a,13)^rotate(a,22))+((a&b)^(a&c)^(b&c)))>>>0;
      h=g; g=f; f=e; e=(d+first)>>>0; d=c; c=b; b=a; a=(first+second)>>>0;
    }
    const result=[a,b,c,d,e,f,g,h];
    for (let i=0; i<8; i++) this.state[i]=(this.state[i]+result[i])>>>0;
  }
  update(data) {
    requireValue(!this.finished && data instanceof Uint8Array, "Invalid SHA-256 input.");
    this.bytes += data.length;
    requireValue(Number.isSafeInteger(this.bytes), "SHA-256 input is too large.");
    let at=0;
    if (this.used) {
      const take=Math.min(64-this.used,data.length);
      this.buffer.set(data.subarray(0,take),this.used);this.used+=take;at+=take;
      if (this.used===64) {this.block(this.buffer,0);this.used=0;}
    }
    while (at+64<=data.length) {this.block(data,at);at+=64;}
    if (at<data.length) {this.buffer.set(data.subarray(at),0);this.used=data.length-at;}
    return this;
  }
  digestHex() {
    requireValue(!this.finished, "SHA-256 already finalized.");
    const tail=new Uint8Array(this.used<56 ? 64 : 128);
    tail.set(this.buffer.subarray(0,this.used));tail[this.used]=0x80;
    const view=new DataView(tail.buffer);
    view.setUint32(tail.length-8,Math.floor(this.bytes/0x20000000),false);
    view.setUint32(tail.length-4,(this.bytes*8)>>>0,false);
    for(let at=0;at<tail.length;at+=64)this.block(tail,at);
    this.finished=true;
    return Array.from(this.state,n=>n.toString(16).padStart(8,"0")).join("");
  }
}
const CRC_TABLE = new Uint32Array(256);
for (let i=0;i<256;i++) {
  let n=i;for(let bit=0;bit<8;bit++)n=(n&1)?0xedb88320^(n>>>1):n>>>1;
  CRC_TABLE[i]=n>>>0;
}
export class Crc32 {
  constructor() {this.value=0xffffffff;}
  update(bytes) {for(const byte of bytes)this.value=CRC_TABLE[(this.value^byte)&255]^(this.value>>>8);return this;}
  digest() {return (this.value^0xffffffff)>>>0;}
}
export function validateSourceUrl(value) {
  requireValue(typeof value==="string", "Missing source URL.");
  const url=new URL(value);
  const segments=url.pathname.split("/");
  requireValue(url.href===value && url.origin==="https://raw.githubusercontent.com" &&
    !url.username && !url.password && !url.search && !url.hash &&
    segments[1]==="maldonado-research" && segments[2]==="HDblast" &&
    COMMITS.includes(segments[3]) && segments[4]==="research" &&
    segments.slice(4).every(s=>s && /^[A-Za-z0-9_.-]+$/.test(s) && s!=="." && s!==".."),
    "Source URL is outside the reviewed immutable repository.");
}
export function zipSize(files) {
  return files.reduce((sum,file)=>sum+file.bytes+76+2*encoder.encode("additions/"+file.filename).length,22);
}
export function validateManifest(manifest) {
  requireValue(manifest.record_id===23114217 && manifest.zip?.filename===BUNDLE_NAME &&
    Array.isArray(manifest.files) && manifest.files.length===13, "Unexpected download manifest.");
  const names=new Set(), urls=new Set();let total=0;
  for(const file of manifest.files) {
    requireValue(typeof file.filename==="string" && NAME_PATTERN.test(file.filename) && !names.has(file.filename) &&
      Number.isSafeInteger(file.bytes) && file.bytes>0 && file.bytes<LIMIT &&
      SHA_PATTERN.test(file.sha256) && Array.isArray(file.parts) && file.parts.length>0 &&
      file.parts.length<=3, "Invalid file pin.");
    names.add(file.filename);total+=file.bytes;let bytes=0;
    for(const part of file.parts) {
      requireValue(Number.isSafeInteger(part.bytes) && part.bytes>0 && part.bytes<=file.bytes &&
        SHA_PATTERN.test(part.sha256) && !urls.has(part.url), "Invalid transport piece pin.");
      validateSourceUrl(part.url);urls.add(part.url);bytes+=part.bytes;
    }
    requireValue(bytes===file.bytes, "Transport piece lengths do not match the full file.");
  }
  requireValue(total===manifest.total_bytes && total<LIMIT &&
    manifest.zip.bytes===zipSize(manifest.files) && manifest.zip.bytes<LIMIT &&
    SHA_PATTERN.test(manifest.zip.sha256) && manifest.zip.sha256!=="0".repeat(64),
    "Invalid complete ZIP pin.");
  return manifest;
}
async function readPart(part, {signal, fetchImpl, onBytes}) {
  checkAbort(signal);validateSourceUrl(part.url);
  const response=await fetchImpl(part.url,{method:"GET",mode:"cors",credentials:"omit",
    redirect:"error",referrerPolicy:"no-referrer",cache:"no-store",signal});
  if(!(response.status===200 && response.url===part.url && !response.redirected &&
    response.type!=="opaque" && response.type!=="opaqueredirect")) {
    try {await response.body?.cancel();}catch{}
    throw new Error("Source request did not return its pinned file.");
  }
  requireValue(response.body && typeof response.body.getReader==="function", "Streaming downloads require a current browser.");
  const reader=response.body.getReader(),hash=new Sha256(),chunks=[];let size=0,done=false,yieldAt=0;
  try {
    while(true) {
      checkAbort(signal);const item=await reader.read();checkAbort(signal);
      if(item.done) {done=true;break;}
      requireValue(item.value instanceof Uint8Array, "Source returned invalid bytes.");
      size+=item.value.length;
      requireValue(size<=part.bytes, "Source exceeded its exact byte limit.");
      hash.update(item.value);onBytes(item.value);chunks.push(item.value);
      if(size>=yieldAt) {yieldAt=size+1024*1024;await new Promise(resolve=>setTimeout(resolve,0));}
    }
    requireValue(size===part.bytes, "Source was shorter than its exact byte pin.");
    requireValue(hash.digestHex()===part.sha256, "Transport piece SHA256 did not match.");
    return new Blob(chunks,{type:"application/octet-stream"});
  } finally {
    if(!done) {try {await reader.cancel();}catch{}}
    reader.releaseLock();
  }
}
export async function downloadFiles(manifest,{signal,fetchImpl=globalThis.fetch,onProgress=()=>{}}={}) {
  validateManifest(manifest);checkAbort(signal);const entries=[];let received=0;
  for(let index=0;index<manifest.files.length;index++) {
    const file=manifest.files[index],hash=new Sha256(),crc=new Crc32(),pieces=[];let fileBytes=0;
    for(let partIndex=0;partIndex<file.parts.length;partIndex++) {
      const part=file.parts[partIndex];
      pieces.push(await readPart(part,{signal,fetchImpl,onBytes:bytes=>{
        hash.update(bytes);crc.update(bytes);fileBytes+=bytes.length;received+=bytes.length;
        onProgress({phase:"download",received,total:manifest.total_bytes,filename:file.filename,
          file_index:index+1,file_count:13,part_index:partIndex+1,part_count:file.parts.length});
      }}));
    }
    checkAbort(signal);
    requireValue(fileBytes===file.bytes && hash.digestHex()===file.sha256, "Complete file SHA256 did not match.");
    const blob=new Blob(pieces,{type:"application/octet-stream"});
    requireValue(blob.size===file.bytes, "Complete file size did not match.");
    entries.push({filename:file.filename,blob,crc32:crc.digest()});
    onProgress({phase:"file_verified",filename:file.filename,file_index:index+1,file_count:13,
      received,total:manifest.total_bytes});
  }
  return entries;
}
export function buildStoredZip(entries) {
  requireValue(entries.length>0 && entries.length<65535, "Invalid ZIP member count.");
  const parts=[],central=[],names=new Set();let offset=0,centralSize=0;
  for(const entry of entries) {
    requireValue(typeof entry.filename==="string" && NAME_PATTERN.test(entry.filename) && !names.has(entry.filename) &&
      entry.blob instanceof Blob && entry.blob.size<LIMIT &&
      Number.isInteger(entry.crc32) && entry.crc32>=0 && entry.crc32<=0xffffffff, "Invalid ZIP member.");
    names.add(entry.filename);
    const name=encoder.encode("additions/"+entry.filename),size=entry.blob.size;
    requireValue(name.length<=65535, "ZIP member name is too long.");
    const local=new Uint8Array(30+name.length),l=new DataView(local.buffer);
    l.setUint32(0,0x04034b50,true);l.setUint16(4,20,true);l.setUint16(6,0x0800,true);
    l.setUint16(12,0x0021,true);l.setUint32(14,entry.crc32,true);
    l.setUint32(18,size,true);l.setUint32(22,size,true);l.setUint16(26,name.length,true);local.set(name,30);
    const directory=new Uint8Array(46+name.length),d=new DataView(directory.buffer);
    d.setUint32(0,0x02014b50,true);d.setUint16(4,20,true);d.setUint16(6,20,true);
    d.setUint16(8,0x0800,true);d.setUint16(14,0x0021,true);d.setUint32(16,entry.crc32,true);
    d.setUint32(20,size,true);d.setUint32(24,size,true);d.setUint16(28,name.length,true);
    d.setUint32(42,offset,true);directory.set(name,46);
    parts.push(local,entry.blob);central.push(directory);
    offset+=local.length+size;centralSize+=directory.length;
    requireValue(offset+centralSize+22<LIMIT, "ZIP exceeds the bounded delivery size.");
  }
  const end=new Uint8Array(22),e=new DataView(end.buffer);
  e.setUint32(0,0x06054b50,true);e.setUint16(8,entries.length,true);e.setUint16(10,entries.length,true);
  e.setUint32(12,centralSize,true);e.setUint32(16,offset,true);
  return new Blob([...parts,...central,end],{type:"application/zip"});
}
export async function hashBlob(blob,{signal,onProgress=()=>{}}={}) {
  requireValue(blob instanceof Blob && blob.size<LIMIT, "Invalid bounded ZIP.");
  requireValue(typeof blob.stream==="function", "ZIP verification requires a current browser.");
  const reader=blob.stream().getReader(),hash=new Sha256();let received=0,done=false;
  let yieldAt=0;
  try {
    while(true) {
      checkAbort(signal);const item=await reader.read();checkAbort(signal);
      if(item.done) {done=true;break;}
      received+=item.value.length;requireValue(received<=blob.size, "ZIP read exceeded its size.");
      hash.update(item.value);onProgress({phase:"verify_zip",received,total:blob.size});
      // Blob streams may resolve repeatedly as microtasks. Yield to UI/cancellation.
      if(received>=yieldAt) {yieldAt=received+1024*1024;await new Promise(resolve=>setTimeout(resolve,0));}
    }
    checkAbort(signal);requireValue(received===blob.size, "Incomplete ZIP read.");
    return hash.digestHex();
  } finally {
    if(!done) {try {await reader.cancel();}catch{}}
    reader.releaseLock();
  }
}
export async function prepareBundle(manifest,options={}) {
  const entries=await downloadFiles(manifest,options);checkAbort(options.signal);
  const blob=buildStoredZip(entries);
  requireValue(blob.size===manifest.zip.bytes, "Delivery ZIP size did not match.");
  const sha256=await hashBlob(blob,options);
  requireValue(sha256===manifest.zip.sha256, "Delivery ZIP SHA256 did not match.");
  checkAbort(options.signal);
  return {blob,sha256,filename:manifest.zip.filename,files:entries.length};
}
