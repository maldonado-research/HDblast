"use strict";
// Offline browser simulations only. No outgoing network or live credentials.
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const assert = require("node:assert/strict");
const crypto = require("node:crypto");

const root = __dirname;
const script = fs.readFileSync(path.join(root, "RESTORE_HDBLAST_METADATA.js"), "utf8");
const bodyLine = script.match(/^  const BODY_TEXT = (.+);$/m);
assert.ok(bodyLine);
const bodyText = JSON.parse(bodyLine[1]);
const wanted = JSON.parse(bodyText);
const inherited = JSON.parse(script.match(/^  const INHERITED = (.+);$/m)[1]);
const api = "https://zenodo.org/api/records/23114217/draft?expand=1";
const self = "https://zenodo.org/api/records/23114217/draft";
const testCookie = "synthetic-browser-session-csrf";
const clone = value => JSON.parse(JSON.stringify(value));
function record() {
  const entries = {};
  for (const item of inherited) entries[item.filename] = {size:item.bytes,checksum:"md5:"+item.md5,status:"completed"};
  return {id:"23114217",parent:{id:"17088132",access:{owned_by:{user:"1386319"}}},
    is_draft:true,is_published:false,revision_id:12,links:{self},
    access:{record:"public",files:"public",embargo:{active:false}},
    files:{enabled:true,count:10,total_bytes:15554644,entries},
    pids:{doi:{identifier:"10.5281/zenodo.23114217",provider:"datacite"}},metadata:{},custom_fields:{}};
}
function fill(data) {
  data.metadata = clone(wanted.metadata);
  data.custom_fields = clone(wanted.custom_fields);
  data.access = { ...clone(wanted.access), status: "open" };
  Object.assign(data.files, clone(wanted.files));
  data.metadata.resource_type.title = {en:"Software"};
  data.metadata.creators[0].role.title = {en:"Researcher"};
  data.metadata.subjects.reverse();
  data.metadata.references.reverse();
  data.metadata.related_identifiers.reverse();
}
async function run(options={}) {
  const before = record();
  if (options.before) options.before(before);
  const requests=[],alerts=[],logs=[],storage=new Map();
  if (options.priorWrite) storage.set("HDBLAST_METADATA_REPAIR_23114217_816a22d1ebfc3b0ef7b831764a37e305866c46485000b3a0cba5d8e4d7b95314", "recorded");
  let data=clone(before),gets=0,reloads=0;
  const location={origin:options.origin||"https://zenodo.org",pathname:options.pathname||"/uploads/23114217",
    reload(){reloads++;}};
  const context={window:{location},document:{cookie:"session=not-readable-to-helper; csrftoken="+testCookie},
    crypto:crypto.webcrypto,TextEncoder,Uint8Array,AbortController,setTimeout,clearTimeout,
    sessionStorage:{getItem:key=>storage.get(key)||null,setItem:(key,value)=>storage.set(key,value)},
    console:{info:(...args)=>logs.push(args)},alert:message=>alerts.push(message),
    fetch:async(url,config)=>{
      requests.push({url,config});
      assert.equal(url,api);
      assert.equal(config.credentials,"same-origin");
      assert.equal(config.mode,"same-origin");
      assert.equal(config.redirect,"error");
      assert.equal(config.cache,"no-store");
      assert.equal(config.headers["X-CSRFToken"],testCookie);
      assert.ok(!("Authorization" in config.headers) && !("Cookie" in config.headers));
      if(config.method==="GET") {
        gets++;
        const status=options.getStatus||200;
        return {status,url:options.responseUrl||api,text:async()=>JSON.stringify(data)};
      }
      assert.equal(config.method,"PUT");
      assert.equal(config.headers["If-Match"],"12");
      assert.equal(config.body,bodyText);
      if(options.putNetworkError) throw new Error("a request failed "+testCookie);
      if(!options.blankAfter && !options.putStatus) fill(data);
      data.revision_id=13;
      if(options.after) options.after(data);
      return {status:options.putStatus||200,url:api,text:async()=>JSON.stringify(data)};
    }};
  await vm.runInNewContext(options.script||script,context,{timeout:1000});
  const result=clone(context.window.HDBLAST_METADATA_REPAIR_RESULT);
  const visible=JSON.stringify({alerts,logs,result,storage:[...storage]});
  assert.ok(!visible.includes(testCookie),"Cookie value surfaced in diagnostic output");
  assert.equal(result.publication_ready,false);
  assert.equal(result.uploads,0);assert.equal(result.publications,0);assert.equal(result.deletions,0);
  assert.equal(context.window.HDBLAST_METADATA_REPAIR_IN_FLIGHT,false);
  return {result,puts:requests.filter(x=>x.config.method==="PUT").length,gets,reloads};
}
(async()=>{
  const cases=[];
  async function check(name,options,status,puts,reason) {
    const observed=await run(options);
    assert.equal(observed.result.status,status,name + " " + JSON.stringify(observed.result));
    assert.equal(observed.puts,puts,name);
    if(reason) assert.equal(observed.result.reason,reason,name);
    assert.equal(observed.reloads,status.startsWith("PASS_")?1:0,name);
    cases.push({name,status,metadata_puts:observed.puts});
  }
  await check("exact approved save and expanded/reordered readback",{},"PASS_BROWSER_SAVED_METADATA_READBACK",1);
  await check("already saved metadata is read-only",{before:fill},"PASS_ALREADY_SAVED_METADATA_READBACK",0);
  await check("approved partial title can be completed",{before:d=>{d.metadata.title=wanted.metadata.title;}},"PASS_BROWSER_SAVED_METADATA_READBACK",1);
  await check("wrong origin stops before network",{origin:"https://github.com"},"BLOCKED_NO_METADATA_WRITE",0,"OPEN_THE_EXISTING_ZENODO_DRAFT_23114217_FIRST");
  await check("wrong editor stops before network",{pathname:"/uploads/23111008"},"BLOCKED_NO_METADATA_WRITE",0,"OPEN_THE_EXISTING_ZENODO_DRAFT_23114217_FIRST");
  await check("wrong family",{before:d=>{d.parent.id="22922927";}},"BLOCKED_NO_METADATA_WRITE",0,"WRONG_DRAFT_OWNER_OR_DOI_FAMILY");
  await check("wrong owner",{before:d=>{d.parent.access.owned_by.user="1";}},"BLOCKED_NO_METADATA_WRITE",0,"WRONG_DRAFT_OWNER_OR_DOI_FAMILY");
  await check("published record",{before:d=>{d.is_published=true;}},"BLOCKED_NO_METADATA_WRITE",0,"DRAFT_STATE_IS_NOT_UNPUBLISHED");
  await check("unexpected self link",{before:d=>{d.links.self="https://example.invalid/draft";}},"BLOCKED_NO_METADATA_WRITE",0,"UNEXPECTED_DRAFT_SELF_LINK");
  await check("changed inherited checksum",{before:d=>{d.files.entries[inherited[0].filename].checksum="md5:wrong";}},"BLOCKED_NO_METADATA_WRITE",0,"INHERITED_FILE_BYTES_OR_CHECKSUM_CHANGED");
  await check("changed file membership",{before:d=>{d.files.count=11;}},"BLOCKED_NO_METADATA_WRITE",0,"INHERITED_FILE_INVENTORY_CHANGED");
  await check("owner's unreviewed title is preserved",{before:d=>{d.metadata.title="Another owner's change";}},"BLOCKED_NO_METADATA_WRITE",0,"UNREVIEWED_EXISTING_METADATA_REQUIRES_REVIEW");
  await check("prior attempt stops repeat",{priorWrite:true},"BLOCKED_NO_METADATA_WRITE",0,"EARLIER_BROWSER_WRITE_RECORDED_RECONCILE_BEFORE_RETRY");
  await check("HTTP200 with empty metadata is failure",{blankAfter:true},"BLOCKED_RECONCILE_BROWSER_WRITE_BEFORE_RETRY",1,"WRITE_READBACK_DID_NOT_VERIFY_ALL_APPROVED_METADATA");
  await check("conditional write conflict",{putStatus:412},"BLOCKED_RECONCILE_BROWSER_WRITE_BEFORE_RETRY",1,"WRITE_READBACK_DID_NOT_VERIFY_ALL_APPROVED_METADATA");
  await check("file mutation after save prevents success",{after:d=>{d.files.entries[inherited[0].filename].size++;}},"BLOCKED_RECONCILE_BROWSER_WRITE_BEFORE_RETRY",1,"INHERITED_FILE_BYTES_OR_CHECKSUM_CHANGED");
  await check("reserved DOI mutation prevents success",{after:d=>{d.pids.doi.identifier="wrong";}},"BLOCKED_RECONCILE_BROWSER_WRITE_BEFORE_RETRY",1,"RESERVED_IDENTIFIERS_CHANGED_REQUIRES_REVIEW");
  await check("network exception does not expose cookie",{putNetworkError:true},"BLOCKED_RECONCILE_BROWSER_WRITE_BEFORE_RETRY",1,"BROWSER_REQUEST_OR_CAPABILITY_ERROR");
  await check("read authentication failure",{getStatus:403},"BLOCKED_NO_METADATA_WRITE",0,"READ_HTTP_403");
  await check("unexpected response URL",{responseUrl:"https://example.invalid"},"BLOCKED_NO_METADATA_WRITE",0,"UNEXPECTED_RESPONSE_URL");
  await check("body mutation prevents all requests",{script:script.replace("HDBLAST: Conditional Source/Operator", "HDBLAST: WRONG Source/Operator")},"BLOCKED_NO_METADATA_WRITE",0,"APPROVED_METADATA_BODY_HASH_MISMATCH");
  const receipt={status:"PASS_OFFLINE_BROWSER_CONTROL_SIMULATIONS",cases:cases.length,
    actual_browser_execution:false,live_network_requests:0,live_metadata_writes:0,
    script_sha256:crypto.createHash("sha256").update(script).digest("hex"),
    metadata_body_sha256:crypto.createHash("sha256").update(bodyText).digest("hex"),results:cases};
  fs.writeFileSync(path.join(root,"OFFLINE_CHECKS.json"),JSON.stringify(receipt,null,2)+"\n",{flag:"wx"});
  console.log(JSON.stringify({status:receipt.status,cases:receipt.cases,actual_browser_execution:false,live_network_requests:0}));
})().catch(error=>{console.error("Offline browser guard simulation failed: " + String(error.stack).replaceAll(testCookie,"[REDACTED_SYNTHETIC_VALUE]"));process.exitCode=1;});
