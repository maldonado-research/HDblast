import {MANIFEST} from "./manifest.mjs";
import {prepareBundle, checkAbort} from "./core.mjs";

const amount = value => (value/1048576).toFixed(1)+" MiB";
export function mount(document, {manifest=MANIFEST, prepare=prepareBundle, urlAPI=URL}={}) {
  const start=document.getElementById("prepare"),cancel=document.getElementById("cancel");
  const link=document.getElementById("download"),status=document.getElementById("status");
  const progress=document.getElementById("progress"),detail=document.getElementById("detail");
  let active=null,objectUrl=null,lastUpdate=0;
  function clearDownload() {
    link.hidden=true;link.removeAttribute("href");
    if(objectUrl) {urlAPI.revokeObjectURL(objectUrl);objectUrl=null;}
  }
  function showProgress(value) {
    if(value.phase==="download") {
      const now=Date.now();if(now-lastUpdate<100 && value.received!==value.total)return;
      lastUpdate=now;status.textContent="Downloading file "+value.file_index+" of 13.";
      detail.textContent=value.filename+" — "+amount(value.received)+" of "+amount(value.total)+
        " collected"+(value.part_count>1?" (piece "+value.part_index+" of "+value.part_count+")":"");
    } else if(value.phase==="verify_zip") {
      status.textContent="Checking the complete ZIP.";
      detail.textContent=amount(value.received)+" of "+amount(value.total)+" checked.";
    } else {
      status.textContent="Verified file "+value.file_index+" of 13.";
      detail.textContent=value.filename;
    }
    progress.hidden=false;progress.max=value.total;progress.value=value.received;
  }
  start.addEventListener("click",async()=>{
    if(active)return;
    clearDownload();lastUpdate=0;
    const controller=new AbortController();active=controller;
    start.disabled=true;cancel.hidden=false;cancel.disabled=false;
    progress.hidden=false;progress.max=manifest.total_bytes;progress.value=0;
    status.textContent="Preparing the files.";detail.textContent="";
    try {
      const result=await prepare(manifest,{signal:controller.signal,onProgress:showProgress});
      checkAbort(controller.signal);
      if(result.sha256!==manifest.zip.sha256 || result.blob.size!==manifest.zip.bytes || result.files!==13)
        throw new Error("The verified ZIP result did not match its pins.");
      objectUrl=urlAPI.createObjectURL(result.blob);link.href=objectUrl;
      link.download=manifest.zip.filename;link.hidden=false;
      progress.max=1;progress.value=1;status.textContent="Verified ZIP ready.";
      detail.textContent="13 files. "+result.blob.size.toLocaleString()+" bytes.";
    } catch(error) {
      clearDownload();progress.hidden=true;
      status.textContent=controller.signal.aborted?"Preparation cancelled.":"Preparation stopped.";
      detail.textContent=controller.signal.aborted?"You can prepare the download again.":String(error.message || "The files could not be verified.");
    } finally {
      active=null;start.disabled=false;cancel.hidden=true;cancel.disabled=false;
    }
  });
  cancel.addEventListener("click",()=>{
    if(active) {cancel.disabled=true;status.textContent="Cancelling preparation.";active.abort();}
  });
}
if(typeof document!=="undefined") mount(document);
