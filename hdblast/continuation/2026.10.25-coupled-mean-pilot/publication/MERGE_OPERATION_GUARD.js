(async()=>{
 const H=eval(load("hash_source")),gate=load("ridge_ci_final_premerge"),review=JSON.parse(load("ridge_final_static_publication_review"));
 const need=(ok,msg)=>{if(!ok)throw Error(msg);};
 need(gate?.status==="PASS_ALL_REQUIRED_AND_ACTUAL_CI"&&gate.all_required_and_actual_success,"all exact-head premerge CI required");
 need(H.sha256(load("ridge_final_static_publication_review"))==="9679297d226f3da86e97054343d130d569d8290a53bf77c712f18ab554a8ef5c","exact static approval");
 need(review.blockers.length===0&&review.status==="APPROVED_STATIC_PUBLICATION_PENDING_EXACT_HEAD_CI","static approval status");
 need(load("ridge_final_pre_reproduction_receipt")?.native?.accuracy_target_passed===true,"successful regenerated scientific target");
 need(load("ridge_final_pre_package_receipt")?.archive_sha256===load("ridge_zip_metadata").sha256,"executed CRC and archive/source equality");
 const rows=[
  {label:"research",repository:"maldonado-research/HDblast",pr:84,base:"2fb40bdc5a1785ef6a55bc50e6c9d81a3e3e9978",head:load("ridge_final_research_head"),tree:load("ridge_final_research_tree"),branch:"research/coupled-mean-pilot-20261025",review:review.reviewed_heads.research},
  {label:"shared",repository:"maldonado-research/maldonado-research.github.io",pr:75,base:"de41839498e54544403dfc6ade9cb68da44342fd",head:load("ridge_final_shared_head"),tree:load("ridge_final_shared_tree"),branch:"continuation/hdblast-coupled-mean-pilot-20261025",review:review.reviewed_heads.shared_site}
 ];
 const get=async url=>{const r=await tools.mcp__codex_apps__github_fetch({url});need(!r.isError,"canonical GET");return JSON.parse(r.structuredContent.content);};
 const done=[];
 for(const row of rows){
  need(row.review.head===row.head&&row.review.base===row.base&&row.review.tree===row.tree,"reviewed exact repository identities");
  const ci=gate.repositories.find(r=>r.label===row.label);need(ci.ok&&ci.head===row.head&&ci.repository===row.repository,"exact successful CI role");
  const prefix="https://api.github.com/repos/"+row.repository;
  const [pr,main,commit,legacy]=await Promise.all([get(prefix+"/pulls/"+row.pr),get(prefix+"/git/ref/heads/main"),get(prefix+"/git/commits/"+row.head),get(prefix+"/commits/"+row.head+"/status?per_page=100")]);
  need(legacy.sha===row.head&&legacy.total_count===legacy.statuses.length&&legacy.statuses.every(v=>v.state==="success"),"complete successful legacy status contexts");
  need(pr.state==="open"&&!pr.draft&&pr.mergeable===true&&pr.head.sha===row.head&&pr.base.sha===row.base&&pr.base.ref==="main"&&pr.head.ref===row.branch&&pr.head.repo.full_name===row.repository&&main.object.sha===row.base&&commit.tree.sha===row.tree,"fresh exact merge guard "+row.label);
  store("ridge_premerge_guard_"+row.label,{pr,main,commit,legacy});
  const action=await tools.mcp__codex_apps__github_merge_pull_request({repository_full_name:row.repository,pr_number:row.pr,expected_head_sha:row.head,merge_method:"merge"});
  store("ridge_merge_action_"+row.label,action);
  need(!action.isError&&action.structuredContent?.merged===true&&/^[a-f0-9]{40}$/.test(action.structuredContent.sha),"confirmed merge receipt "+row.label+"; inspect canonical state before retry if uncertain");
  const actual=await get(prefix+"/git/commits/"+action.structuredContent.sha);
  need(actual.tree.sha===row.tree&&actual.parents.length===2&&actual.parents[0].sha===row.base&&actual.parents[1].sha===row.head,"actual ordered merge parents and exact tree");
  store("ridge_canonical_merge_"+row.label,actual);
  done.push({label:row.label,repository:row.repository,pr:row.pr,url:"https://github.com/"+row.repository+"/pull/"+row.pr,base:row.base,head:row.head,tree:row.tree,merge:actual.sha});
  store("ridge_actual_merges",done);
  text({status:"MERGED_EXACT_REVIEWED_TREE",...done.at(-1)});
 }
 store("ridge_both_merges_complete",true);
})();