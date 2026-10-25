#!/usr/bin/env node
"use strict";
// Deterministically export independent exact polygon/endpoint fixtures.
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const source = fs.readFileSync(path.join(__dirname, "independent_kernel_checks.js"), "utf8");
const marker = 'return{status:"PASS_EXACT_REFLECTED_RIDGE_KERNEL_CONTROLS"';
if (source.split(marker).length !== 2) throw Error("unique fixture seam required");
const inject = "(() => {\n const pick=[0,8,17,29,37,46,57,58];\n const encode=P=>({polynomial:P.p.map(str),hinges:P.hinges.map(row=>row.map(str))});\n const cases=[];\n for(let si=0;si<sources.length;si++)for(const J of[-1,1])for(const i of pick){\n  const[x,r]=points[i], [F,G]=sources[si], value=directK(F,G,x,r,J), d=endpointDerivs(F,G,x,r,J);\n  cases.push({source:si,reflection:J,x:str(x),r:str(r),value:str(value),dx:str(d[0]),dr:str(d[1])});\n }\n return {schema:\"HDBLAST_REFLECTED_RIDGE_EXACT_FIXTURES_V1\",provenance:\"Exact polygon integration and endpoint integrals in independent_kernel_checks.js; no numerical quadrature.\",sources:sources.map(([F,G])=>({F:encode(F),G:encode(G)})),cases};\n})()";
let output = "";
vm.runInNewContext(source.replace(marker, 'return{fixtures:' + inject + ',status:"PASS_EXACT_REFLECTED_RIDGE_KERNEL_CONTROLS"'),
  {process: {stdout: {write: value => { output += value; }}}}, {timeout: 30000});
process.stdout.write(JSON.stringify(JSON.parse(output).fixtures, null, 2) + "\n");
