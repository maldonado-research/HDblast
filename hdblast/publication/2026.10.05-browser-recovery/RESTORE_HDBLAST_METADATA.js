/* HDBLAST existing-draft metadata recovery. Run inside the signed-in Zenodo tab.
 * Supported protocol: Invenio's editor uses same-origin session credentials,
 * csrftoken / X-CSRFToken, application/json and ?expand=1 for a draft PUT.
 * The exact approved body is embedded and SHA256 checked before any request.
 */
(async () => {
  "use strict";
  const ORIGIN = "https://zenodo.org";
  const PAGE_PATH = "/uploads/23114217";
  const API = ORIGIN + "/api/records/23114217/draft?expand=1";
  const SELF = ORIGIN + "/api/records/23114217/draft";
  const BODY_SHA = "816a22d1ebfc3b0ef7b831764a37e305866c46485000b3a0cba5d8e4d7b95314";
  const BODY_TEXT = "{\n  \"metadata\": {\n    \"resource_type\": {\n      \"id\": \"software\"\n    },\n    \"creators\": [\n      {\n        \"person_or_org\": {\n          \"type\": \"personal\",\n          \"name\": \"Maldonado, Ricardo\",\n          \"given_name\": \"Ricardo\",\n          \"family_name\": \"Maldonado\",\n          \"identifiers\": [\n            {\n              \"identifier\": \"0009-0009-3937-6527\",\n              \"scheme\": \"orcid\"\n            }\n          ]\n        },\n        \"role\": {\n          \"id\": \"researcher\"\n        },\n        \"affiliations\": [\n          {\n            \"name\": \"Independent Researcher\"\n          }\n        ]\n      }\n    ],\n    \"title\": \"HDBLAST: Conditional Source/Operator Enclosures and Preserved Metric-Response Ledger Diagnostics\",\n    \"publisher\": \"Zenodo\",\n    \"publication_date\": \"2026-10-03\",\n    \"subjects\": [\n      {\n        \"subject\": \"HDBLAST\"\n      },\n      {\n        \"subject\": \"higher-dimensional-origin hypothesis\"\n      },\n      {\n        \"subject\": \"higher-dimensional cosmology\"\n      },\n      {\n        \"subject\": \"Einstein-scalar system\"\n      },\n      {\n        \"subject\": \"open-FRW cosmology\"\n      },\n      {\n        \"subject\": \"hyperbolic spatial geometry\"\n      },\n      {\n        \"subject\": \"kinetic-tail asymptotics\"\n      },\n      {\n        \"subject\": \"cosmological perturbation theory\"\n      },\n      {\n        \"subject\": \"transverse-traceless tensor perturbations\"\n      },\n      {\n        \"subject\": \"principal L2(H4) scalar perturbations\"\n      },\n      {\n        \"subject\": \"spectral wave packets\"\n      },\n      {\n        \"subject\": \"frequency-uniform energy estimate\"\n      },\n      {\n        \"subject\": \"quiescent singularity\"\n      },\n      {\n        \"subject\": \"logarithmic memory\"\n      },\n      {\n        \"subject\": \"singular scattering\"\n      },\n      {\n        \"subject\": \"Bessel asymptotics\"\n      },\n      {\n        \"subject\": \"Sobolev trace law\"\n      },\n      {\n        \"subject\": \"validated numerics\"\n      },\n      {\n        \"subject\": \"computer-assisted mathematics\"\n      },\n      {\n        \"subject\": \"negative results\"\n      },\n      {\n        \"subject\": \"reproducible computational physics\"\n      },\n      {\n        \"subject\": \"claim-safe release\"\n      },\n      {\n        \"subject\": \"characteristic initial value problem\"\n      },\n      {\n        \"subject\": \"Goursat problem\"\n      },\n      {\n        \"subject\": \"Riemann function\"\n      },\n      {\n        \"subject\": \"Ferrers functions\"\n      },\n      {\n        \"subject\": \"adjoint equations\"\n      },\n      {\n        \"subject\": \"compact constraint data\"\n      },\n      {\n        \"subject\": \"weighted contraction bounds\"\n      },\n      {\n        \"subject\": \"rank-two response matrix\"\n      },\n      {\n        \"subject\": \"scalar conjugation\"\n      },\n      {\n        \"subject\": \"exact arithmetic\"\n      },\n      {\n        \"subject\": \"registered numerical experiment\"\n      },\n      {\n        \"subject\": \"metric response\"\n      },\n      {\n        \"subject\": \"conservation diagnostic\"\n      },\n      {\n        \"subject\": \"oscillatory quadrature\"\n      },\n      {\n        \"subject\": \"failure evidence\"\n      },\n      {\n        \"subject\": \"Duhamel operators\"\n      },\n      {\n        \"subject\": \"interval arithmetic\"\n      },\n      {\n        \"subject\": \"Cauchy remainder bounds\"\n      },\n      {\n        \"subject\": \"source forcing\"\n      },\n      {\n        \"subject\": \"Arb ball arithmetic\"\n      },\n      {\n        \"subject\": \"exact rational defect\"\n      },\n      {\n        \"subject\": \"registered source/operator certificate\"\n      }\n    ],\n    \"languages\": [\n      {\n        \"id\": \"eng\"\n      }\n    ],\n    \"related_identifiers\": [\n      {\n        \"identifier\": \"10.5281/zenodo.17069900\",\n        \"scheme\": \"doi\",\n        \"relation_type\": {\n          \"id\": \"issupplementto\"\n        },\n        \"resource_type\": {\n          \"id\": \"dataset\"\n        }\n      },\n      {\n        \"identifier\": \"10.5281/zenodo.16937520\",\n        \"scheme\": \"doi\",\n        \"relation_type\": {\n          \"id\": \"issupplementto\"\n        },\n        \"resource_type\": {\n          \"id\": \"dataset\"\n        }\n      },\n      {\n        \"identifier\": \"10.5281/zenodo.16907982\",\n        \"scheme\": \"doi\",\n        \"relation_type\": {\n          \"id\": \"issupplementto\"\n        },\n        \"resource_type\": {\n          \"id\": \"dataset\"\n        }\n      },\n      {\n        \"identifier\": \"10.5281/zenodo.17088133\",\n        \"scheme\": \"doi\",\n        \"relation_type\": {\n          \"id\": \"issupplementto\"\n        },\n        \"resource_type\": {\n          \"id\": \"dataset\"\n        }\n      },\n      {\n        \"identifier\": \"10.5281/zenodo.17547897\",\n        \"scheme\": \"doi\",\n        \"relation_type\": {\n          \"id\": \"issupplementto\"\n        },\n        \"resource_type\": {\n          \"id\": \"dataset\"\n        }\n      },\n      {\n        \"identifier\": \"https://arxiv.org/abs/2401.08437v1\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"publication-preprint\"\n        }\n      },\n      {\n        \"identifier\": \"10.5281/zenodo.22347452\",\n        \"scheme\": \"doi\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"10.5281/zenodo.22922927\",\n        \"scheme\": \"doi\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"publication\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/commit/a8394d0127e58200cd4b9a87c14ae63a9dd69f02\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/pull/20\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/commit/06984aa6b142499c592850e6b4afd49c28ced394\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/pull/18\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/commit/79342726f097affa0c38d7726105d109f9b52a5c\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/commit/78f2ec797b5c3e4d1073f297b60ff946ad8bde71\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/commit/ead803c9ba575877fbda000cf2452f3ea35e272b\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/commit/5a851f22c84f14b454d398cebc6a4cbfa21d19da\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/blob/5a851f22c84f14b454d398cebc6a4cbfa21d19da/research/HDBLAST_CHECKPOINT_20261003_VERIFIED_SOURCE_OPERATOR/RESULTS.md\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      },\n      {\n        \"identifier\": \"https://github.com/maldonado-research/HDblast/blob/5a851f22c84f14b454d398cebc6a4cbfa21d19da/research/HDBLAST_VERIFIED_SOURCE_OPERATOR_PACKAGE_20261003/FRESH_REPRODUCTION.json\",\n        \"scheme\": \"url\",\n        \"relation_type\": {\n          \"id\": \"references\"\n        },\n        \"resource_type\": {\n          \"id\": \"software\"\n        }\n      }\n    ],\n    \"rights\": [\n      {\n        \"id\": \"cc-by-4.0\"\n      }\n    ],\n    \"copyright\": \"Copyright \u00a9 2026 Ricardo Maldonado. Creative Commons Attribution 4.0 International (CC BY 4.0)\",\n    \"description\": \"<p>This HDBLAST main-program methods checkpoint encloses the registered source and Duhamel operators in a specified fixed-geometry model. On a=-9/2, b=-7/2 and eta in [a,b], z=eta+4, B=exp(-z^2/(1-z^2)), h=B or h=zB, and L=-1/eta, the forcing is g=4L^2 h-2L h\u2032-h\u2032\u2032. It encloses integral g, integral exp(2ik(b-s))g(s), integral Phi_k(b-s)g(s), and the separate primitive integral Lg, with Phi_k(t)=(exp(2ikt)-1)/(2ik) and Phi_0(t)=t. Integral Lg is not the action-derived pressure/contact work.</p><p>Two independently implemented routes return 1,170 moment rows and 130 source-primitive rows each, across 64 panels, two source profiles, nine exact registered momenta and whole-interval reductions. All corresponding real and imaginary enclosures overlap. The conservative maximum whole-interval L1 hull radius is exactly 1459446311/10633823966279326983230456482242756608, strictly below 1.373e-28 and the unchanged 1e-26 absolute model-unit gate. The bound includes source approximation, phase or ODE defect, arithmetic, error propagation and outward endpoint serialization. Signed residuals or cross-method agreement are not substitutes for these bounds.</p><p>The primary route uses 256-bit Arb balls, degree 24 source jets and degree 96 entire local phase kernels. The independent route uses its own exact-rational source recurrences, degree 200 scalar-exponential enclosure, 512-bit chosen coefficients and degree 96 polynomial ODE defect. Each constructs the two real source profiles on 64 registered panels exactly once, records every source attempt before its callback and completes within the fixed 120-second/128-MiB resource caps. No retained quantum-state or numerical archive arrays are decoded.</p><p>The analytic source and local kernel bounds hold throughout the registered interior interval and real k in [0,256]. The complete computed numerical widths cover only nine exact probes: 0, 2^-40, 2^-12, 1/4, 1, 16, 64, 128 and 256. They do not certify a continuum momentum integral or include the nonanalytic support endpoints. The conditional conclusion depends on the archived analytic identities and proofs, numerical-library contracts and custodian public-freeze/readback chronology. All 155 prospective registered public files were verified before the first real-source run. A clean standalone ZIP extraction and optimized replay reproduces the complete numerical payload and stable scientific evidence exactly, with all frozen payload bytes unchanged.</p><p>The preceding registered metric-response experiments retain scientific FAIL, including 29 conservation endpoint and 30 refinement checks outside their unchanged gates. The source-free and active-source saved-data diagnostics retain LEDGER_ERROR_DEMONSTRATED: all twelve consistency cases pass and eight meet their unchanged empirical attribution criterion. Those diagnoses do not supply the new source/operator error certificate or repair every historical metric failure. This deposit preserves the original execution failure, disclosed pre-execution repair chronology, raw failure evidence and prior diagnostic packages.</p><p>All ten inherited files from semantic version 2026.09.05-v24 are retained unchanged as historical material. The earlier fixed-background five-field result remains a distinct conditional rank-two certificate under its stated model, background-root, Green/source, stability and gluing premises. This checkpoint neither withdraws that certificate nor extends it to coupled nonlinear cosmology. The scalar-junction/static-branch companion, record 23111008 in concept family 10.5281/zenodo.22922927, remains a separate unchanged study.</p><p>The intended complete main-family draft contains 23 files: ten inherited historical files, the twelve previously reviewed ledger/metric additions, and one new standalone registered-source/operator ZIP. The new manuscript is inside that ZIP; its external fresh-reproduction receipt is linked through the immutable public GitHub source and remains outside the archive it verifies. Reports/data use CC BY 4.0; bundled code retains its included MIT license and third-party material retains its included licenses. The established main concept DOI remains 10.5281/zenodo.17088132.</p><p>The full twelve-case action-derived pressure/contact conservation certificate remains unresolved. Incoming state accuracy, momentum integration, renormalization, geometry response, coupled backreaction, nonlinear persistence, stability, physical source-energy normalization, heating and thermalization, a hot radiation era, agreement with cosmological data and a higher-dimensional cause of the Big Bang remain not established. No new law of nature is claimed; external mathematical novelty is not assessed. AI-assisted internal mathematical and computational reviews are not external peer review or proof-assistant formalization.</p>\",\n    \"additional_descriptions\": [\n      {\n        \"description\": \"<p>Semantic version 2026.10.03-verified-source-operator is dated October 3, 2026, Pacific time; it is separate from the display ordinal Zenodo assigns. Start with the new ZIP\u2019s RESULTS.md and the immutable GitHub research manuscript: https://github.com/maldonado-research/HDblast/blob/5a851f22c84f14b454d398cebc6a4cbfa21d19da/research/HDBLAST_CHECKPOINT_20261003_VERIFIED_SOURCE_OPERATOR/RESULTS.md. The prospective source freeze is commit ead803c9ba575877fbda000cf2452f3ea35e272b; the complete scientific result and reviewed package are pinned to commit 5a851f22c84f14b454d398cebc6a4cbfa21d19da.</p><p>The registered source/operator ZIP is 5,331,624 bytes, SHA256 dfc8c6c6144fa1d9039390bcbaa7a9bd062355863c0fe1494d63824721193cd7, with 175 files including MANIFEST.json. Its whole-payload manifest SHA256 is d870a76e15fe49da9f624a55258f3f6aa93933ea19860b9bff7c4f2da3ec93f6. The external fresh-reproduction receipt is https://github.com/maldonado-research/HDblast/blob/5a851f22c84f14b454d398cebc6a4cbfa21d19da/research/HDBLAST_VERIFIED_SOURCE_OPERATOR_PACKAGE_20261003/FRESH_REPRODUCTION.json, SHA256 672300640f4256597eede8ccf10f625903deb3ddbba316b14813c48cee95602c. Both archived and fresh scientific states are PASS_UNIFORM_MODEL_AND_REGISTERED_PROBE_CERTIFICATE; the fresh replay status is PASS_EXACT_FRESH_REPLAY. These labels do not establish the unresolved pressure/contact ledger or physical Big Bang mechanism.</p><p>All twelve preceding additions retain their exact filenames, byte lengths, MD5 and SHA256. Older V25 and date-stamped archive/report filenames are historical preparation labels. The source-free ZIP remains 48,470,467 bytes, SHA256 b98291b7d8a2a40d5ddf79d877bb0177c03b748becc726d63a94456adc4830a0; the active-source ZIP remains 32,195,260 bytes, SHA256 0cd3dd2da604444541c06e98167f537f769c0ce430108bf02354e5090f328c4f. Both retain their separately delivered replay receipts and registered LEDGER_ERROR_DEMONSTRATED classification.</p><p>The inherited ten files describe historical v24 record 22347452, semantic version 2026.09.05-v24. Record 23112891 was an incomplete inherited-only test upload removed by its owner and now returns a removal tombstone. Preserve concept 10.5281/zenodo.17088132 and reuse existing owned draft 23114217. Do not create a second standalone main record or another draft. The separate companion record 23111008, concept 22922927, retains its original two files.</p><p>The existing draft was freshly authenticated as unsubmitted with ten unchanged inherited files. Earlier metadata save returned HTTP 500 and binary upload returned HTTP 400/Bad Content-Length; saved metadata and file readback did not verify a completed update. No write retry, upload, release, deletion, new version or publication was performed while preparing this revised candidate. Automatic GitHub archiving remains OFF for all eight public repositories. This prepared browser/API packet is not proof that the new version has been published. Publication requires exact saved metadata, all 23 filename/size/checksum pins, the intended family and direct public HTTP200 confirmation. Sealed archive publication-status text records its preparation snapshot, not current Zenodo state.</p>\",\n        \"type\": {\n          \"id\": \"notes\"\n        },\n        \"lang\": {\n          \"id\": \"eng\"\n        }\n      },\n      {\n        \"description\": \"<p>Historical v24 summary: the following text describes the unchanged inherited five-field linear-response checkpoint, published as record 22347452 (semantic version 2026.09.05-v24). It is retained for historical context. The main Description above describes the current source/operator update.</p><p>Historical v24 plain-language summary</p>\\n<p>This checkpoint completes an error-controlled calculation across the whole domain of a registered fixed-background five-field linear mathematical model. It joins 176 calculation regions with exactly matching boundaries, then includes bounds for the equation error, background uncertainty, boundary-data error and source integration. The resulting two-by-two response matrix has a determinant between 2.12278154373842557755059 and 2.53743323308364125340774. Because this interval is entirely positive, the two selected input pulses produce linearly independent output patterns in the two chosen mathematical readouts, within the model's stated assumptions.</p>\\n<p>The advance is a complete certification argument for this adjoint calculation over the full domain. Earlier project checkpoints already reported rank-two results by other arguments. This release does not claim a more precise determinant estimate or external mathematical novelty.</p>\\n<p>The conclusion depends on the archived background branch and inherited equations, source identities and stability bounds. A separate rigorous error enclosure for the independent forward calculation remains open. The detector normalizations are mathematical coefficient conventions; they do not establish physical measurement sensitivity or observational accessibility.</p>\\n<p>These results do not show that a higher-dimensional blast caused the Big Bang. Nonlinear persistence, agreement with cosmological observations and the proposed physical mechanism remain unestablished. The reviews and replays are internal checks, not external peer review or proof-assistant formalization.</p>\",\n        \"type\": {\n          \"id\": \"other\"\n        },\n        \"lang\": {\n          \"id\": \"eng\"\n        }\n      }\n    ],\n    \"references\": [\n      {\n        \"reference\": \"Randall & Sundrum (1999), RS2\"\n      },\n      {\n        \"reference\": \"Shiromizu, Maeda & Sasaki (2000), effective brane equations\"\n      },\n      {\n        \"reference\": \"Maartens review on brane-world gravity\"\n      },\n      {\n        \"reference\": \"NIST Digital Library of Mathematical Functions, \u00a715.5, equation 15.5.4 (Gauss hypergeometric differentiation identity). https://dlmf.nist.gov/15.5.E4\"\n      },\n      {\n        \"reference\": \"Forrester, P. J. (2018), Meet Andr\u00e9ief, Bordeaux 1886, and Andreev, Kharkov 1882\u201383. arXiv:1806.10411. https://arxiv.org/abs/1806.10411\"\n      },\n      {\n        \"reference\": \"Johansson, F. (2017). Arb: Efficient Arbitrary-Precision Midpoint-Radius Interval Arithmetic. IEEE Transactions on Computers. https://doi.org/10.1109/TC.2017.2690633\"\n      },\n      {\n        \"reference\": \"Johansson, F. (2018). Numerical Integration in Arbitrary-Precision Ball Arithmetic. https://doi.org/10.1007/978-3-319-96418-8_30\"\n      }\n    ],\n    \"version\": \"2026.10.03-verified-source-operator\"\n  },\n  \"access\": {\n    \"record\": \"public\",\n    \"files\": \"public\",\n    \"embargo\": {\n      \"active\": false,\n      \"reason\": null\n    }\n  },\n  \"custom_fields\": {\n    \"code:programmingLanguage\": [\n      {\n        \"id\": \"python\"\n      },\n      {\n        \"id\": \"shell\"\n      },\n      {\n        \"id\": \"c++\"\n      }\n    ],\n    \"code:developmentStatus\": {\n      \"id\": \"active\"\n    }\n  },\n  \"files\": {\n    \"enabled\": true,\n    \"default_preview\": \"00_READ_FIRST.md\",\n    \"order\": []\n  }\n}\n";
  const INHERITED = [{"filename":"HDBLAST_CONDITIONAL_GLOBAL_ADJOINT_RESPONSE_CERTIFICATION_20260904.zip","bytes":15455879,"md5":"c63b6c68e5e5c0dea5f21db36f6b09b1"},{"filename":"HDBLAST_CONDITIONAL_GLOBAL_ADJOINT_RESPONSE_CERTIFICATION_20260904.sha256","bytes":435,"md5":"972b0f8c3bce439201c3dfd88c7c1628"},{"filename":"ZENODO_RELEASE_RECOMMENDATION_AND_NOTES.md","bytes":4953,"md5":"4397859f8c7d9b332eb2f90999bfec1d"},{"filename":"PUBLICATION_README.md","bytes":1675,"md5":"0d417c564aafc7d8ec00fe57c0b07687"},{"filename":"DELIVERY_SHA256SUMS.txt","bytes":1053,"md5":"12b33265994b7cf88b4f73f7741bd977"},{"filename":"METHODS_AND_LITERATURE_REVIEW.md","bytes":4838,"md5":"eed0ca61772d65bd7694ba03ab07a67d"},{"filename":"00_READ_FIRST.md","bytes":7875,"md5":"fa13c49aae6dafb81e336075e3645235"},{"filename":"HDBLAST_CONDITIONAL_GLOBAL_ADJOINT_RESPONSE_CERTIFICATION_20260904_REPLAY.json","bytes":4461,"md5":"5508d0d2908699497628fda36fce02f9"},{"filename":"GLOBAL_RESPONSE_CERTIFICATE.json","bytes":71963,"md5":"ca71277a5809614b9649504c0b25e0f0"},{"filename":"HDBLAST_CONDITIONAL_GLOBAL_ADJOINT_RESPONSE_CERTIFICATION_20260904_DELIVERY_AUDIT.json","bytes":1512,"md5":"8575d1e8d2c3d77fe48876dc498663d8"}];
  const STORAGE_KEY = "HDBLAST_METADATA_REPAIR_23114217_" + BODY_SHA;
  const LOCK = "HDBLAST_METADATA_REPAIR_IN_FLIGHT";
  const result = { record_id: 23114217, status: "STARTED", publication_ready: false,
    uploads: 0, publications: 0, deletions: 0, metadata_write_attempts: 0 };
  if (window[LOCK]) {
    console.info("HDBLAST: a metadata recovery is already running.");
    return;
  }
  window[LOCK] = true;
  let writeStarted = false;
  let reloadAfterSuccess = false;
  function fail(code) {
    const error = new Error(code);
    error.hdblastCode = code;
    throw error;
  }
  function canonical(value) {
    if (Array.isArray(value)) return JSON.stringify(value.map(x => JSON.parse(canonical(x))));
    if (value && typeof value === "object") {
      const sorted = {};
      for (const key of Object.keys(value).sort()) sorted[key] = JSON.parse(canonical(value[key]));
      return JSON.stringify(sorted);
    }
    return JSON.stringify(value);
  }
  function empty(value) {
    return value == null || value === "" ||
      (typeof value === "object" && Object.keys(value).length === 0);
  }
  // Server vocabulary expansion may add labels/links. Expected leaves must match.
  function matches(actual, expected) {
    if (Array.isArray(expected)) {
      if (!Array.isArray(actual) || actual.length !== expected.length) return false;
      const unused = actual.slice();
      for (const item of expected) {
        const at = unused.findIndex(value => matches(value, item));
        if (at === -1) return false;
        unused.splice(at, 1);
      }
      return true;
    }
    if (expected && typeof expected === "object") {
      return actual && typeof actual === "object" && !Array.isArray(actual) &&
        Object.keys(expected).every(key => Object.hasOwn(actual, key) && matches(actual[key], expected[key]));
    }
    return actual === expected;
  }
  function approvedPartial(actual, expected) {
    if (empty(actual)) return true;
    if (Array.isArray(actual)) {
      if (!Array.isArray(expected) || actual.length > expected.length) return false;
      const unused = expected.slice();
      for (const item of actual) {
        const at = unused.findIndex(value => approvedPartial(item, value));
        if (at === -1) return false;
        unused.splice(at, 1);
      }
      return true;
    }
    if (actual && typeof actual === "object") {
      if (!expected || typeof expected !== "object" || Array.isArray(expected)) return false;
      const expanded = new Set(["title", "props", "links", "icon"]);
      return Object.keys(actual).every(key => Object.hasOwn(expected, key)
        ? approvedPartial(actual[key], expected[key])
        : (Object.hasOwn(expected, "id") && expanded.has(key)) || empty(actual[key]));
    }
    return actual === expected;
  }
  function checkRecord(data) {
    if (String(data?.id) !== "23114217" || String(data?.parent?.id) !== "17088132" ||
        String(data?.parent?.access?.owned_by?.user) !== "1386319") fail("WRONG_DRAFT_OWNER_OR_DOI_FAMILY");
    if (data.is_draft !== true || data.is_published !== false) fail("DRAFT_STATE_IS_NOT_UNPUBLISHED");
    if (data.links?.self !== SELF) fail("UNEXPECTED_DRAFT_SELF_LINK");
    if (data.access?.record !== "public" || data.access?.files !== "public" ||
        data.access?.embargo?.active === true) fail("VISIBILITY_REQUIRES_REVIEW");
    const entries = data.files?.entries;
    if (data.files?.count !== 10 || data.files?.total_bytes !== 15554644 ||
        !entries || Object.keys(entries).length !== INHERITED.length) fail("INHERITED_FILE_INVENTORY_CHANGED");
    for (const item of INHERITED) {
      const actual = entries[item.filename];
      if (!actual || actual.size !== item.bytes || actual.checksum !== "md5:" + item.md5)
        fail("INHERITED_FILE_BYTES_OR_CHECKSUM_CHANGED");
    }
    if (!Number.isSafeInteger(data.revision_id) || data.revision_id < 0) fail("INVALID_DRAFT_REVISION");
  }
  async function request(method, headers, body) {
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), 45000);
    try {
      const response = await fetch(API, { method, headers, body, credentials: "same-origin",
        mode: "same-origin", redirect: "error", cache: "no-store", signal: controller.signal });
      if (response.url !== API) fail("UNEXPECTED_RESPONSE_URL");
      const text = await response.text();
      if (text.length > 4 * 1024 * 1024) fail("RESPONSE_TOO_LARGE");
      let data;
      try { data = JSON.parse(text); } catch { fail("RESPONSE_WAS_NOT_JSON"); }
      return { status: response.status, data };
    } finally {
      clearTimeout(timer);
    }
  }
  async function read(headers) {
    const response = await request("GET", headers);
    if (response.status !== 200) fail("READ_HTTP_" + response.status);
    checkRecord(response.data);
    return response.data;
  }
  function mismatches(data, wanted) {
    const values = [];
    for (const key of Object.keys(wanted.metadata)) {
      if (!matches(data.metadata?.[key], wanted.metadata[key])) values.push("metadata." + key);
    }
    for (const key of Object.keys(wanted.custom_fields)) {
      if (!matches(data.custom_fields?.[key], wanted.custom_fields[key])) values.push("custom_fields." + key);
    }
    if (!matches(data.access, wanted.access)) values.push("access");
    for (const key of Object.keys(wanted.files)) {
      if (!matches(data.files?.[key], wanted.files[key])) values.push("files." + key);
    }
    return values;
  }
  try {
    if (window.location.origin !== ORIGIN || window.location.pathname.replace(/\/$/, "") !== PAGE_PATH)
      fail("OPEN_THE_EXISTING_ZENODO_DRAFT_23114217_FIRST");
    if (!crypto?.subtle || typeof TextEncoder !== "function") fail("BROWSER_CRYPTO_UNAVAILABLE");
    const digest = await crypto.subtle.digest("SHA-256", new TextEncoder().encode(BODY_TEXT));
    const hex = Array.from(new Uint8Array(digest), x => x.toString(16).padStart(2, "0")).join("");
    if (hex !== BODY_SHA) fail("APPROVED_METADATA_BODY_HASH_MISMATCH");
    const wanted = JSON.parse(BODY_TEXT);
    const headers = { Accept: "application/vnd.inveniordm.v1+json", "Content-Type": "application/json" };
    const cookie = document.cookie.split(";").map(value => value.trim()).find(value => value.startsWith("csrftoken="));
    if (cookie) {
      try { headers["X-CSRFToken"] = decodeURIComponent(cookie.slice("csrftoken=".length)); }
      catch { fail("BROWSER_CSRF_COOKIE_FORMAT_REQUIRES_REVIEW"); }
    }
    const before = await read(headers);
    result.revision_before = before.revision_id;
    if (mismatches(before, wanted).length === 0) {
      result.status = "PASS_ALREADY_SAVED_METADATA_READBACK";
      result.revision_after = before.revision_id;
      reloadAfterSuccess = true;
    } else {
      if (sessionStorage.getItem(STORAGE_KEY)) fail("EARLIER_BROWSER_WRITE_RECORDED_RECONCILE_BEFORE_RETRY");
      for (const group of ["metadata", "custom_fields"]) {
        for (const key of Object.keys(before[group] || {})) {
          if (!empty(before[group][key]) &&
              (!Object.hasOwn(wanted[group], key) || !approvedPartial(before[group][key], wanted[group][key])))
            fail("UNREVIEWED_EXISTING_" + group.toUpperCase() + "_REQUIRES_REVIEW");
        }
      }
      sessionStorage.setItem(STORAGE_KEY, JSON.stringify({ phase: "WRITE_STARTED", revision: before.revision_id, body_sha256: BODY_SHA }));
      writeStarted = true;
      result.metadata_write_attempts = 1;
      const saved = await request("PUT", { ...headers, "If-Match": String(before.revision_id) }, BODY_TEXT);
      result.write_http_status = saved.status;
      const after = await read(headers);
      result.revision_after = after.revision_id;
      if (canonical(before.pids || {}) !== canonical(after.pids || {})) fail("RESERVED_IDENTIFIERS_CHANGED_REQUIRES_REVIEW");
      result.mismatch_fields = mismatches(after, wanted);
      if (saved.status !== 200 || result.mismatch_fields.length) fail("WRITE_READBACK_DID_NOT_VERIFY_ALL_APPROVED_METADATA");
      result.status = "PASS_BROWSER_SAVED_METADATA_READBACK";
      sessionStorage.setItem(STORAGE_KEY, JSON.stringify({ phase: "VERIFIED_METADATA", revision: after.revision_id, body_sha256: BODY_SHA }));
      reloadAfterSuccess = true;
    }
    result.checked_inherited_files = 10;
    result.body_sha256 = BODY_SHA;
    window.HDBLAST_METADATA_REPAIR_RESULT = result;
    console.info("HDBLAST metadata recovery result:", result);
    alert("HDBLAST metadata is saved and its browser readback matches the approved fields.\n\nClick OK to refresh this editor, then return to Codex for an independent saved-record check. The thirteen additions still need uploading after that check.");
  } catch (error) {
    result.status = writeStarted ? "BLOCKED_RECONCILE_BROWSER_WRITE_BEFORE_RETRY" : "BLOCKED_NO_METADATA_WRITE";
    result.reason = error?.hdblastCode || "BROWSER_REQUEST_OR_CAPABILITY_ERROR";
    window.HDBLAST_METADATA_REPAIR_RESULT = result;
    console.info("HDBLAST metadata recovery result:", result);
    const next = writeStarted
      ? "A save may have applied. Close this editor tab without clicking Save draft or Publish, and return to Codex with this message. Reopen it after the saved state is reconciled."
      : "Return to Codex with this message. Leave the draft unpublished; a recorded write must be reconciled before another attempt.";
    alert("HDBLAST recovery stopped.\n\n" + result.reason + "\n\n" + next);
  } finally {
    window[LOCK] = false;
  }
  if (reloadAfterSuccess) window.location.reload();
})()
