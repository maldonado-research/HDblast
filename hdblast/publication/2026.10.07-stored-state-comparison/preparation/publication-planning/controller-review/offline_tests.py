#!/usr/bin/env python3
"""Offline stateful Zenodo fixtures. Uses fake bytes/tokens; forbids sockets.

The injected synthetic baseline has 27 small manufactured files.
All 29 edition files undergo fresh fake content streams in both states.
These are test data; their total is not the live 448387915-byte inventory.
These fixtures test control flow, not real Zenodo byte contents or availability.
"""
from __future__ import annotations
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import socket
import sys
import tempfile
import unittest
from urllib.parse import urlsplit

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / "new_edition_controller.py"
spec = importlib.util.spec_from_file_location("edition_under_test", SOURCE)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)
M = mod
DRAFT = "24000000"
FAKE_TOKEN = "OFFLINE_FAKE_TOKEN_ONLY_71d2"

def clone(x):
    return copy.deepcopy(x)

def pin(name, data):
    return {"filename": name, "bytes": len(data), "md5": hashlib.md5(data).hexdigest(),
            "sha256": hashlib.sha256(data).hexdigest()}

def response(obj=None, status=200, headers=None, raw=None):
    return M.Response(status, headers or {}, raw if raw is not None else json.dumps(obj).encode())

class FakeServer:
    """In-memory stateful routes, with deliberately ambiguous critical POSTs."""
    def __init__(self, old, new, metadata, payloads):
        self.old, self.new, self.metadata = clone(old), clone(new), clone(metadata)
        self.payloads = clone(payloads)
        self.prior_metadata = clone(metadata)
        self.prior_metadata["metadata"]["title"] = "Prior HDBLAST edition"
        self.prior_metadata["metadata"]["version"] = M.PRIOR_SEMANTIC_VERSION
        self.prior_metadata["metadata"]["description"] = "<p>Prior description.</p>"
        self.exists = False
        self.published = False
        self.draft_metadata = clone(self.prior_metadata)
        self.files = {}
        self.calls, self.streams = [], []
        self.create_form = "original"
        self.create_fault = self.publish_fault = None
        self.create_reject_status = None
        self.metadata_error_response = None
        self.import_drop = self.wrong_family = False
        self.latest_redirect = True
        self.owner_pages = None
        self.owner_latest_only = False
        self.draft_etag = '"7"'
        self.draft_revision = 7
        self.old_files = {p["filename"]: self.entry(p, old=True) for p in old}

    def entry(self, p, old=False):
        name = p["filename"]
        prefix = f"{M.ORIGIN}/api/records/{M.PRIOR_ID if old else DRAFT}"
        if not old:
            prefix += "/draft"
        prefix += "/files/" + name
        return {"key": name, "size": p["bytes"], "checksum": "md5:" + p["md5"],
                "status": "completed", "file_id": "file-" + name,
                "version_id": "version-" + name,
                "links": {"content": prefix + "/content", "commit": prefix + "/commit"}}

    def document(self, prior=False, public=False):
        ident = M.PRIOR_ID if prior else DRAFT
        base = f"{M.ORIGIN}/api/records/{ident}" + ("" if prior or public else "/draft")
        files = clone(self.old_files if prior else self.files)
        if public and not prior:
            for name, f in files.items():
                f["links"]["content"] = f"{M.ORIGIN}/api/records/{DRAFT}/files/{name}/content"
        obj = clone(self.prior_metadata if prior else self.draft_metadata)
        obj.update({"id": ident, "parent": {"id": "999" if self.wrong_family and not prior else M.PARENT_ID, "access":{"owned_by":{"user":str(M.EXPECTED_OWNER)}}},
                    "revision_id": self.draft_revision,
                    "is_published": prior or public, "is_draft": not (prior or public),
                    "files": {"entries": files},
                    "links": {"self": base, "files": base + "/files",
                              "latest": f"{M.ORIGIN}/api/records/{M.PRIOR_ID}/versions/latest",
                              "publish": base + "/actions/publish"}})
        return obj

    def owned(self):
        out = [{"id": int(M.PRIOR_ID), "conceptrecid": M.PARENT_ID,
                "owner": M.EXPECTED_OWNER, "submitted": True}]
        if self.exists:
            out.append({"id": int(DRAFT), "conceptrecid": M.PARENT_ID,
                        "owner": M.EXPECTED_OWNER, "submitted": self.published})
            if self.owner_latest_only:
                out = out[-1:]
        return out

    def critical_fault(self, kind):
        fault = getattr(self, kind + "_fault")
        setattr(self, kind + "_fault", None)
        if fault == "unknown_after":
            raise OSError("simulated socket loss; offline fake only")
        if fault == "malformed_after":
            return response(raw=b"{", status=201)
        return None

    def request(self, method, url, data=None, headers=None):
        M.safe_url(url, query=method == "GET")
        self.calls.append((method, url, clone(headers or {})))
        path = urlsplit(url).path
        priorbase = "/api/records/" + M.PRIOR_ID
        draftbase = "/api/records/" + DRAFT + "/draft"
        publicbase = "/api/records/" + DRAFT
        if method == "GET":
            if path == "/api/deposit/depositions":
                if self.owner_pages is not None:
                    from urllib.parse import parse_qs
                    page = int(parse_qs(urlsplit(url).query)["page"][0])
                    return response(self.owner_pages.get(page, []))
                return response(self.owned())
            if path == "/api/deposit/depositions/" + M.PRIOR_ID:
                return response({"id": int(M.PRIOR_ID), "conceptrecid": M.PARENT_ID,
                                 "submitted": True, "links": {"newversion": M.ORIGIN + path + "/actions/newversion"}})
            if path == priorbase:
                return response(self.document(prior=True), headers={"ETag": '"prior"'})
            if path == priorbase + "/versions/latest":
                target = publicbase if self.published else priorbase
                if self.latest_redirect:
                    return response({}, status=301, headers={"Location": M.ORIGIN + target})
                return response(self.document(public=True) if self.published else self.document(prior=True))
            if path == priorbase + "/files":
                return response({"entries": clone(self.old_files)})
            if path == draftbase and self.exists and not self.published:
                return response(self.document(), headers={"ETag": self.draft_etag})
            if path == draftbase + "/files" and self.exists and not self.published:
                return response({"entries": clone(self.files)})
            if path == publicbase and self.published:
                return response(self.document(public=True), headers={"ETag": '"public-1"'})
            if path == publicbase + "/files" and self.published:
                return response(self.document(public=True)["files"])
            return response({}, status=404)
        if method == "POST" and path == "/api/deposit/depositions/" + M.PRIOR_ID + "/actions/newversion":
            if (headers or {}).get("Accept") != "application/json":
                return response({}, status=406)
            if self.create_reject_status is not None:
                return response({}, status=self.create_reject_status)
            if self.create_fault == "reject_before":
                self.create_fault = None
                return response({}, status=409)
            if self.exists:
                raise AssertionError("second newversion POST")
            self.exists = True
            self.files = clone(self.old_files)
            for name, f in self.files.items():
                f["links"]["content"] = f"{M.ORIGIN}{draftbase}/files/{name}/content"
                f["links"]["commit"] = f"{M.ORIGIN}{draftbase}/files/{name}/commit"
            if self.import_drop:
                self.files.pop(self.old[0]["filename"])
            failure = self.critical_fault("create")
            if failure:
                return failure
            if self.create_form == "original":
                return response({"id": int(M.PRIOR_ID), "links": {"latest_draft": f"{M.ORIGIN}/api/deposit/depositions/{DRAFT}"}}, status=201)
            return response({"id": int(DRAFT)}, status=201)
        if method == "POST" and path == draftbase + "/files":
            if (headers or {}).get("Accept") != "application/json":
                return response({}, status=406)
            for row in json.loads(data):
                name = row["key"]
                p = next(p for p in self.new if p["filename"] == name)
                self.files[name] = self.entry(p)
                self.files[name].update({"size": 0, "checksum": "", "status": "pending"})
            return response({"entries": clone(self.files)}, status=201)
        if path.startswith(draftbase + "/files/"):
            if method in ("PUT", "POST") and (headers or {}).get("Accept") != "application/json":
                return response({}, status=406)
            name, action = path[len(draftbase + "/files/"):].rsplit("/", 1)
            if method == "PUT" and action == "content":
                self.payloads[name] = bytes(data)
                return response({}, status=200)
            if method == "POST" and action == "commit":
                b = self.payloads[name]
                self.files[name].update({"size": len(b), "checksum": "md5:" + hashlib.md5(b).hexdigest(), "status": "completed"})
                return response(self.files[name], status=200)
        if method == "PUT" and path == draftbase:
            if (headers or {}).get("Accept") != M.VENDOR:
                return response({}, status=406)
            if (headers or {}).get("If-Match") != str(self.draft_revision):
                return response({"errors":[{"field":"if_match","messages":["Not a valid integer."]}]},status=400)
            if self.metadata_error_response is not None:
                return self.metadata_error_response
            self.draft_metadata = json.loads(data)
            return response(self.document(), status=200)
        if method == "POST" and path == draftbase + "/actions/publish":
            if (headers or {}).get("Accept") != M.VENDOR:
                return response({}, status=406)
            if (headers or {}).get("If-Match") != str(self.draft_revision):
                raise AssertionError("publication header is not the expected decimal revision")
            if self.publish_fault == "reject_before":
                self.publish_fault = None
                return response({}, status=409)
            if self.published:
                raise AssertionError("second publish POST")
            self.published = True
            failure = self.critical_fault("publish")
            return failure or response(self.document(public=True), status=202)
        raise AssertionError(f"unexpected offline route {method} {path}")

    def stream(self, url, headers=None):
        M.safe_url(url)
        self.streams.append((url, clone(headers or {})))
        name = urlsplit(url).path.rsplit("/", 2)[-2]
        if name not in self.payloads:
            raise AssertionError("attempt to allocate/stream synthetic large identity-only row")
        data = self.payloads[name]
        for offset in range(0, len(data), 7):
            yield data[offset:offset + 7]

class ControllerFixtures(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="edition-offline-")
        self.root = Path(self.tmp.name)
        self.payloads = {f"old-{i:02d}.txt": f"old fixture {i}\n".encode() for i in range(27)}
        old = [pin(name, data) for name, data in self.payloads.items()]
        baseline = {"files": old, "concept_id": 17088132, "current_record": 23225288, "expected_bytes":sum(x["bytes"] for x in old)}
        prior = {"rows": [{"filename": p["filename"], "status": "PASS_PRIOR_COMPLETE_BYTES_WITH_FRESH_IMMUTABLE_IDENTITY",
                           "observed": {k: p[k] for k in ("bytes", "md5", "sha256")},
                           "remote_file_id": "file-" + p["filename"], "remote_version_id": "version-" + p["filename"]} for p in old]}
        new = []
        for i in range(2):
            name, b = f"new-{i}.txt", f"New frozen fixture addition {i}.\n".encode()
            asset = self.root / name
            asset.write_bytes(b)
            self.payloads[name] = b
            new.append(dict(pin(name, b), local_path=str(asset)))
        self.metadata = {"metadata": {"title": "HDBLAST final fixture edition", "version": "2026.10.07-offline-fixture",
                                      "description": "<p>Final fixture description.</p>", "creators": [{"person_or_org": {"name": "Fixture Author", "type": "personal"}}],
                                      "rights": [{"id": "cc-by-4.0"}], "resource_type": {"id": "software"},
                                      "publisher": "Zenodo", "copyright": "Fixture", "languages": [{"id": "eng"}], "publication_date":"2026-10-07", "additional_descriptions":[{"description":"<p>New notes.</p>","lang":{"id":"eng"},"type":{"id":"notes"}},{"description":"<p>Preserved old notes.</p>","type":{"id":"notes"}}],"related_identifiers":[{"identifier":"10.5281/zenodo.23114217","scheme":"doi","relation_type":{"id":"references"},"resource_type":{"id":"software"}}]},
                         "custom_fields": {}, "access": {"record": "public", "files": "public", "embargo": {"active": False}}}
        self.inventory = {"sealed_upload_manifest": True, "concept_doi": "10.5281/zenodo.17088132",
                          "inherited_files": old, "new_files": new, "total_bytes_final": sum(p["bytes"] for p in old + new),
                          "expected_total_files_final": 29, "scientific_replay_status": M.SCIENTIFIC_REPLAY_STATUS}
        self.baseline, self.prior = baseline, prior
        self.server = FakeServer(old, new, self.metadata, self.payloads)
        self.controller = self.make_controller()
        self.original_socket = socket.socket
        socket.socket = lambda *a, **k: (_ for _ in ()).throw(AssertionError("real sockets forbidden in offline fixtures"))

    def tearDown(self):
        socket.socket = self.original_socket
        self.tmp.cleanup()

    def make_controller(self):
        return M.Controller(self.root / "evidence", self.server, token=FAKE_TOKEN,
                            inventory=clone(self.inventory), metadata=clone(self.metadata),
                            inventory_sha=M.sha(M.canonical_bytes(self.inventory)), metadata_sha=M.sha(M.canonical_bytes(self.metadata)),
                            baseline=clone(self.baseline), prior_receipt=clone(self.prior),
                            published_baseline=clone(self.server.prior_metadata))

    def posts(self, suffix):
        return [c for c in self.server.calls if c[0] == "POST" and c[1].endswith(suffix)]

    def stop(self, code, fn):
        with self.assertRaises(M.Stop) as ctx:
            fn()
        self.assertEqual(str(ctx.exception), code)

    def prepared(self):
        self.controller.create()
        result = self.controller.prepare()
        self.assertEqual(result["file_count"], 29)
        self.assertFalse(result["published"])
        self.assertEqual(len(self.posts("/commit")), 2)
        return result

    def test_both_newversion_response_forms_and_latest_resolution(self):
        for form in ("original", "newdraft"):
            with self.subTest(form=form):
                if form == "newdraft":
                    self.server = FakeServer(self.baseline["files"], self.inventory["new_files"], self.metadata, self.payloads)
                    (self.root / "evidence").rename(self.root / "evidence-preserved-first-form")
                    self.controller = self.make_controller()
                self.server.create_form = form
                result = self.controller.create()
                self.assertEqual(result["draft_id"], DRAFT)
                self.assertEqual(len(self.posts("/actions/newversion")), 1)
                self.assertEqual(len(self.server.files), 27)
                self.assertTrue(all(headers.get("_public") for method, url, headers in self.server.calls if "/versions/latest" in url))

    def test_wrong_family_blocks_adoption(self):
        self.server.wrong_family = True
        self.stop("DRAFT_FAMILY_OR_STATE_DIFFERS", self.controller.create)
        self.assertEqual(len(self.posts("/actions/newversion")), 1)
        self.assertIsNone(self.controller.state["draft_id"])

    def test_latest_only_owner_listing_complete_create_prepare_publish_flow(self):
        self.server.owner_latest_only = True
        self.prepared()
        self.assertEqual(len(self.posts("/actions/newversion")), 1)
        self.assertEqual(self.server.owned()[0]["id"], int(DRAFT))
        self.assertFalse(self.server.owned()[0]["submitted"])
        result = self.controller.publish()
        self.assertTrue(result["published"])
        self.assertTrue(self.server.owned()[0]["submitted"])
        self.controller = self.make_controller()
        self.assertTrue(self.controller.reconcile()["published"])
        self.assertTrue(self.controller.create()["published"])
        self.assertEqual(len(self.posts("/actions/newversion")), 1)
        self.assertEqual(len(self.posts("/actions/publish")), 1)

    def test_file_and_record_actions_use_supported_accept_media_types(self):
        self.server.owner_latest_only = True
        self.prepared()
        self.controller.publish()
        files = [(method, url, headers) for method, url, headers in self.server.calls
                 if method in ("POST", "PUT") and "/draft/files" in url]
        self.assertEqual(len(files), 6)
        for method, url, headers in files:
            self.assertEqual(headers["Accept"], "application/json")
            self.assertEqual(headers["Content-Type"], "application/octet-stream" if method == "PUT" else "application/json")
        records = [(method, url, headers) for method, url, headers in self.server.calls
                   if (method == "PUT" and url.endswith("/draft")) or url.endswith("/actions/publish")]
        self.assertEqual(len(records), 2)
        self.assertTrue(all(headers["Accept"] == M.VENDOR for method, url, headers in records))
        self.assertTrue(all(headers["If-Match"] == "7" for method, url, headers in records))

    def test_raw_verified_etag_preserved_and_publish_header_is_decimal(self):
        result=self.prepared()
        self.assertEqual(result["etag"],'"7"')
        self.controller.publish()
        self.assertEqual(self.posts("/actions/publish")[0][2]["If-Match"],"7")

    def test_invalid_or_mismatched_metadata_etag_stops_before_put(self):
        for index,tag in enumerate(('W/"7"','"-7"','"seven"','7','"07"','"8"')):
            with self.subTest(tag=tag):
                self.controller.create()
                self.server.draft_etag=tag
                reason="ETAG_REVISION_ID_DIFFERS" if tag=='"8"' else "ETAG_NOT_CANONICAL_QUOTED_DECIMAL"
                self.stop(reason,self.controller.prepare)
                self.assertFalse(any(method=="PUT" and url.endswith("/draft") for method,url,headers in self.server.calls))
                self.assertFalse(self.posts("/actions/publish"))
                self.reset_synthetic_case("invalid-etag-"+str(index))

    def test_invalid_publish_etag_stops_before_attempt_latch_and_post(self):
        self.prepared()
        self.server.draft_etag='W/"7"'
        self.controller.verify()
        self.stop("ETAG_NOT_CANONICAL_QUOTED_DECIMAL",self.controller.publish)
        self.assertFalse(self.controller.state.get("publish_post_attempted",False))
        self.assertFalse(self.posts("/actions/publish"))

    def test_metadata_400_json_saved_redacted_before_stop_and_no_publish(self):
        self.controller.create()
        self.server.metadata_error_response = response({"message":"Provider validation error " + FAKE_TOKEN,
             "errors":[{"field":"metadata.example", "message":"Unknown field"}],
             "access_token":FAKE_TOKEN,
             "detail":"See https://zenodo.org/api/records/1?credential=OTHER_PRIVATE_VALUE now"},status=400)
        self.stop("WRITE_REJECTED_400",self.controller.prepare)
        p = next((self.root / "evidence").glob("error_save_metadata_*.json"))
        text = p.read_text();error = json.loads(text)
        self.assertEqual(error["format"],"JSON")
        self.assertEqual(error["server_error"]["errors"][0]["message"],"Unknown field")
        self.assertNotIn(FAKE_TOKEN,text)
        self.assertNotIn("OTHER_PRIVATE_VALUE",text)
        self.assertIsNone(self.controller.state["pending"])
        self.assertEqual(len(self.server.files),29)
        self.assertFalse(self.posts("/actions/publish"))

    def test_metadata_only_retry_keeps_two_completed_uploads_and_exact_body(self):
        self.controller.create()
        self.server.metadata_error_response = response({"message":"Validation error"},status=400)
        self.stop("WRITE_REJECTED_400",self.controller.prepare)
        self.controller = self.make_controller()
        self.stop("WRITE_REJECTED_400",self.controller.prepare)
        self.assertEqual(sum(method == "PUT" and url.endswith("/content") for method,url,headers in self.server.calls),2)
        self.assertEqual(len(self.posts("/commit")),2)
        self.assertEqual(len([x for x in self.server.calls if x[0]=="POST" and x[1].endswith("/draft/files")]),2)
        journal=[json.loads(x) for x in (self.root / "evidence/JOURNAL.jsonl").read_text().splitlines()]
        intents=[x for x in journal if x["kind"]=="WRITE_INTENT" and x.get("operation")=="save_metadata"]
        self.assertEqual(len(intents),2)
        self.assertEqual(intents[0]["body_sha256"],intents[1]["body_sha256"])
        self.assertFalse(self.posts("/actions/publish"))

    def test_invalid_error_json_preserves_definite_400_without_raw_body(self):
        bodies=(b"<html>"+FAKE_TOKEN.encode()+b"</html>",b"["*70+b"0"+b"]"*70,
                b'{"value":1e309}',b"x"*(M.MAX_JSON+1))
        for index,body in enumerate(bodies):
            with self.subTest(index=index):
                self.controller.create()
                self.server.metadata_error_response = response(raw=body,status=400)
                self.stop("WRITE_REJECTED_400",self.controller.prepare)
                p=next((self.root / "evidence").glob("error_save_metadata_*.json"))
                error=json.loads(p.read_text())
                self.assertNotIn("server_error",error)
                self.assertNotIn(FAKE_TOKEN,p.read_text())
                self.assertIsNone(self.controller.state["pending"])
                self.assertFalse(self.posts("/actions/publish"))
                self.reset_synthetic_case("invalid-error-"+str(index))

    def test_error_500_diagnostic_keeps_pending_and_blocks_metadata_retry(self):
        self.controller.create()
        self.server.metadata_error_response = response({"message":"Uncertain server outcome"},status=500)
        self.stop("WRITE_OUTCOME_UNKNOWN_RECONCILE_ONLY",self.controller.prepare)
        self.assertEqual(self.controller.state["pending"]["kind"],"save_metadata")
        self.controller = self.make_controller()
        self.stop("COMPLETE_METADATA_READBACK_DIFFERS",self.controller.prepare)
        self.assertEqual(sum(method=="PUT" and url.endswith("/draft") for method,url,headers in self.server.calls),1)
        self.assertIsNotNone(self.controller.state["pending"])
        self.assertFalse(self.posts("/actions/publish"))

    def test_latest_only_unjournaled_draft_cannot_authorize_creation_or_adoption(self):
        self.controller.create()
        self.server.owner_latest_only = True
        self.controller.state["draft_id"] = None
        self.controller.save()
        self.controller = self.make_controller()
        self.stop("PUBLISHED_PRIOR_NOT_IN_OWNED_LIST", self.controller.create)
        self.assertEqual(len(self.posts("/actions/newversion")), 1)

    def test_latest_only_wrong_owner_family_state_or_id_are_rejected(self):
        cases = (("owner", M.EXPECTED_OWNER + 1, "OWNED_FAMILY_OWNER_DIFFERS"),
                 ("conceptrecid", "999", "PUBLISHED_PRIOR_NOT_IN_OWNED_LIST"),
                 ("submitted", True, "PUBLISHED_PRIOR_NOT_IN_OWNED_LIST"),
                 ("submitted", 0, "OWNED_FAMILY_STATE_UNKNOWN"),
                 ("id", int(DRAFT) + 1, "PUBLISHED_PRIOR_NOT_IN_OWNED_LIST"),
                 ("is_draft", False, "PUBLISHED_PRIOR_NOT_IN_OWNED_LIST"))
        for index, (field, value, reason) in enumerate(cases):
            with self.subTest(field=field, value=value):
                self.controller.create()
                owned = self.server.owned()[-1]
                owned[field] = value
                self.server.owner_pages = {1: [owned]}
                self.stop(reason, self.controller.reconcile)
                self.assertEqual(len(self.posts("/actions/newversion")), 1)
                self.assertFalse(self.posts("/actions/publish"))
                self.reset_synthetic_case("latest-invalid-" + str(index))

    def test_latest_only_published_record_requires_matching_submitted_state(self):
        self.server.owner_latest_only = True
        self.prepared()
        self.controller.publish()
        owned = self.server.owned()[0]
        owned["submitted"] = False
        self.server.owner_pages = {1: [owned]}
        self.stop("PUBLISHED_PRIOR_NOT_IN_OWNED_LIST", self.controller.reconcile)
        self.assertEqual(len(self.posts("/actions/newversion")), 1)
        self.assertEqual(len(self.posts("/actions/publish")), 1)

    def test_latest_only_unknown_publish_recovers_only_from_public_gets(self):
        self.server.owner_latest_only = True
        self.prepared()
        self.server.publish_fault = "unknown_after"
        self.stop("WRITE_OUTCOME_UNKNOWN_RECONCILE_ONLY", self.controller.publish)
        self.assertEqual(self.controller.state["pending"]["kind"], "publish")
        self.controller = self.make_controller()
        self.assertTrue(self.controller.reconcile()["published"])
        self.assertIsNone(self.controller.state["pending"])
        self.assertTrue(self.controller.verify()["published"])
        self.assertEqual(len(self.posts("/actions/publish")), 1)
        self.assertEqual(len(self.posts("/actions/newversion")), 1)

    def test_latest_only_unknown_create_without_journaled_id_remains_blocked(self):
        self.server.owner_latest_only = True
        self.server.create_fault = "unknown_after"
        self.stop("WRITE_OUTCOME_UNKNOWN_RECONCILE_ONLY", self.controller.create)
        self.controller = self.make_controller()
        self.stop("PUBLISHED_PRIOR_NOT_IN_OWNED_LIST", self.controller.create)
        self.assertIsNotNone(self.controller.state["pending"])
        self.assertIsNone(self.controller.state["draft_id"])
        self.assertEqual(len(self.posts("/actions/newversion")), 1)

    def test_missing_inherited_pin_blocks_adoption(self):
        self.server.import_drop = True
        self.stop("FILE_MEMBERSHIP_DIFFERS", self.controller.create)
        self.assertIsNone(self.controller.state["draft_id"])

    def test_complete_prepare_publish_once_and_persisted_redaction(self):
        result = self.prepared()
        self.assertEqual(sum(x["basis"] == "FRESH_COMPLETE_CONTENT_STREAM" for x in result["content"]), 29)
        self.assertEqual(sum(c[0] == "PUT" and c[1].endswith("/content") for c in self.server.calls), 2)
        final = self.controller.publish()
        self.assertTrue(final["published"])
        self.assertEqual(len(self.posts("/actions/publish")), 1)
        self.stop("EDITION_ALREADY_PUBLISHED_NO_REPEAT_POST", self.controller.publish)
        self.controller = self.make_controller()
        self.stop("EDITION_ALREADY_PUBLISHED_NO_REPEAT_POST", self.controller.publish)
        self.assertEqual(len(self.posts("/actions/publish")), 1)
        for path in (self.root / "evidence").glob("*"):
            self.assertNotIn(FAKE_TOKEN, path.read_text())

    def test_metadata_mismatch_and_extra_file_block_verify(self):
        self.prepared()
        self.server.draft_metadata["metadata"]["description"] = "<p>Changed content.</p>"
        self.stop("COMPLETE_METADATA_READBACK_DIFFERS", self.controller.verify)
        self.server.draft_metadata = clone(self.metadata)
        self.server.files["unexpected.txt"] = self.server.entry(pin("unexpected.txt", b"extra"))
        self.stop("FILE_MEMBERSHIP_DIFFERS", self.controller.verify)
        self.assertFalse(self.posts("/actions/publish"))

    def test_changed_inherited_id_requires_complete_fresh_stream(self):
        self.controller.create()
        name = self.baseline["files"][0]["filename"]
        self.server.files[name]["file_id"] = "changed-immutable-id"
        result = self.controller.prepare()
        row = next(x for x in result["content"] if x["filename"] == name)
        self.assertEqual(row["basis"], "FRESH_COMPLETE_CONTENT_STREAM")
        self.assertTrue(any(name + "/content" in url for url, headers in self.server.streams))
        self.server.payloads[name] += b"corrupted"
        self.stop("CONTENT_EXCEEDS_FROZEN_SIZE", self.controller.verify)

    def test_original_public_immutable_identity_change_blocks_all_writes(self):
        name = self.baseline["files"][0]["filename"]
        for field in ("file_id", "version_id"):
            with self.subTest(field=field):
                original = self.server.old_files[name][field]
                self.server.old_files[name][field] = "unexpected-original-public-identity"
                self.stop("ORIGINAL_PUBLISHED_IMMUTABLE_FILE_IDENTITY_CHANGED", self.controller.create)
                self.assertFalse(any(method != "GET" for method, url, headers in self.server.calls))
                self.server.old_files[name][field] = original

    def test_ambiguous_create_recovers_with_gets_and_no_second_post(self):
        self.server.create_fault = "unknown_after"
        self.stop("WRITE_OUTCOME_UNKNOWN_RECONCILE_ONLY", self.controller.create)
        self.assertEqual(self.controller.state["pending"]["kind"], "create_version")
        self.controller = self.make_controller()
        result = self.controller.create()
        self.assertEqual(result["draft_id"], DRAFT)
        self.assertIsNone(self.controller.state["pending"])
        self.assertEqual(len(self.posts("/actions/newversion")), 1)

    def test_malformed_success_create_recovers_without_repeat(self):
        self.server.create_fault = "malformed_after"
        with self.assertRaises(json.JSONDecodeError):
            self.controller.create()
        self.controller = self.make_controller()
        self.assertEqual(self.controller.create()["draft_id"], DRAFT)
        self.assertEqual(len(self.posts("/actions/newversion")), 1)

    def test_rejected_create_latch_survives_restart(self):
        self.server.create_fault = "reject_before"
        self.stop("WRITE_REJECTED_409", self.controller.create)
        self.controller = self.make_controller()
        self.stop("NEWVERSION_POST_ALREADY_ATTEMPTED_RECONCILE_ONLY", self.controller.create)
        self.assertEqual(len(self.posts("/actions/newversion")), 1)

    def initial_legacy_accept_406(self):
        """Reproduce the prior vendor-Accept error against the header-aware server."""
        original = self.controller.write
        def prior_header(kind, url, data=None, headers=None, **context):
            h = dict(headers or {})
            if kind == "create_version":
                h["Accept"] = M.VENDOR
            return original(kind, url, data=data, headers=h, **context)
        self.controller.write = prior_header
        self.stop("WRITE_REJECTED_406", self.controller.create)
        self.assertFalse(self.server.exists)
        self.assertTrue(self.controller.state["create_post_attempted"])
        self.assertIsNone(self.controller.state["pending"])
        self.controller = self.make_controller()

    def reset_synthetic_case(self, label):
        # Retain each synthetic journal; never touch the live execution directory.
        (self.root / "evidence").rename(self.root / ("evidence-preserved-" + label))
        self.server = FakeServer(self.baseline["files"], self.inventory["new_files"], self.metadata, self.payloads)
        self.controller = self.make_controller()

    def test_legacy_newversion_uses_json_accept(self):
        self.controller.create()
        self.assertEqual(self.posts("/actions/newversion")[0][2]["Accept"], "application/json")

    def test_definitive_406_one_explicit_corrected_retry_after_reconciliation(self):
        self.initial_legacy_accept_406()
        result = self.controller.reconcile()
        self.assertTrue(result["owner_listing_exhausted"])
        self.assertIsNone(result["draft_id"])
        self.assertEqual(len(self.posts("/actions/newversion")), 1)
        self.stop("NEWVERSION_POST_ALREADY_ATTEMPTED_RECONCILE_ONLY", self.controller.create)
        result = self.controller.create(retry_create_406=True)
        self.assertEqual(result["draft_id"], DRAFT)
        self.assertTrue(self.controller.state["create_post_attempted"])
        self.assertTrue(self.controller.state["create_406_retry_attempted"])
        self.assertEqual([x[2]["Accept"] for x in self.posts("/actions/newversion")], [M.VENDOR, "application/json"])
        self.controller = self.make_controller()
        self.assertEqual(self.controller.create(retry_create_406=True)["draft_id"], DRAFT)
        self.assertEqual(len(self.posts("/actions/newversion")), 2)

    def test_rejected_corrected_406_retry_latch_survives_restart(self):
        self.initial_legacy_accept_406()
        self.server.create_reject_status = 409
        self.stop("WRITE_REJECTED_409", lambda: self.controller.create(retry_create_406=True))
        self.controller = self.make_controller()
        self.stop("CREATE_406_RETRY_ALREADY_ATTEMPTED", lambda: self.controller.create(retry_create_406=True))
        self.assertEqual(len(self.posts("/actions/newversion")), 2)
        self.assertTrue(self.controller.state["create_post_attempted"])
        self.assertTrue(self.controller.state["create_406_retry_attempted"])

    def test_non406_rejections_never_qualify_for_explicit_retry(self):
        for status in (400, 409, 429):
            with self.subTest(status=status):
                self.server.create_reject_status = status
                self.stop("WRITE_REJECTED_" + str(status), self.controller.create)
                self.controller = self.make_controller()
                self.stop("CREATE_406_RETRY_REQUIRES_DEFINITIVE_LAST_406", lambda: self.controller.create(retry_create_406=True))
                self.assertEqual(len(self.posts("/actions/newversion")), 1)
                self.assertFalse(self.controller.state.get("create_406_retry_attempted", False))
                self.reset_synthetic_case(str(status))

    def unknown_creation_before_remote_effect(self, transport_error=False):
        original = self.server.request
        def unknown(method, url, data=None, headers=None):
            if method == "POST" and url.endswith("/actions/newversion"):
                self.server.calls.append((method, url, clone(headers or {})))
                if transport_error:
                    raise OSError("offline ambiguous transport")
                return response({}, status=500)
            return original(method, url, data=data, headers=headers)
        self.server.request = unknown

    def test_unknown_initial_500_and_transport_remain_blocked(self):
        for transport_error in (False, True):
            with self.subTest(transport_error=transport_error):
                self.unknown_creation_before_remote_effect(transport_error)
                reason = "WRITE_OUTCOME_UNKNOWN_RECONCILE_ONLY"
                self.stop(reason, self.controller.create)
                self.controller = self.make_controller()
                self.stop("AMBIGUOUS_WRITE_NOT_RECONCILED", lambda: self.controller.create(retry_create_406=True))
                self.assertIsNotNone(self.controller.state["pending"])
                self.assertEqual(len(self.posts("/actions/newversion")), 1)
                self.assertFalse(self.controller.state.get("create_406_retry_attempted", False))
                self.reset_synthetic_case("unknown-" + str(transport_error))

    def test_unknown_corrected_retry_preserves_both_latches_and_never_third_post(self):
        for transport_error in (False, True):
            with self.subTest(transport_error=transport_error):
                self.initial_legacy_accept_406()
                self.unknown_creation_before_remote_effect(transport_error)
                self.stop("WRITE_OUTCOME_UNKNOWN_RECONCILE_ONLY", lambda: self.controller.create(retry_create_406=True))
                self.controller = self.make_controller()
                self.stop("AMBIGUOUS_WRITE_NOT_RECONCILED", lambda: self.controller.create(retry_create_406=True))
                self.assertEqual(len(self.posts("/actions/newversion")), 2)
                self.assertTrue(self.controller.state["create_post_attempted"])
                self.assertTrue(self.controller.state["create_406_retry_attempted"])
                self.assertIsNotNone(self.controller.state["pending"])
                self.reset_synthetic_case("retry-unknown-" + str(transport_error))

    def test_retry_rejects_missing_duplicate_or_intervening_create_evidence(self):
        for case in ("missing", "duplicate", "unknown"):
            with self.subTest(case=case):
                self.initial_legacy_accept_406()
                p = self.root / "evidence/JOURNAL.jsonl"
                rows = [json.loads(line) for line in p.read_text().splitlines()]
                intent = next(x for x in rows if x.get("kind") == "WRITE_INTENT")
                if case == "missing":
                    rows.remove(intent)
                elif case == "duplicate":
                    rows.append(clone(intent))
                else:
                    rows.append({"kind":"WRITE_OUTCOME_UNKNOWN", "operation":"create_version"})
                p.write_text("".join(json.dumps(x) + "\n" for x in rows))
                self.stop("CREATE_406_RETRY_REQUIRES_DEFINITIVE_LAST_406", lambda: self.controller.create(retry_create_406=True))
                self.assertEqual(len(self.posts("/actions/newversion")), 1)
                self.reset_synthetic_case("journal-" + case)

    def test_retry_rejects_wrong_original_url_or_payload_pin(self):
        for field, value in (("url", M.ORIGIN + "/api/deposit/depositions/999/actions/newversion"),
                             ("body_sha256", "0" * 64), ("body_bytes", 1), ("body_bytes", 0.0)):
            with self.subTest(field=field, value=value):
                self.initial_legacy_accept_406()
                p = self.root / "evidence/JOURNAL.jsonl"
                rows = [json.loads(line) for line in p.read_text().splitlines()]
                next(x for x in rows if x.get("kind") == "WRITE_INTENT")[field] = value
                p.write_text("".join(json.dumps(x) + "\n" for x in rows))
                self.stop("CREATE_406_RETRY_ORIGINAL_INTENT_DIFFERS", lambda: self.controller.create(retry_create_406=True))
                self.assertEqual(len(self.posts("/actions/newversion")), 1)
                self.reset_synthetic_case("pin-" + field + "-" + str(value).replace("/", "_")[-10:])

    def test_ambiguous_publication_recovers_from_public_gets(self):
        self.prepared()
        self.server.publish_fault = "unknown_after"
        self.stop("WRITE_OUTCOME_UNKNOWN_RECONCILE_ONLY", self.controller.publish)
        self.controller = self.make_controller()
        result = self.controller.reconcile()
        self.assertTrue(result["published"])
        self.assertIsNone(self.controller.state["pending"])
        self.stop("EDITION_ALREADY_PUBLISHED_NO_REPEAT_POST", self.controller.publish)
        self.assertTrue(self.controller.verify()["published"])
        self.assertEqual(len(self.posts("/actions/publish")), 1)

    def test_rejected_publication_latch_survives_restart(self):
        self.prepared()
        self.server.publish_fault = "reject_before"
        self.stop("WRITE_REJECTED_409", self.controller.publish)
        self.controller = self.make_controller()
        self.stop("PUBLISH_POST_ALREADY_ATTEMPTED_RECONCILE_ONLY", self.controller.publish)
        self.assertEqual(len(self.posts("/actions/publish")), 1)

    def edit_after_publish_verification(self, change):
        """Make the final GET observe an edit after the full verify returns."""
        original = self.controller.verify
        def interleaved(published=None):
            result = original(published)
            change()
            return result
        self.controller.verify = interleaved

    def test_metadata_edit_between_verify_and_final_get_blocks_publish(self):
        self.prepared()
        def change():
            self.server.draft_metadata["metadata"]["description"] = "<p>Concurrent metadata edit.</p>"
            self.server.draft_etag = '"draft-edited"'
        self.edit_after_publish_verification(change)
        self.stop("COMPLETE_METADATA_READBACK_DIFFERS", self.controller.publish)
        self.assertFalse(self.posts("/actions/publish"))
        self.assertFalse(self.controller.state.get("publish_post_attempted", False))

    def test_completed_file_edit_between_verify_and_final_get_blocks_publish(self):
        self.prepared()
        name = self.inventory["new_files"][0]["filename"]
        def change():
            replacement = b"X" * len(self.server.payloads[name])
            self.server.payloads[name] = replacement
            self.server.files[name]["checksum"] = "md5:" + hashlib.md5(replacement).hexdigest()
            self.server.draft_etag = '"draft-edited"'
        self.edit_after_publish_verification(change)
        self.stop("FILE_MD5_DIFFERS", self.controller.publish)
        self.assertFalse(self.posts("/actions/publish"))
        self.assertFalse(self.controller.state.get("publish_post_attempted", False))

    def test_etag_only_edit_between_verify_and_final_get_blocks_publish(self):
        self.prepared()
        self.edit_after_publish_verification(lambda: setattr(self.server, "draft_etag", '"draft-edited"'))
        self.stop("DRAFT_CHANGED_SINCE_COMPLETE_VERIFICATION", self.controller.publish)
        self.assertFalse(self.posts("/actions/publish"))
        self.assertFalse(self.controller.state.get("publish_post_attempted", False))

    def test_publish_requires_saved_verified_receipt_and_untampered_file(self):
        self.controller.create()
        self.stop("VERIFY_STAGE_REQUIRED_BEFORE_PUBLISH", self.controller.publish)
        self.controller.prepare()
        receipt = self.root / "evidence/VERIFIED_DRAFT.json"
        receipt.write_text(receipt.read_text() + " ")
        self.stop("VERIFICATION_RECEIPT_CHANGED", self.controller.publish)
        self.assertFalse(self.posts("/actions/publish"))

    def test_numeric_bool_and_unrecognized_decoration_are_rejected(self):
        actual = clone(self.metadata)
        actual["access"]["embargo"]["active"] = 0
        self.stop("ACCESS_READBACK_DIFFERS", lambda: M.metadata_equal(actual, self.metadata))
        actual = clone(self.metadata)
        actual["metadata"]["creators"][0]["person_or_org"].update({"id": "arbitrary", "title": "unrecognized decoration"})
        self.stop("COMPLETE_METADATA_READBACK_DIFFERS", lambda: M.metadata_equal(actual, self.metadata))
        files = clone(self.server.old_files)
        files[self.baseline["files"][0]["filename"]]["size"] = float(self.baseline["files"][0]["bytes"])
        self.stop("FILE_SIZE_DIFFERS", lambda: M.assert_pins(files, self.baseline["files"]))

    def test_owner_reconciliation_reads_later_page_before_adoption(self):
        self.server.exists = True
        self.server.files = clone(self.server.old_files)
        for name,f in self.server.files.items():f["links"]["content"]=f"{M.ORIGIN}/api/records/{DRAFT}/draft/files/{name}/content"
        self.server.owner_pages = {1: [self.server.owned()[0]] + [
            {"id": 25000000 + i, "conceptrecid": "unrelated", "submitted": True}
            for i in range(99)], 2: [self.server.owned()[1]]}
        result = self.controller.reconcile()
        self.assertEqual(result["draft_id"], DRAFT)
        pages = [url for method, url, headers in self.server.calls if method == "GET" and "/api/deposit/depositions?" in url]
        self.assertEqual(len(pages), 2)
        self.assertFalse(any(method != "GET" for method, url, headers in self.server.calls))

    def test_owner_pagination_repeat_blocks_creation(self):
        first = [self.server.owned()[0]] + [
            {"id": 25000000 + i, "conceptrecid": "unrelated", "submitted": True}
            for i in range(99)]
        self.server.owner_pages = {1: first, 2: first}
        self.stop("OWNER_PAGINATION_REPEATED_IDS", self.controller.create)
        self.assertFalse(self.posts("/actions/newversion"))

    def test_owner_pagination_bound_never_implies_absence(self):
        self.server.owner_pages = {page: [
            {"id": 25000000 + page * 100 + i, "conceptrecid": "unrelated", "submitted": True}
            for i in range(100)] for page in range(1, 21)}
        self.stop("OWNER_PAGINATION_NOT_EXHAUSTED", self.controller.create)
        pages = [url for method, url, headers in self.server.calls if method == "GET" and "/api/deposit/depositions?" in url]
        self.assertEqual(len(pages), 20)
        self.assertFalse(self.posts("/actions/newversion"))

if __name__ == "__main__":
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ControllerFixtures)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    report = {"status": "PASS" if result.wasSuccessful() else "FAIL", "tests_run": result.testsRun,
              "failures": [{"test": str(t), "traceback": trace} for t, trace in result.failures],
              "errors": [{"test": str(t), "traceback": trace} for t, trace in result.errors],
              "skipped": list(result.skipped), "controller_sha256": M.sha(SOURCE.read_bytes()),
              "suite_sha256": M.sha(Path(__file__).read_bytes()), "real_network": "FORBIDDEN_BY_SOCKET_SENTINEL",
              "remote_mutations": 0, "credentials": "FAKE_TOKEN_ONLY", "baseline": "SYNTHETIC_INJECTED_27_SMALL_FILES_ALL29_FAKE_STREAMED"}
    (HERE / "OFFLINE_TEST_RESULTS.json").write_text(json.dumps(report, indent=2) + "\n")
    raise SystemExit(0 if result.wasSuccessful() else 1)
