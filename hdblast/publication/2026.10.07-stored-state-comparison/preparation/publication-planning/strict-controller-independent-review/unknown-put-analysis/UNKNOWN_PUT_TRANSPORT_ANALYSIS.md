# Unknown PUT transport analysis

The installed CPython 3.12.14 request path can lose access to an early definitive HTTP response when a request-body send fails. The socket-free reproduction passes 15 cases in normal Python and the same 15 under `-O`. This demonstrates a runtime mechanism; it does not prove that the real PUT was rejected, identify a status, or establish the real server/proxy cause.

## Source-supported mechanism

`urllib.request.AbstractHTTPHandler.do_open()` calls `h.request(...)` before `h.getresponse()` (`urllib/request.py`, lines 1342–1351). `HTTPConnection.request()` sends headers and the entire body before it returns. `_send_output()` sends a bytes-like body as a single buffer to `send()`, and `send()` calls the retained socket's `sendall()` (`http/client.py`, lines 1077–1112, 1132–1184, 1370–1432). `SSLSocket.sendall()` itself loops until all TLS writes finish (`ssl.py`, lines 1211–1225). No concurrent HTTP response parser runs during that synchronous send.

If that body send raises `BrokenPipeError`, `ConnectionResetError`, `SSLEOFError`, a timeout, or another `OSError`, `do_open()` wraps it in `URLError`, closes the connection, and raises before `getresponse()` can parse a response. A response already received into the underlying connection can therefore go unread. The current controller's `Transport.request()` catches that exception in its general exception branch and raises `Stop('TRANSPORT_OUTCOME_UNKNOWN')`; only a parsed `HTTPError` reaches its branch that returns a definitive status and bounded error body. `Controller.write()` already saved the pending intent before transport, and preserves it when the transport raises. This matches an absent `WRITE_RESPONSE` but does not determine the remote effect.

A 240,135,519-byte body makes this mechanism relevant because the send runs over many lower-level writes. Size and a roughly two-second failure alone cannot distinguish early rejection, TLS/socket trouble, proxy trouble, a network reset, or other failures. The controller deliberately suppresses raw exception details; no retrospective HTTP status can be reconstructed from this error name.

## Manufactured reproduction

`offline_early_response_probe.py` authenticates the exact current controller source and compiles only `Transport.request` plus its required nonmutating helpers into an isolated harness; it never imports the controller module or runs its main, write, reconcile, or publication methods. A lazy anonymous mapping supplies an exact 240,135,519-byte memoryview. The fake socket measures its length without reading/uploading it. Socket construction and connection APIs are patched to reject all real network calls.

For each of four manufactured send exceptions, an already queued HTTP 413 is never read by urllib; the isolated source request method reports `TRANSPORT_OUTCOME_UNKNOWN`; and a direct same-connection response read recovers that manufactured 413 without another PUT. Additional controls expose 413 after a successful body send, reject redirect following, and keep the outcome unknown when no response exists. Both 15-case runs pass. No real TLS handshake, ZIP read, secret access, request, localhost connection, or remote operation occurred. The fake queued 413 is an example, not evidence of a real 413.

## Source-supported alternative without repeating PUT

A narrowly written direct `http.client.HTTPSConnection` transport can send headers once, stream body chunks, and stop immediately on a body-send error. If the original verified TLS socket is still retained, it can call `getresponse()` once on that same connection before closing it. Installed `getresponse()` reads from the existing socket and never calls `connect()`, `send()`, or `request()`. This has a source-supported no-repeat property and the manufactured recovery passes. A failed recovery remains unknown; a TCP reset, failed TLS state, missing response, or partial response may prevent recovery entirely. It cannot recover the already closed real urllib connection after the fact.

Use a default verified TLS context, hostname checking, the fixed HTTPS origin and exact URL, no redirects or authentication negotiation, no debug/raw credential logs, bounded response capture and operation deadlines, and only an in-memory bearer header. Explicitly distinguish body-send errors from local file reads or pre-request connection errors. After explicit TLS setup, disable `auto_open`; never reconnect or restart the request in a recovery branch. Preserve the existing persistent attempt latch and GET hash reconciliation. This is a prospective transport design, not authorization to repeat the pending request.

A transport that reads responses while streaming the body can expose an early rejection sooner. Curl/libcurl provides that kind of transfer loop, but the proposed strict single-PUT route needs additional scrutiny below.

## Installed curl: Expect route is uncertified

Installed `/usr/bin/curl` reports curl/libcurl 8.14.1, Debian security build 8.14.1-2+deb13u4, OpenSSL 3.5.7, and HTTP/2 support. It links the pinned `/usr/lib/x86_64-linux-gnu/libcurl.so.4.8.0`. The installed library contains these exact static diagnostic strings:

- `Got HTTP failure 417 while waiting for a 100`
- `Got HTTP failure 417 while sending data`
- `Need to rewind upload for next request`
- `REFUSED_STREAM, retrying a fresh connect`
- `Connection died, retrying a fresh connect (retry count: %d)`

`INSTALLED_CURL_EVIDENCE.json` retains complete binary hashes and byte offsets for these strings. They establish that this installed build includes relevant 417, rewind and implicit retry machinery. They do not prove the branch conditions or whether a particular invocation would repeat a PUT. The captured local manual describes `--retry 0` as disabling the tool's transient-transfer retry feature; it does not state that it disables every libcurl internal resend. Local implementation source was not available for proving the exact 417 branch. A persistent process-level attempt latch cannot prevent an internal repeat within that process. The proposed `Expect: 100-continue` route therefore remains uncertified for the strict one-wire-PUT requirement.

To narrow known resend routes, use a fresh process, `-q`/`--disable` as the first argument, exactly one fixed URL, explicit `--http1.1`, `--retry 0`, no redirects, no challenge-based origin/proxy authentication options, and suppress Expect with the empty `Expect:` header. A fresh one-URL process removes prior libcurl connection reuse; `--no-keepalive` only disables TCP probes. These constraints narrow the concern but do not prove all internal resend branches absent. If the environment requires its proxy, retaining proxy environment settings is an explicit scope limit: this review cannot establish the proxy's upstream retry behavior. Avoid `--location`, insecure TLS flags, tracing/verbose logs, and implicit curlrc options. Stream a pinned file through the upload interface rather than assuming `--data-binary @file` is a streaming interface. Keep the token in a properly escaped stdin configuration, never arguments, URLs or saved configuration, and redact it from retained response diagnostics. No sidecar was run or reviewed here.

## What remains unknown

The parent reports that the real draft is 23228395, the original PUT lacks `WRITE_RESPONSE`, and fresh GET diagnostics leave the pending intent in place with no completed size/checksum and unavailable content. Those are parent-provided observations, not independent network evidence from this review. None of the source or manufactured results establishes whether the remote service rejected or partially received the real upload. Root's GET reconciliation, request authorization, transport selection and any independently latched attempt remain outside this read-only analysis.

The controller source pin is `aa7e6cfe271559731ed3b45c5339f8eb612ebe6f017ddafa5fdadad9579e084b`. Complete artifact and local runtime source pins are recorded in `TRANSPORT_ANALYSIS_RECEIPT.json`.
