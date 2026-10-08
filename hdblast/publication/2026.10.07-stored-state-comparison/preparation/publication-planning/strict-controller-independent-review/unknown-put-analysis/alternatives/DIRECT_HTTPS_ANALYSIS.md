# Direct HTTPS early-response recovery: local-source analysis

Scope: installed CPython 3.12.14 standard-library source only. No network, controller import or execution, repository/controller/state edits, or credential access occurred. This is a source-supported design option, not an empirical claim that a particular server response remains readable after a TLS write failure.

## Finding

A fresh `http.client.HTTPSConnection` can preserve an early final response that `urllib.request` would discard when the request-body send raises `OSError`. After stopping the body immediately, one `getresponse()` call on the same existing connection can attempt to read that response. The call does not connect, send, authenticate again, redirect, or repeat the PUT. An unsuccessful recovery remains an unknown operation outcome.

The installed `http/client.py` establishes the relevant sequence:

- `endheaders()` sets the connection state to `Request-sent` before `_send_output()` (lines 1370–1381).
- `_send_output()` sends headers and then the body chunks (lines 1132–1184).
- `send()` calls the retained socket's `sendall()` and catches only `TypeError`; it does not close/reset the socket on `OSError` (lines 1077–1112).
- `getresponse()` requires `Request-sent`, builds the response from `self.sock`, and calls `response.begin()`; there is no call to `connect()`, `send()`, or `request()` in it (lines 1434–1495).
- `HTTPResponse` reads from `sock.makefile("rb")` (line 281). Its parser treats EOF without a status as `RemoteDisconnected`, validates the status line, then reads headers (lines 302–365).

By contrast, installed `urllib/request.py` lines 1342–1351 convert an `OSError` from `h.request()` into `URLError`, close the connection, and raise before reaching `h.getresponse()`.

## Required boundaries

Use a new connection for exactly one request, complete the verified TLS connection before headers, and distinguish the body-send phase explicitly. Separate header sending (`endheaders()` without a body) and local file reading from the catch around each `conn.send(chunk)`. An `OSError` raised by opening/reading the local file is not evidence of a server rejection. Do not resume a failed chunk or continue remaining body bytes.

Retain the socket identity after the TLS connection. Before recovery, require that the same socket is still present; never call `connect()`, `request()`, or `send()` from a recovery branch. Setting `conn.auto_open = False` after explicit connection setup provides an additional guard against a later accidental send reconnect. Attempt `getresponse()` once, then close the connection in all cases. A close-delimited response remains readable through `HTTPResponse`'s separate file reference when `getresponse()` closes its connection handle; this is normal library behavior, not a new connection.

Use the default verified HTTPS context, or an explicit `ssl.create_default_context()` with `CERT_REQUIRED`, hostname checking, and ALPN restricted to `http/1.1`. Installed `http/client.py` lines 845–855 and 1507–1528 select the verified default context and pass the origin hostname to `wrap_socket`; installed `ssl.py` lines 694–724 create a client context, require certificates, check the hostname, and load trusted roots. Do not weaken verification after an error. Direct `HTTPSConnection` has no redirect or HTTP authentication-challenge handler; preserve that property. Reject redirection responses instead of following them.

Keep the bearer token solely in the in-memory Authorization header. Do not put it in the URL, command line, diagnostics, or files. Keep `debuglevel` zero: installed `send()` prints the raw data at a positive debug level. Existing audit hooks can receive raw send buffers, so a transport cannot claim to conceal the header from arbitrary in-process auditing. Avoid TLS key logging: `ssl.create_default_context()` can honor `SSLKEYLOGFILE` unless Python ignores the environment (`ssl.py` lines 725–729); isolated Python or explicitly clearing the context's keylog filename prevents that particular output path.

Keep finite socket timeouts and an independent whole-operation deadline. A socket timeout is not a whole-response deadline; a slowly progressing peer can extend total time. Bound captured headers/body and parse the API response under the existing strict receipt rules. Library header limits are 65,536 bytes per line, 100 header lines, and 100 interim `100` responses in this installed build (`http/client.py` lines 111–119). These are not a substitute for a tighter application receipt bound. Accept only a valid final status (at least 200), with all required receipt bytes and headers captured; a parsed status with an incomplete required body is not a complete definitive receipt.

## Limits

The socket may have been reset, closed, timed out, or left in an unusable TLS state; the peer's response may never have arrived or may no longer be recoverable. A second read can fail or time out. `SSLSocket.sendall()` itself loops over partial TLS writes (`ssl.py` lines 1211–1225), but does not expose how many body bytes reached the server when it fails. Neither a byte counter nor the exception proves the application accepted or rejected the operation.

The recovery does not make an unknown PUT safe to repeat. Even a valid HTTP rejection needs the application's specific definitive-rejection policy before it could authorize any separately journaled attempt. A complete successful response recovered after a send error must likewise be assessed using the normal success receipt and readback rules. No change to a publisher's unknown-operation latch follows from this source analysis alone.

## Source pins

- `/opt/codex/runtimes/codex-primary-runtime/dependencies/python/lib/python3.12/http/client.py`: SHA-256 `c2c0694155837d72fc710a4138e5d1e53ae1a652b117a2e23a5b1b0db527a0f9`.
- `/opt/codex/runtimes/codex-primary-runtime/dependencies/python/lib/python3.12/ssl.py`: SHA-256 `db62eeaae59c6e7087dcafb28087cd0bb092087381edc556f97c8cb614484380`.
- `/opt/codex/runtimes/codex-primary-runtime/dependencies/python/lib/python3.12/urllib/request.py`: SHA-256 `6ffc7d993de64319a40fad631966bd1fd643d355a5c154c1a8b1821f49e5c470`.
