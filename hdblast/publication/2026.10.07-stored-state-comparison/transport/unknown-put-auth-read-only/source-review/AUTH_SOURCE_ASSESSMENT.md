# Saved PUT 401: source and read-only evidence assessment

The accepted helper's curl configuration has no identified quoting/parser fault. The read-only GET comparison does not reproduce an authentication effect from removing `ZENODO_ACCESS_TOKEN` from curl's child environment. The exact cause of the saved PUT's managed-integration 401 remains unresolved.

Accepted source SHA256: `39ec701ae3133b699bbad1f74031c549a0aedbc2116109cf5be199d3c7b0916d` (19,883 bytes). `token_config` validates a nonempty printable ASCII value, escapes backslashes and then double quotes, places the bearer in one double-quoted curl config header, and ends that option with one newline. Installed curl documentation explicitly supports this quoting and both escape sequences, and `--config -` reads stdin. The command puts `-q` first, supplies one fixed URL, and contains no other Authorization option. Curl's child environment removes `ZENODO_ACCESS_TOKEN` and `SSLKEYLOGFILE`; the parent token variable remains unchanged for generating that stdin configuration. Credential values, proxy URLs and Authorization headers were not inspected or printed.

The parent-run comparison uses the same backslash-then-quote config expression. Its three saved receipts independently confirm:

| Read-only variant | HTTP | Body bytes | Body SHA256 |
|---|---:|---:|---|
| urllib vendor GET | 200 | 42,180 | `fc5fe357c2171567de78d5baa874a420e3c1a5782508097aebafe1997cb86035` |
| curl vendor GET, credential variable retained | 200 | 42,180 | same |
| curl vendor GET, credential variable removed | 200 | 42,180 | same |

Both curl GETs retain proxy/trust settings, remove TLS key logging, deliver the header via stdin, report no redirects, TLS verification result zero and curl exit zero. Thus neither config quoting nor child credential-variable removal produced a differential failure in the demonstrated GET path. This observation does not test the large streaming PUT path or prove its integration behaves identically.

The saved one-attempt PUT receipt and raw failure evidence independently confirm final HTTP 401, curl exit zero, a 76-byte body matching `Unauthorized: authentication failed for integration codex-secret-zenodo.org`, `x-at-upstream-error: false`, 24,051,375 bytes in curl's uploaded counter, and 4.407293 seconds. Exit zero is consistent with curl's default behavior for a received HTTP error; it is not upload success. The counter describes curl's transfer observation, not proof of complete upstream receipt. The header's integration-specific meaning and exact origin are not independently documented here, so the flag does not establish what Zenodo received or processed.

The integration-named rejection is evidence of an authentication failure along the managed upload route. It does not prove the supplied bearer was malformed, the integration credential was invalid, removal of an environment variable caused it, or a particular upstream component rejected it. The matching successful GETs support leaving the accepted config unchanged unless additional evidence isolates a PUT-specific fault.

The accepted helper, original state, original pending intent and recovery latch remain unchanged by this review. Its saved result retains `pending_cleared: false`, `commits: 0`, `publishes: 0`, and reconciliation-required status. This review performed no network request, helper/controller main or credential access. It read only pinned source and existing selected diagnostic evidence. Full source/evidence pins are in `AUTH_SOURCE_ASSESSMENT.json`.
