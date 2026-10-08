# Transport-options source audit

Candidate SHA-256 verified: `4db55b72c33287c04e686365af4cd9bbe33b7ad734b0bbc5d516c8985733db8c`.

Scope: source-only inspection of this snapshot and previously captured installed curl documentation. No helper/controller imports or execution, network calls, environment values, state, or assets were read.

## Correct controls

Lines 182–192 put `-q` first, use exactly one fixed HTTPS URL, force HTTP/1.1, suppress Expect with `Expect:`, set zero redirects/retries, and contain no `-L`, insecure TLS flag, HTTP authentication negotiation, token-bearing URL, or token-bearing argv. Explicit stdin config at lines 194–196 escapes backslashes and quotes and excludes CR/LF/control/non-ASCII values, so the config does not permit a newline-injected option. TLS verification remains enabled; inherited CA environment may select trust roots.

These controls avoid the documented redirect/negotiation and Expect-related concerns. Installed curl CLI documentation still does not prove a universal prohibition on every libcurl internal resend; this audit does not expand that claim.

## Credential-safety blockers

1. **Redaction after JSON encoding is unsound for token values that the helper permits.** Lines 63–66 serialize first and then replace literal token bytes. `token_config` permits quote/backslash characters. For example, a body string containing a token with a quote is JSON-escaped by `json.dumps`, so the serialized bytes no longer contain the literal token and the replacement misses it. `FRESH_MISSING_CONTENT.json` stores an observed response string through this function (line 279). Redact strings/keys recursively before serialization and retain the existing literal-byte defense, or constrain the accepted token grammar appropriately. A regression should cover quotes and backslashes.

2. **The curl child inherits the bearer environment variable and optional TLS key logging.** `Popen` at line 212 has no explicit `env`, so the token remains available as `ZENODO_ACCESS_TOKEN` in addition to the stdin header. The captured installed curl manual lines 6789–6796 says `SSLKEYLOGFILE` writes TLS secrets and explicitly includes this OpenSSL backend. Thus stdin is not the only credential path, and curl can write session secrets to an inherited file path. Construct a child environment that removes the bearer variable and `SSLKEYLOGFILE`, while preserving the specifically required proxy and trust configuration. No environment values were inspected here; this is a source-permitted path, not a claim that key logging is currently configured.

3. **Persisted raw stderr does not cover inherited proxy credentials.** Lines 241 and 293–298 replace only the literal Zenodo token before saving curl body/stderr/header bytes. Inherited proxy URLs may contain independent credentials; curl errors are not guaranteed to omit those URLs. JSON-escaped token echoes also evade the literal raw-byte replacement. Preserve only needed sanitized diagnostics/metrics or explicitly redact the known proxy credential representations in memory before writing. Do not expose proxy environment values during review. The installed manual lines 6676–6681 confirms proxy environment variables act like `--proxy`; lines 6740–6787 confirm inherited CA configuration also takes effect.

The first issue is directly demonstrable from the allowed token grammar and serialization order. The latter issues identify missing guarantees in the environment/output paths; their current environmental activation was intentionally not inspected.
