# Original passing audit with doubled canonical label

The script and normal/optimized JSON outputs here are the actual original bytes,
preserved before the factor correction. Original source SHA-256 is
`833188aed2d0f6e41c21f7d6f5ca26f6dbc0651cfdd6746da04179db82eb9e9d`. Both original audits passed.

`ORIGINAL_NORMAL_TOOL_STDOUT.log` and `ORIGINAL_OPTIMIZED_TOOL_STDOUT.log`
are verbatim transcriptions of stdout returned by the actual exec tool calls;
they were saved afterward, not shell-redirected during those executions. The
exit files record the actual tool-return exit code 0. They are not new reruns.

The scalar called `canonical_energy_integral` in these original records was
the variation of |v'|^2+k^2|v|^2, twice the canonical Hamiltonian. V2 adds the
physical factor 1/2 and reports 3.153086115117045e-17 instead of
6.30617223023409e-17. The frozen producer, stress results, data and all
acceptance outcomes are unchanged. This was a post-run label/normalization
correction to a passing audit, not a failed physical experiment.
