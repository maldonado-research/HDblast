"""Verify short quoted passages against retrieved text; no numerical sources."""
import hashlib
import json
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent
QUOTES = [
    ("raw/flint_v3_6_0_arb.rst", "opening enclosure contract", "The result of an (approximate) operation done on :type:`arb_t` variables is a ball which contains the result of the (mathematically exact) operation applied to any choice of points in the input balls.", "Enclosure applies to declared input balls, not an unrepresented ideal input."),
    ("raw/flint_v3_6_0_arb.rst", "arb_set_d note", "Be cautious when using :func:`arb_set_d` as it does not impose any error bounds and will only convert a ``double`` to an ``arb_t``.", "Native-double promotion does not add a decimal-real input error bound."),
    ("raw/flint_v3_6_0_acb_calc.rst", "acb_calc_integrate order=1 contract", "*func* must verify that *f* is holomorphic on this domain (and output a non-finite value if it is not).", "An integrand's domain validation is a caller obligation."),
    ("raw/flint_v3_6_0_acb_calc.rst", "integration accuracy goals", "These parameters are only guidelines; the cumulative error may be larger than both the prescribed absolute and relative error goals, depending on the number of subdivisions, cancellation between segments of the integral, and numerical errors in the evaluation of the integrand.", "Accept the actual returned enclosure radius, not the configured tolerance alone."),
    ("raw/johansson_integration_1802_07942v1.txt", "Sec.3.2, printed p.6", "Tail bounds must then be added based on symbolic knowledge about f .", "Manual truncation requires independently justified omitted-integral bounds."),
    ("raw/mottola_ward_2607_18180v1.txt", "Sec.II, printed p.7, after Eq.(2.21)", "which contains a second δ-function contact term, in addition to that in (2.13).", "The complete response has more than the nonlocal stress correlator."),
    ("raw/borinsky_collider_2609_10673v1.txt", "Sec.3.2, printed p.14, footnote12", "We use deterministic tensor-product tanh–sinh quadrature, monitor convergence under successive refinement, and compare with the independent tropical Monte Carlo calculation of Appendix B at moderate squeezing.", "The inspected validation protocol is empirical, not an interval certificate."),
    ("raw/chattopadhyay_brane_2609_37346v1.txt", "Sec.8.1, printed p.25", "This is not by itself a prediction of the reheating temperature: the actual value depends on the brane decay width, the initial abundance of the resonance, and the subsequent thermalization history.", "The recent thick-brane reheating estimate has explicit missing cosmological premises."),
]

def normalize(text):
    return " ".join(text.split())

def main():
    entries = []
    words_by_source = {}
    for relative, locator, quote, support in QUOTES:
        path = ROOT / relative
        body = path.read_bytes()
        if normalize(quote) not in normalize(body.decode("utf-8")):
            raise ValueError(f"quoted passage absent from retrieved text: {relative}, {locator}")
        count = len(quote.split())
        words_by_source[relative] = words_by_source.get(relative, 0) + count
        entries.append({"source_path": relative, "source_sha256": hashlib.sha256(body).hexdigest(), "locator": locator, "short_excerpt": quote, "quote_normalization": "Whitespace only; punctuation and document markup retained", "word_count": count, "supported_point": support})
    if any(n > 250 for n in words_by_source.values()):
        raise ValueError("short quotation budget exceeded")
    output = {"schema_version": 1, "status": "PASS_SHORT_PASSAGES_MATCH_RETRIEVED_TEXT", "scope": "A few supported passages from the selected review, not full-paper or API correctness verification.", "passages": entries, "quoted_words_by_source": words_by_source}
    (ROOT / "SELECTED_PASSAGES.json").write_text(json.dumps(output, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps({"status": output["status"], "quotes": len(entries), "maximum_quoted_words_per_source": max(words_by_source.values())}))

if __name__ == "__main__":
    main()
