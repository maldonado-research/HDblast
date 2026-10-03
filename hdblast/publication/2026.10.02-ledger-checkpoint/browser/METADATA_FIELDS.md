# Copy these fields into the existing Zenodo draft

Open **https://zenodo.org/uploads/23114217** (or https://zenodo.org/deposit/23114217), sign in as the owner, and edit this existing main-family draft. It belongs to concept **10.5281/zenodo.17088132** and currently has ten inherited files. The preserved published v24 record is **22347452**. Record **23112891** is removed. Keep companion **23111008**, concept **22922927**, unchanged.

The corrected API metadata attempt returned HTTP500, and its authenticated readback still has the old metadata and ten files. No additions or publication were attempted. This browser packet is the next concrete way to complete the same draft.

1. Preserve the existing draft and all ten inherited files. Use its existing DOI reservation; do not enter an invented DOI or create another record/version.
2. Enter every field below. [metadata.txt](metadata.txt) contains all values together; [METADATA_REQUEST.json](../METADATA_REQUEST.json) is the exact approved machine-readable payload. The displayed version number assigned by Zenodo is separate from the semantic Version field below.
3. Save the metadata and leave the draft unpublished. The fresh authenticated GET readback must match the pinned approved metadata and pass `PASS_EXACT_SAVED_DRAFT_METADATA` before uploading additions; matching visible text alone is insufficient. If Zenodo rewrites rich-text HTML, stop before uploads until that serialization difference has a reviewed metadata-format fix and the fresh guard check passes. Literal HTML tags must not appear as visible prose.
4. After the metadata check passes, upload exactly the twelve files in [ZENODO_UPLOAD_FILES.md](ZENODO_UPLOAD_FILES.md). Preserve all ten inherited files so the completed draft contains exactly 22 files.
5. Save again and leave the draft unpublished until its complete metadata and all 22 filename/byte/checksum pins are verified. Publication success requires the live public record as well as the authenticated saved record.

| Form field | Exact value |
|---|---|
| Title | HDBLAST: Registered Metric-Response Failures and Source-Free/Active-Source Ledger Diagnostics |
| Version | `2026.10.02-ledger-checkpoint` |
| Publication date | `2026-10-02` |
| Resource type | Software |
| Access | Open |
| License | Creative Commons Attribution 4.0 International (CC BY 4.0) |
| Language | English |
| Creator | Maldonado, Ricardo |
| Creator given / family names | Ricardo / Maldonado |
| Affiliation | Independent Researcher |
| ORCID | `0009-0009-3937-6527` |

Keep the existing creator; do not add a duplicate. Reports/data retain CC BY 4.0, bundled code retains its included MIT license, and third-party material retains its included licenses, as stated in the description.

For Description and Notes/Additional notes, use the exact [DESCRIPTION.html](DESCRIPTION.html) and [NOTES.html](NOTES.html) content. If the rich-text editor offers HTML/source mode, paste the corresponding file contents there. Otherwise open each file in a browser and copy its rendered paragraphs into the matching editor. The plain-text versions are also in metadata.txt. All scientific limitations and historical failures must remain present.

Enter all **37 keywords** from [KEYWORDS.txt](KEYWORDS.txt), one keyword per entry, and retain all **5 references** from [REFERENCES.txt](REFERENCES.txt). Add or preserve every related identifier below; retain the exact relationship and resource type. [RELATED_IDENTIFIERS.tsv](RELATED_IDENTIFIERS.tsv) is a spreadsheet-friendly copy.

| Identifier | Relation in form | Resource type | Scheme |
|---|---|---|---|
| 10.5281/zenodo.17069900 | Is supplement to (`isSupplementTo`) | Dataset (`dataset`) | DOI |
| 10.5281/zenodo.16937520 | Is supplement to (`isSupplementTo`) | Dataset (`dataset`) | DOI |
| 10.5281/zenodo.16907982 | Is supplement to (`isSupplementTo`) | Dataset (`dataset`) | DOI |
| 10.5281/zenodo.17088133 | Is supplement to (`isSupplementTo`) | Dataset (`dataset`) | DOI |
| 10.5281/zenodo.17547897 | Is supplement to (`isSupplementTo`) | Dataset (`dataset`) | DOI |
| https://arxiv.org/abs/2401.08437v1 | References (`references`) | Preprint (`publication-preprint`) | URL |
| 10.5281/zenodo.22347452 | References (`references`) | Software (`software`) | DOI |
| 10.5281/zenodo.22922927 | References (`references`) | Publication (`publication`) | DOI |
| https://github.com/maldonado-research/HDblast/commit/a8394d0127e58200cd4b9a87c14ae63a9dd69f02 | References (`references`) | Software (`software`) | URL |
| https://github.com/maldonado-research/HDblast/pull/20 | References (`references`) | Software (`software`) | URL |
| https://github.com/maldonado-research/HDblast/commit/06984aa6b142499c592850e6b4afd49c28ced394 | References (`references`) | Software (`software`) | URL |
| https://github.com/maldonado-research/HDblast/pull/18 | References (`references`) | Software (`software`) | URL |
| https://github.com/maldonado-research/HDblast/commit/79342726f097affa0c38d7726105d109f9b52a5c | References (`references`) | Software (`software`) | URL |
| https://github.com/maldonado-research/HDblast/commit/78f2ec797b5c3e4d1073f297b60ff946ad8bde71 | References (`references`) | Software (`software`) | URL |

The prior V25 overview filename is retained as a historical nine-addition snapshot. Start with HDBLAST_ACTIVE_SOURCE_LEDGER_REPORT_20261003.md for the completed active-source round. All date-stamped filenames remain unchanged for checksum reproduction.

This public packet contains exact metadata and source links, not twelve hosted release downloads. GitHub binary upload failed with HTTP400 before any asset was saved; the empty unpublished release draft was removed. Use the preserved public scientific-source/package instructions in ZENODO_UPLOAD_FILES.md. Full prepared files that require multipart reassembly remain local until reconstructed and verified.
