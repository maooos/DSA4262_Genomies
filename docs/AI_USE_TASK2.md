# AI use for Task 2

Draft record for the assignment; the team should review and combine it with
its records for Task 1 and any other AI assistance.

| AI tool used | Prompt and output | How the output was used |
| --- | --- | --- |
| OpenAI Codex, 12 October 2026 | Asked for Task 2 ideas given an existing model, then “Can you do it”. Produced the SG-NEx inventory/downloader, batch inference wrapper, descriptive analysis, report generator, tests, and interactive explorer. | Implemented in this repository; H04 was reused without retraining. Tests and actual SG-NEx runs were used to check the implementation. |
| OpenAI Codex, same session | Examined the handout, local inference code, SG-NEx documentation, live S3 listing and sample metadata; produced scope notes and a results-section draft. | Sources and limitations are included in the guide, report and explorer. The team must review interpretation and adapt the draft to the combined report's page limit. |
| OpenAI Codex, 19 October 2026 | “Help me do task 2, using docker … make the interactive visualisation for researchers easy to use and easy to interpret/understand.” Recovered the earlier analysis code, added Docker services, rebuilt the dashboard and added verification. | Existing predictions were retained; all 23 inputs were verified inside Docker, and exported summaries/transcripts were checked against the prediction files. UI text separates model scores, coverage and replicate support. |
| OpenAI Codex, 19 October 2026 | “Can you make it like user can enter any number of reads?” Added a validated number field and local coverage-query service. | All saved predictions are bundled for exact custom filters, including 1–9 reads; scientific comparisons are checked against the original analysis functions. |

The October 12 entries below were recovered from the repository's earlier AI-use
record, rather than independently reconstructed in the October 19 session.

## Issues found in the Docker/dashboard session

1. The first nginx configuration served `.mjs` files as `application/octet-stream`,
   so the browser did not execute the imported module. HTTP header inspection
   identified the cause. An explicit JavaScript MIME type fixed it, and the
   dashboard loaded successfully in Chrome.
2. The recovered prediction wrapper assumed that an input's absolute path never
   changes. Docker relocates it under `/work`, invalidating that assumption.
   SHA256 verification against original download audits now permits relocation.
   Tests reject missing audits, wrong sample identities and same-size corruption.

## Assumptions corrected during verification

1. **The metadata protocol was spelled `direct-RNA`.** An initial exploratory
   query used that spelling and incorrectly returned no direct-RNA samples.
   Inspecting the protocol column's value counts showed that the actual value
   is `directRNA` (56 rows). The inventory filter and reconciliation use the
   verified spelling, and 33 metadata runs without processed JSON are recorded
   separately. The older linked processed-sample list was also checked against
   the live S3 listing rather than used as a complete inventory.
2. **Rounding scores for a smaller browser export was harmless.** Inspection
   found that rounding can move a score across the saved decision threshold.
   The exporter now preserves full-precision scores. A regression test checks
   scores immediately below, exactly at, and immediately above the threshold.
3. **Transcript IDs were directly comparable across every sample.** The first
   full comparison showed zero MCF7 replicate overlap. Inspecting the original
   identifiers revealed versioned versus unversioned Ensembl names. Their
   969,788 shared raw coordinate pairs all had identical 7-mers after name
   normalization. The analysis now preserves original IDs, harmonizes only
   unambiguous version suffixes, and still requires coordinate/sequence matches.

Additional checks cover missing observations, sequence matching, low coverage,
technical runs that belong to the same replicate, and preservation of original
H04 scores. These checks verify computation; they do not validate the biological
truth of the predicted modifications.
