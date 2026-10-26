# SG-NEx Atlas: a guide for first-time users

SG-NEx Atlas is an interactive website for exploring **predicted m6A RNA modifications**. It helps you move from a broad question, such as “Do these cell lines show similar patterns?”, to the individual sequencing observations behind a prediction.

The application uses a saved machine-learning model called **H04**. You can compare its predictions, inspect support from repeated samples, and download results for further exploration.

**The most important distinction:** a prediction is evidence from a model. A score of 0.8 does not mean that 80% of the RNA is modified, and a missing observation does not mean that a site is unmodified.

## Contents

- [Open the application](#open-the-application)
- [A few useful terms](#a-few-useful-terms)
- [Take a five-minute tour](#take-a-five-minute-tour)
- [Controls used across the application](#controls-used-across-the-application)
- [Study overview](#study-overview)
- [Compare cell lines](#compare-cell-lines)
- [Explore a transcript](#explore-a-transcript)
- [How to read the data](#how-to-read-the-data)
- [Download and share your findings](#download-and-share-your-findings)
- [Common questions and problems](#common-questions-and-problems)

## Open the application

If a teammate has already started the application for you, open the address they provide.

To run it on your own computer, follow the [setup instructions in the README](../README.md#start-the-website). You need Docker and both complete dashboard data folders from the same analysis run:

```text
task2_site/dist/data/
task2_site/query_data/
```

These large data folders are shared separately from GitHub. Cloning the repository alone does not include them.

Once the folders are in place and Docker is running, open a terminal in the project folder and run:

```bash
docker compose up --build -d
```

Then open **http://localhost:8080**. Use the running website instead of double-clicking `index.html`, because the application needs to load data through a web server.

The bundled analysis covers **23 sequencing runs across seven cell lines**: A549, H9, HEYA8, Hct116, HepG2, K562, and MCF7. This is the available processed collection used by this project, rather than every raw SG-NEx experiment. The badge at the top shows how many expected runs are included.

## A few useful terms

| Term | Plain-language meaning |
| --- | --- |
| m6A | A chemical modification of RNA that this project tries to predict. |
| Cell line | A named group of cells being studied, such as A549 or HepG2. |
| Transcript | An RNA sequence identified here by an Ensembl ID such as `ENST00000451270`. |
| Position | A location along a transcript. These numbers are transcript coordinates, not chromosome coordinates. |
| Site | A transcript ID, a position, and the surrounding seven-base sequence, considered together. |
| 7-mer | The seven-base sequence around a site. It helps ensure that comparisons refer to matching sequence contexts. |
| Read support / coverage | How many sequencing reads cover a site within one run. |
| Sequencing run | One sequencing measurement. One replicate group may contain several technical runs. |
| Replicate group | A named repeat-sample group in the dataset. Its qualifying technical runs are averaged before comparing groups. |
| Model score | A number from 0 to 1. Higher scores indicate stronger model evidence for modification. |
| Prediction threshold | The score cutoff for calling a prediction positive. A score equal to the threshold is positive. |
| Eligible observation | A site in a run that meets the selected minimum read support. |

A site observed in three runs contributes **three run–site observations**. Therefore, an observation count is not necessarily a count of distinct sites.

## Take a five-minute tour

1. Open **Study overview** and leave **Minimum read support** at **20 reads**.
2. Look at the dots in **How often is a site called modified?** Hover over or click a dot to inspect a run.
3. Select a square in **Which cell lines share predictions?** to open that pair in **Compare cell lines**. Alternatively, choose A549 and HepG2 from the comparison dropdowns.
4. Scroll to **Sites with consistent differences**. Click a transcript link in a candidate row to inspect that exact position.
5. In **Explore a transcript**, click the score overview to select a position, then inspect its linked heatmap column and detail table. Select **Expand tracks** for individual run scores/read support. You can also enter `ENST00000451270` directly or use a **Try a candidate** button.
6. Change the read support to **50 reads** and see whether the pattern still has support. A site disappearing means it no longer meets the filters; it does not prove that the prediction was wrong.
7. Download the relevant CSV and use **How to read the data** to check the interpretation rules before describing a finding.

## Controls used across the application

### Minimum read support

This controls which observations are included. It is applied to each run separately, before scores are averaged.

- **20 reads** is the default.
- **10 reads** admits more observations with lower read support.
- **50 reads** is stricter and generally retains fewer observations.
- To choose another cutoff, enter a positive whole number and select **Apply**, or press Enter.

For example, at a cutoff of 20, an observation with 18 reads is excluded and one with 25 reads is included. Counts from separate runs are not added together to help a low-coverage observation pass this filter.

Typing a custom value alone does not apply it. Invalid values show an error and leave the previous filter in place. The preset buttons apply immediately.

The 10/20/50 settings use prepared results. Other settings calculate summaries from saved predictions through the local query service. Very low cutoffs can take a minute or more. Changing this control does not retrain H04 or generate new model scores.

### Reset view

Select **Reset view** to restore 20 reads, the saved prediction threshold, and the default A549/HepG2 pair. It also clears candidate filters and transcript selection, and restores the transcript cell-line filter to all cell lines. You remain on the current page.

### Navigation and reading help

Use the left-hand navigation to switch between the four pages. The **Reading guide** button and the footer's **Methods & interpretation** link open **How to read the data**.

The application supports keyboard navigation. Candidate transcript links provide a keyboard-accessible way to explore sites without clicking scatter-plot points.

## Study overview

Use this page to understand the collection before looking at individual sites.

### Summary cards

| Card | What it tells you |
| --- | --- |
| Cell lines | How many cell lines are included. |
| Sequencing runs | How many runs are included, keeping technical runs separate. |
| Eligible observations | Total run–site observations passing the selected read-support filter. |
| Transcript catalogue | Number of transcripts in the prepared lookup catalogue, which contains sites exported at at least 10 reads. This is not a count recalculated for every filter. |

### How often is a site called modified?

Each dot represents one sequencing run. The horizontal line connects the lowest and highest predicted-positive fractions within that cell line. Hover over a dot for details, or click it for a short run summary.

The plotted fraction is:

```text
sites predicted positive in the run / sites meeting read support in the run
```

For example, 200 positive sites out of 1,000 eligible sites gives 20%. This is a percentage of eligible sites called positive by the model, not a measured percentage of modified RNA in the cells.

The **Prediction threshold** dropdown offers 0.15, the saved H04 threshold (approximately **0.2084**), 0.30, and 0.50. Raising it makes the positive-call rule stricter.

**This dropdown changes only the run-level summaries and their CSV export.** The overlap heatmap, replicate agreement, cell-line comparisons, and transcript calls continue to use the saved H04 threshold.

Select **Download CSV** to save the run summaries for the current read support and overview threshold.

### Which cell lines share predictions?

This heatmap compares pairs of cell lines on sites eligible in both. Select a square to open that pair in **Compare cell lines**.

The number is **Jaccard similarity**:

```text
sites positive in both groups / sites positive in either group
```

For example, 30 shared positive sites out of 60 sites positive in either group gives 0.5. A value of 1 means complete positive-set overlap; 0 means no positive overlap when at least one positive site exists. If neither group has positive sites, the value is undefined and shown as a dash.

Each pair can have a different set of shared eligible sites. Read the overlap alongside the number of sites supporting that comparison.

### Do repeat samples tell a similar story?

This table compares named replicate groups within each cell line. It shows the groups compared, shared eligible sites, score-ranking agreement, and positive-site overlap.

**Spearman correlation** describes whether the groups rank sites similarly: values near 1 indicate similar rankings, near 0 indicate little ranking agreement, and near -1 indicate opposite rankings. A dash means the statistic could not be calculated, for example because scores are constant or there are too few observations.

Technical runs are averaged within each replicate group first. Agreement is useful supporting evidence, but shared sequence features can also produce agreement; it is not independent proof of biological accuracy.

## Compare cell lines

Use this page to compare matching sites in two cell lines and find candidates for follow-up.

### Choose a pair

Choose different values for **Cell line A** and **Cell line B**. The **swap button** exchanges them, including the plot axes and the sign of the score difference.

A site enters a comparison only when its transcript ID, position, and seven-base sequence match, and each cell line has at least **two qualifying replicate groups** at that site.

The four cards show the total shared eligible sites, sites positive in both cell lines, positive-site overlap, and number of consistent candidates. These counts use the full eligible comparison.

### Read the scatter plot

Each point represents a matched site:

- The horizontal axis is the mean score in cell line A.
- The vertical axis is the mean score in cell line B.
- Points near the diagonal have similar scores in both.
- Points above the diagonal have higher scores in B; points below it have higher scores in A.
- The horizontal and vertical threshold lines divide sites into positive in both, positive only in A, positive only in B, or below threshold in both. The panel beside the plot gives the color key and full counts.

Hover over a point to see its transcript, position, and scores. Click it to open that position in the transcript explorer. The adjacent panel also shows score-ranking correlation.

The plot displays **up to 2,000 sampled sites** for readability. Its visible point count can be much smaller than the total shown in the cards.

### Find consistent candidates

The candidate table uses a stricter rule than simply comparing two means: **all observed qualifying replicate scores must be at or above the threshold in one cell line and below it in the other**, with at least two qualifying replicate groups per cell line.

For an invented example using a threshold of about 0.2084, replicate scores of 0.60 and 0.70 in A versus 0.10 and 0.15 in B satisfy this rule. If one qualifying replicate in A instead scores 0.10, the site no longer has a consistent difference.

Available controls:

| Control | What it does |
| --- | --- |
| Find a transcript | Filters by all or part of a transcript ID; matching ignores letter case. |
| Higher score in | Keeps candidates with higher scores in A, B, or either cell line. |
| Sort by | Orders by largest absolute score difference, transcript position, or transcript ID. |
| Previous / Next | Moves through the table, 20 candidates per page. |
| Transcript link | Opens the transcript and selected candidate position. |
| Download filtered CSV | Exports every row matching the current candidate filters within the available candidate subset, across all table pages. |

The table includes transcript ID, position, 7-mer, both cell-line mean scores, **Difference A - B**, and the number of qualifying replicates in each cell line. A positive difference means A has the higher score; a negative difference means B does.

**The table contains at most the top 500 candidates ranked by absolute score difference.** Search, sorting, and CSV download operate within that subset. The summary card may report more than 500 candidates. A transcript missing from this table may still be available through direct transcript lookup.

These are exploratory candidates. The application does not attach p-values or claim statistically significant differences.

## Explore a transcript

Use this page to inspect the observations behind a transcript's predictions.

### Find a transcript

Enter an Ensembl transcript ID and select **Explore transcript**, or press Enter. Suggestions appear as you type an ID prefix. The two **Try a candidate** buttons provide ready-made examples.

Known versioned IDs are supported. If the app resolves one to a shared identifier, it displays a message explaining which ID is shown. Conflicting versions are kept separate, and original identifiers are retained in the downloads.

Search uses transcript IDs. Gene symbols are not mapped in this dataset. Transcripts outside the prepared catalogue can be checked through the query service, so the suggestions are not a complete list of everything searchable.

### Choose cell lines

Use **Show runs from** to select one or more cell lines. **All cell lines** selects or clears all.

The selection and minimum read support apply to every transcript view and the CSV export. Changing filters may remove positions that no longer qualify.

### Transcript overview: predicted scores

The overview stays visible when tracks are collapsed. Bars show the highest score **among qualifying runs** at each position; the dark line marks the selected site. A high bar may reflect just one run, rather than agreement across replicates.

Click to select the nearest eligible position. The overview, heatmap and detail table update together. The dark line cannot be dragged.

Coordinate ticks cover the observed span, which is not necessarily the full transcript length. The summary counts all eligible positions under your current filters.

### Individual sequencing runs

Select **Expand tracks** to reveal per-run plots and navigation controls. Runs are labelled **Cell line · replicate / run** (Eg. **A549 · 5 / 1**):

| Panel | What it shows |
| --- | --- |
| Score | Run-level scores from 0 to 1. The orange dashed line marks the saved H04 threshold. |
| Reads | Read counts at scored sites, on a shared log(1 + reads) scale—not continuous sequencing coverage. |

Hover for exact values; click a site to update the linked views.

Use **Move**, **Zoom**, **Full observed span**, or **From/To** → **View range** to navigate. You can also pan by dragging the run plots or the **green overview rectangle**, and zoom with **Ctrl + wheel**. Keyboard shortcuts appear in the expanded section.

Collapsing preserves the viewing range. The viewing status counts positions within that range, while the overview summary counts all eligible positions.

**Zoom changes only the track display**. The heatmap, detail table and CSV continue to use all observations qualifying under your filters.

### Read the transcript summary and heatmap

The cards show observed positions, qualifying runs, eligible observations, and positive observations under the current transcript filters. A qualifying run has at least one eligible observation; a positive observation is a run-level score at or above the saved threshold.

In the heatmap, rows are sequencing runs and columns are observed transcript positions. Cell color represents the score, with darker green indicating higher scores. Hover over a populated cell to see its score and read count.

Use **Show runs from** to choose one cell line or all cell lines. Use **Positions** arrows to browse up to 24 observed positions at a time. These columns need not represent consecutive positions along the transcript. Selecting a site in the overview or run tracks opens its heatmap page automatically.

A dash means that an observation is missing or does not meet the current read-support filter. It is never treated as a negative prediction.

### Inspect a position

Select a position heading, populated heatmap cell, overview site or individual track site to update **Selected position**. The overview's dark line follows the same selection.

The cell-line summaries show the mean score, qualifying replicate count, number of positive replicates, and sequence context. One qualifying replicate is labeled as limited support; it can be inspected here but is insufficient for pairwise comparisons.

The detailed table shows each qualifying run's cell line, replicate/run identifiers, original transcript ID, 7-mer, read count, score, and prediction label.

The displayed cell-line mean is calculated in two steps:

1. Average qualifying technical-run scores within each replicate group.
2. Average those replicate scores with equal weight.

For example, if one replicate has two runs scoring 0.2 and 0.4, its score is 0.3. If another replicate scores 0.7, the cell-line mean is `(0.3 + 0.7) / 2 = 0.5`. The replicate with two technical runs does not get twice the weight. Different sequence contexts are kept separate.

### Download transcript observations

**Download run-level CSV** exports observations for the current transcript, read support, and cell-line filter. It includes all qualifying positions, not just the visible position page or selected position.

The export retains original transcript IDs, run identity, sequence, read counts, full stored scores, predictions, threshold, and model identifier. The screen rounds scores for readability.

## How to read the data

This page provides the built-in glossary, comparison rules, study scope, and interpretation notes. It also offers:

| Feature | What you get |
| --- | --- |
| Open results report | A separate report covering the data, predicted-positive fractions, agreement, overlap, candidate differences, sensitivity, and limitations. It is a saved report and does not change with your current dashboard filters. |
| Data source & access | The SG-NEx data-access documentation. The sidebar also links to the SG-NEx project. |
| Download study metadata | A JSON file describing the study bundle, including its saved summaries and provenance. |
| Model and data provenance | An expandable record of the analysis date, model identifier, threshold, expected/missing runs, identifier handling, and sequence conflicts. |

Keep these interpretation points in mind:

- Run-level positive fractions use different observed site sets, so they do not directly rank true m6A abundance between cell lines.
- Read depth, expression, batch effects, and shared sequence features can affect the patterns.
- Named replicate groups should not automatically be assumed to be experimentally independent.
- The available training metadata does not resolve possible overlap with SG-NEx samples. This application is not an independent accuracy benchmark.
- Different transcript isoforms may refer to the same genomic nucleotide. Transcript sites should not automatically be counted as unique genomic sites.

## Download and share your findings

| Download | Scope |
| --- | --- |
| Overview CSV | Run-level summaries at the current read support and overview threshold. |
| Candidate CSV | Filtered candidates across all pages of the available top-500 subset. |
| Transcript CSV | All qualifying observations for the current transcript and cell-line filter. |
| Study metadata JSON | The loaded study bundle; it is not a save file for every current interface setting. |

CSV files can be opened in spreadsheet software. They include the model identifier and relevant filter settings so the results are easier to trace.

You can also copy the browser address to share a view. The address records the page and read support, plus the overview threshold, comparison pair, or transcript and selected position where applicable.

Candidate search text, candidate direction/sort, table page, and the transcript cell-line filter are **not** saved in the address. Tell your friend about these settings separately.

An address beginning with `localhost` refers to the computer opening it. Your friend needs their own running copy with the same data, or access to a shared running deployment. Sending the link does not send the dataset.

When describing a finding, record the cell-line pair, transcript ID, position, read-support cutoff, threshold, and supporting replicate counts. This makes it easier for someone else to inspect the same evidence.

## Common questions and problems

| What you see | Explanation or next step |
| --- | --- |
| Data fails to load after cloning | Restore both complete data folders, then rebuild with `docker compose up --build -d`. See the [README](../README.md). |
| The page fails when opened as a local file | Start Docker and open `http://localhost:8080` instead of double-clicking the HTML file. |
| Presets work but a custom cutoff fails | Custom cutoffs need the query service. Start both services through Docker; check them with `docker compose ps`. |
| A custom cutoff is taking a long time | Low read cutoffs can include many more observations. Allow the calculation to finish; a preset is a useful starting point. |
| No observations meet the filters | Try a lower read cutoff or another cell line. This means no eligible observations, not proof of no modification. |
| A transcript cannot be found | Check the Ensembl transcript ID and avoid gene symbols. The transcript may be absent from this study. |
| A candidate disappears at 50 reads | It may have lost eligible observations, sufficient replicate support, or the consistent-difference pattern after filtering. |
| The candidate total is larger than the table/export | The total uses all consistent candidates; the table and export use the top-500 subset. |
| The scatter plot has fewer points than the site count | The plot samples at most 2,000 sites; the count uses the full eligible comparison. |
| Changing the overview threshold leaves comparisons unchanged | This is expected. Comparisons and transcript calls always use the saved H04 threshold. |
| Transcript data reports a decompression error | Use an up-to-date browser and keep the supplied `.json.gz` files intact. |
| The query service reports mismatched snapshots | Replace both data folders with copies from the same analysis run, then rebuild. |
| Individual run plots or zoom controls are not visible | Select **Expand tracks** below the always-visible transcript overview. |

For setup and Docker troubleshooting, see the [README](../README.md). For analysis methods, validation, and regeneration commands, see the [Task 2 team guide](TASK2.md). For scoring new input files with H04 outside the dashboard, see the [inference guide](INFERENCE.md).
