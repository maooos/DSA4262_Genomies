# SG-NEx Atlas

An interactive website for exploring H04-predicted m6A RNA modifications across
**23 SG-NEx sequencing runs and seven cell lines**.

**For teammates: clone the code, get the two dashboard data folders, then start
Docker.** The datasets are shared separately because they are too large to keep
in this Git repository. A fresh clone does not include them.

## Start the website

### 1. Get the code and Docker

- Clone or download this repository and open its folder in VS Code.
- Install and start [Docker Desktop](https://docs.docker.com/get-started/get-docker/).
- Open a terminal in the project folder: the one containing `Dockerfile`,
  `compose.yaml` and this `README.md`.

You do not need to install Python or download the original 33.4 GB dataset to
use the website with the shared dashboard data.

### 2. Put the dashboard data in these exact folders

Ask a teammate who already has the data for **both complete folders**, from the
same analysis run:

| Folder to copy | Where to put it in your project | Approximate size |
| --- | --- | --- |
| Dashboard `data` folder | `task2_site/dist/data/` | 112 MB |
| `query_data` folder | `task2_site/query_data/` | 340 MB |

After copying, your project should look like this:

```text
DSA4262_Genomies/                 ← your project folder
├── README.md
├── Dockerfile
├── compose.yaml
└── task2_site/
    ├── dist/
    │   ├── index.html           ← website code, already in the repository
    │   ├── app.js
    │   └── data/               ← copy the dashboard data here
    │       ├── summary.json
    │       ├── transcript_index.json
    │       ├── transcript_aliases.json
    │       ├── pairs/          ← copy the whole folder
    │       └── transcripts/    ← copy the whole folder, including .json.gz files
    └── query_data/             ← copy the query data here
        ├── manifest.json
        ├── identifiers.parquet
        ├── SGNex_A549_directRNA_replicate5_run1.parquet
        └── ...                ← all 23 sample .parquet files
```

**Copy the whole folders, not just the JSON files open in VS Code.** Avoid an
extra folder level such as `task2_site/dist/data/data/`. The dashboard `data/`
folder belongs under `task2_site/dist/`, not at the project root.

If you receive a `task2-data.tar.gz` archive prepared using the instructions below,
place it in the project folder and extract it there:

```bash
tar -xzf task2-data.tar.gz
```

This puts both folders in the correct locations. Keep the `.json.gz` files inside
`transcripts/` compressed; the website reads them directly.

### 3. Start Docker and open the website

From the project folder, run:

```bash
docker compose up --build -d
```

Open **http://localhost:8080** in your browser.

The first build needs internet to download dependencies. Once built, the website
uses local data. Docker starts both the website and its custom-filter query
service automatically.

## What you can do on the website

New to the application? Read the [first-time user guide](docs/USER_GUIDE.md) for a
five-minute tour, explanations of every page and control, and help interpreting
and sharing your results.

| Page or control | Use it to… |
| --- | --- |
| Study overview | See per-run predictions, cell-line overlap and replicate agreement. |
| Compare cell lines | Find matched sites with different predicted scores in two cell lines. |
| Explore a transcript | Inspect scores, read counts and missing observations at individual positions. |
| How to read the data | Understand the methods, limitations and results report. |
| Minimum read support | Enter a whole number ≥1 and click **Apply**, or use the 10/20/50 presets. |

Custom read filters recalculate from saved predictions; they do not retrain H04.
Very low cutoffs can take a minute or more. A model score is not a measured
fraction of modified RNA, and a missing observation does not mean unmodified.

## Everyday Docker commands

Run these from the project folder:

```bash
docker compose up -d                      # Start again after the first build
docker compose down                       # Stop the website
docker compose up -d                      # Restart after a crash or computer sleep
docker compose ps                         # Check whether services are healthy
docker compose logs dashboard queries     # Read startup or calculation errors
docker compose up --build -d              # Rebuild after code or data changes
```

| Problem | What to check |
| --- | --- |
| Docker cannot connect to its engine | Open Docker Desktop and wait until it is running. |
| Build fails because a data file or folder is missing | Check the folder tree above, especially `dist/data/summary.json` and `query_data/manifest.json`. |
| Query service reports that data snapshots differ | Copy both complete dashboard folders from the same analysis run, then rebuild. |
| Port 8080 is already in use | Create a `.env` file in the project folder containing `PORT=8081`, run the startup command again, then open http://localhost:8081. |
| Website still shows older results | Rebuild with `docker compose up --build -d`, then refresh your browser. |

## Share the dashboard data with a teammate

On a machine that already has both data folders, run this from the project folder:

```bash
tar -czf task2-data.tar.gz task2_site/dist/data task2_site/query_data
```

Transfer the archive through your team's shared storage or a direct file transfer.
Your teammate extracts it using the command in step 2. The two folders total about
452 MB before archiving. The archive and dataset folders are ignored by Git;
commit the code, Docker setup, model and documentation to GitHub.

## Optional: download or rerun the full analysis

**Skip this section if you only want to use the website and already have the two
dashboard data folders.**

If nobody has the dashboard data, the pipeline can download the processed SG-NEx
inputs, run H04 and generate both folders automatically:

```bash
docker compose run --rm --build analysis -m scripts.task2_pipeline
docker compose up --build -d
```

This downloads about **33.4 GB**. Allow at least **45 GB of free disk space** and
start with **8 GB or more of Docker memory**. Full analysis may need more memory.

If transferring existing research inputs or outputs instead, use these locations:

| Data you have | Exact location in the project | Needed just to open the website? |
| --- | --- | --- |
| Processed SG-NEx input for a sample | `data/sgnex/raw/<sample_id>/data.json` | No |
| That sample's download audit | `data/sgnex/raw/<sample_id>/download.json` | No; keep it with the input when transferring. |
| SG-NEx manifest and inventory | `data/sgnex/manifest.csv` and `data/sgnex/inventory.json` | No |
| Saved full prediction outputs | `outputs/task2/predictions/` | No; keep each sample's CSV, Parquet and metadata JSON together. |
| Generated analysis tables and figures | `outputs/task2/analysis/` | No |
| Original data0 input used by the inference notebook | `data/raw/data0/dataset0.json.gz` | No |

For example, one SG-NEx input belongs at
`data/sgnex/raw/SGNex_A549_directRNA_replicate5_run1/data.json`.
Keep the sample folder names unchanged. The pipeline creates the required folders
when downloading or generating results; these local data and output folders are
excluded from Git.

See the [Task 2 team guide](docs/TASK2.md) for rerunning selected stages,
validation and the comparison rules. See the [H04 inference guide](docs/INFERENCE.md)
for using the notebook or scoring your own input files.
