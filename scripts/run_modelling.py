"""Run the pooled modelling notebook's Python cells from the command line.

Uses the active Python environment, with no extra Jupyter execution dependency.
The source notebook stays unexecuted; tables/models and plot PNGs go to the new
run directory. --prepare-only stops after saving X/y/metadata and the split.
"""

import argparse
import json
import os
from pathlib import Path


def run_notebook(path, *, prepare_only=False, namespace=None):
    notebook = json.loads(Path(path).read_text())
    stop_tag = "prepared-inputs"
    if prepare_only and not any(stop_tag in c.get("metadata", {}).get("tags", []) for c in notebook["cells"]):
        raise ValueError("Notebook has no prepared-inputs boundary; refusing to run preparation.")
    namespace = {} if namespace is None else namespace
    namespace.setdefault("__name__", "__main__")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    previous_show = plt.show
    plot_number = 0

    def save_plots(*args, **kwargs):
        nonlocal plot_number
        for figure_number in plt.get_fignums():
            figure = plt.figure(figure_number)
            if namespace.get("SAVE_OUTPUTS") and "OUTPUT_DIR" in namespace:
                plot_number += 1
                figure.savefig(Path(namespace["OUTPUT_DIR"]) / f"plot_{plot_number:02d}.png", dpi=140, bbox_inches="tight")
            plt.close(figure)

    plt.show = save_plots
    try:
        for number, cell in enumerate(notebook["cells"], 1):
            if cell["cell_type"] != "code":
                continue
            print(f"[notebook] cell {number}/{len(notebook['cells'])}", flush=True)
            exec(compile("".join(cell["source"]), f"{path}:cell_{number}", "exec"), namespace)
            if prepare_only and stop_tag in cell.get("metadata", {}).get("tags", []):
                print("Preparation complete; no model fitting was run.")
                break
    finally:
        plt.show = previous_show
        plt.close("all")
    return namespace


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-only", action="store_true", help="Save split and features, then stop before CV/training")
    parser.add_argument("--processed-dir", type=Path, default=Path("data/processed"))
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("--threads must be positive")
    root = Path(__file__).resolve().parents[1]
    processed_dir = args.processed_dir.resolve()
    previous_cwd = Path.cwd()
    try:
        os.chdir(root)
        run_notebook(root / "notebooks/data0_modelling_pipeline.ipynb",
                     prepare_only=args.prepare_only,
                     namespace={"PROCESSED_DIR": processed_dir, "MODEL_THREADS": args.threads})
    finally:
        os.chdir(previous_cwd)


if __name__ == "__main__":
    main()
