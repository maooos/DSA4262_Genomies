"""Run the Task 2 workflow in order, using a frozen inventory by default."""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--stages', nargs='+', choices=['download', 'predict', 'analyse', 'validate', 'report'],
                        default=['download', 'predict', 'analyse', 'validate', 'report'])
    parser.add_argument('--refresh-inventory', action='store_true', help='Discover the current public inventory instead of the bundled snapshot')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--batch-size', type=int, default=200_000)
    args = parser.parse_args()
    if min(args.threads, args.workers, args.batch_size) < 1:
        parser.error('threads, workers and batch-size must be positive')
    folder = Path('data/sgnex')
    folder.mkdir(parents=True, exist_ok=True)
    if not (folder / 'manifest.csv').exists():
        for name in ['manifest.csv', 'inventory.json', 'directRNA_without_processed_json.csv']:
            target = folder / name
            if not target.exists():
                shutil.copy2(Path('reports/task2') / name, target)
    commands = {
        'download': ['scripts.sgnex_download', '--download', '--workers', str(args.workers),
                     *([] if args.refresh_inventory else ['--reuse-inventory'])],
        'predict': ['scripts.task2_predict', '--threads', str(args.threads), '--batch-size', str(args.batch_size)],
        'analyse': ['scripts.task2_analyse'],
        'validate': ['scripts.task2_validate'],
        'report': ['scripts.task2_report', '--publish-directory', 'task2_site/dist/report'],
    }
    for stage in ['download', 'predict', 'analyse', 'validate', 'report']:
        if stage in args.stages:
            print(f'\nTask 2 / {stage}', flush=True)
            subprocess.run([sys.executable, '-m', *commands[stage]], check=True)
            if stage == 'analyse':
                subprocess.run([sys.executable, '-m', 'scripts.task2_bundle_queries'], check=True)
    print('\nSelected stages complete. Rebuild the dashboard image to display regenerated results.', flush=True)


if __name__ == '__main__':
    main()
