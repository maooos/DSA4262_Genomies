"""Standard-library integration checks; never load or execute the RNA model."""

import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from scripts.export_submission import export_submission
from scripts.run_pipeline import load_manifest, run_pipeline


class PipelineTests(unittest.TestCase):
    def test_manifest_rejects_missing_or_extra_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "reads.json").write_text("placeholder")
            manifest = root / "samples.csv"
            for row in ["sample1", "sample1,reads.json,unexpected"]:
                with self.subTest(row=row):
                    manifest.write_text(f"sample_id,input_path\n{row}\n")
                    with self.assertRaisesRegex(ValueError, "exactly two fields"):
                        load_manifest(manifest)

    def test_failed_sample_is_recorded_and_later_samples_continue(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, model, manifest = root / "reads.json", root / "model.joblib", root / "samples.csv"
            source.write_text("placeholder")
            model.write_text("placeholder")
            manifest.write_text("sample_id,input_path\nbroken,reads.json\nworking,reads.json\n")

            def fake_predict(command, **kwargs):
                output = Path(command[command.index("--output") + 1])
                if output.parent.name.startswith(".broken.work-"):
                    raise subprocess.CalledProcessError(1, command)
                output.write_text("transcript_id,transcript_position,score,prediction\nT1,2,0.5,1\n")
                output.with_suffix(".metadata.json").write_text(json.dumps({"site_count": 1}))

            with patch("scripts.run_pipeline.subprocess.run", side_effect=fake_predict) as process:
                self.assertEqual(run_pipeline(manifest, model, root / "out"), 1)
                self.assertEqual(process.call_count, 2)
            failed = root / "out/broken"
            self.assertTrue((failed / "failed.json").is_file())
            self.assertFalse((failed / "complete.json").exists())
            self.assertTrue((root / "out/working/complete.json").is_file())
            with self.assertRaises(FileExistsError):
                run_pipeline(manifest, model, root / "out", resume=True)
            def retry_predict(command, **kwargs):
                output = Path(command[command.index("--output") + 1])
                output.write_text("transcript_id,transcript_position,score,prediction\nT1,2,0.5,1\n")
                output.with_suffix(".metadata.json").write_text("{}")
            with patch("scripts.run_pipeline.subprocess.run", side_effect=retry_predict):
                self.assertEqual(run_pipeline(manifest, model, root / "out", retry_failed=True), 0)
            self.assertTrue((failed / "complete.json").is_file())

    def test_export_preserves_scores_and_rejects_duplicates(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, target = root / "scores.csv", root / "submission.csv"
            source.write_text("transcript_id,transcript_position,score,prediction\nT1,2,0.12345678901234567,0\n")
            self.assertEqual(export_submission(source, target), 1)
            self.assertEqual(target.read_text(), "transcript_id,transcript_position,score\nT1,2,0.12345678901234567\n")
            with self.assertRaises(FileExistsError):
                export_submission(source, target)
            source.write_text("transcript_id,transcript_position,score\nT1,2,0.1\nT1,2,0.2\n")
            with self.assertRaises(ValueError):
                export_submission(source, root / "bad.csv")
            self.assertFalse((root / "bad.csv").exists())

    def test_invalid_scores_never_publish(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for score in ["nan", "inf", "-0.1", "1.1"]:
                source = root / "scores.csv"
                source.write_text(f"transcript_id,transcript_position,score\nT1,2,{score}\n")
                with self.assertRaises(ValueError):
                    export_submission(source, root / "bad.csv")
                self.assertFalse((root / "bad.csv").exists())

    def test_manifest_and_dry_run_write_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "reads.json").write_text("placeholder")
            (root / "model.joblib").write_text("placeholder")
            manifest = root / "samples.csv"
            manifest.write_text("sample_id,input_path\nsample1,reads.json\n")
            self.assertEqual(load_manifest(manifest)[0][1], root / "reads.json")
            with patch("scripts.run_pipeline.subprocess.run") as process:
                self.assertEqual(run_pipeline(manifest, root / "model.joblib", root / "out", dry_run=True), 0)
                process.assert_not_called()
            self.assertFalse((root / "out").exists())
            manifest.write_text("sample_id,input_path\n../unsafe,reads.json\n")
            with self.assertRaises(ValueError):
                load_manifest(manifest)

    def test_mocked_prediction_resume_and_changed_input(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source, model, manifest = root / "reads.json", root / "model.joblib", root / "samples.csv"
            source.write_text("placeholder")
            model.write_text("placeholder")
            manifest.write_text("sample_id,input_path\nsample1,reads.json\n")

            def fake_predict(command, **kwargs):
                output = Path(command[command.index("--output") + 1])
                output.write_text("transcript_id,transcript_position,score,prediction\nT1,2,0.5,1\n")
                output.with_suffix(".metadata.json").write_text(json.dumps({"site_count": 1}))

            with patch("scripts.run_pipeline.subprocess.run", side_effect=fake_predict) as process:
                self.assertEqual(run_pipeline(manifest, model, root / "out"), 0)
                self.assertEqual(run_pipeline(manifest, model, root / "out", resume=True), 0)
                self.assertEqual(process.call_count, 1)
                source.write_text("changed")
                with self.assertRaises(FileExistsError):
                    run_pipeline(manifest, model, root / "out", resume=True)

    def test_output_lock_rejects_concurrent_run(self):
        from scripts.run_pipeline import output_lock
        with tempfile.TemporaryDirectory() as directory:
            with output_lock(Path(directory)):
                with self.assertRaisesRegex(RuntimeError, "already holds"):
                    with output_lock(Path(directory)):
                        pass


if __name__ == "__main__":
    unittest.main()
