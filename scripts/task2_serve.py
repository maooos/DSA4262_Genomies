"""Serve the local dashboard and its exact custom-coverage API."""
import argparse
import json
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

from src.task2_query import CoverageQuery, parse_reads


class DashboardHandler(SimpleHTTPRequestHandler):
    extensions_map = {**SimpleHTTPRequestHandler.extensions_map, ".mjs": "application/javascript"}

    def do_GET(self):
        url = urlsplit(self.path)
        if url.path == "/healthz":
            return self.send_json({"status": "ok"})
        if not url.path.startswith("/api/"):
            return super().do_GET()
        try:
            params = parse_qs(url.query, keep_blank_values=True)
            if any(len(v) != 1 for v in params.values()):
                raise ValueError("Each query parameter must appear once.")
            params = {k: v[0] for k, v in params.items()}
            # Validate BEFORE queueing expensive work.
            reads = parse_reads(params.get("reads", "20"))
            with self.server.query_lock:
                engine = self.server.engine
                if url.path == "/api/overview":
                    result = engine.overview(reads)
                elif url.path == "/api/pair":
                    result = engine.pair(reads, params.get("a"), params.get("b"))
                elif url.path == "/api/transcript":
                    result = engine.transcript(params.get("id", ""))
                else:
                    return self.send_json({"error": "Unknown endpoint."}, 404)
            self.send_json(result)
        except ValueError as exc:
            self.send_json({"error": str(exc)}, 400)
        except (BrokenPipeError, ConnectionResetError):
            pass
        except Exception as exc:
            self.log_error("Query failed: %s", exc)
            self.send_json({"error": "Calculation failed. Retry or check the query service logs."}, 500)

    def send_json(self, value, status=200):
        content = json.dumps(value, allow_nan=False, separators=(",", ":")).encode()
        self.send_response(status)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(content)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.end_headers()
        self.wfile.write(content)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--data", type=Path, default=Path("task2_site/query_data"))
    parser.add_argument("--site", type=Path, default=Path("task2_site/dist"))
    args = parser.parse_args()
    engine = CoverageQuery(args.data, args.site / "data")
    server = ThreadingHTTPServer((args.host, args.port), partial(DashboardHandler, directory=str(args.site)))
    server.engine, server.query_lock = engine, threading.Lock()
    print(f"SG-NEx Atlas: http://{args.host}:{args.port}", flush=True)
    try:
        server.serve_forever()
    finally:
        server.server_close()
        engine.close()


if __name__ == "__main__":
    main()
