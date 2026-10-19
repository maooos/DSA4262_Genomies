FROM nginx:1.28-alpine@sha256:a8b39bd9cf0f83869a2162827a0caf6137ddf759d50a171451b335cecc87d236 AS dashboard
COPY docker/nginx.conf /etc/nginx/conf.d/default.conf
COPY task2_site/dist/ /usr/share/nginx/html/
EXPOSE 8080
HEALTHCHECK --interval=15s --timeout=3s --start-period=5s \
  CMD wget -q -O /dev/null http://127.0.0.1:8080/healthz || exit 1

FROM python:3.12-slim@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1 AS queries
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
WORKDIR /work
COPY requirements-dashboard.txt ./
RUN pip install --no-cache-dir -r requirements-dashboard.txt
COPY src/task2_query.py src/task2_query.py
COPY scripts/task2_serve.py scripts/task2_serve.py
COPY task2_site/query_data/ task2_site/query_data/
COPY task2_site/dist/data/summary.json task2_site/dist/data/summary.json
USER 65534:65534
EXPOSE 8081
HEALTHCHECK --interval=15s --timeout=3s --start-period=10s \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8081/healthz', timeout=2)" || exit 1
CMD ["python", "-m", "scripts.task2_serve", "--host", "0.0.0.0", "--port", "8081"]

FROM python:3.12-slim@sha256:a6e34c598f2467ed0e9a8d349809fcd8b5c603269512df273a0bb1784edc11b1 AS analysis
ENV PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/matplotlib
WORKDIR /work
COPY requirements.txt requirements-task2.txt requirements-dashboard.txt ./
RUN pip install --no-cache-dir -r requirements-task2.txt
COPY src/ src/
COPY scripts/ scripts/
COPY models/ models/
COPY examples/ examples/
COPY reports/task2/manifest.csv reports/task2/inventory.json reports/task2/directRNA_without_processed_json.csv reports/task2/
ENTRYPOINT ["python"]
CMD ["-m", "scripts.task2_pipeline", "--help"]
