# Minimal container for MPPH. Build: docker build -t mpph .
# Run:   docker run --rm -v "$PWD:/work" -w /work mpph Prochlorococcus --cluster
FROM python:3.11-slim

WORKDIR /app
COPY pyproject.toml README.md ./
COPY mpph ./mpph
RUN pip install --no-cache-dir .

# KEGG data is fetched at runtime and is NOT bundled; see DATA_SOURCES.md.
WORKDIR /work
ENTRYPOINT ["mpph"]
CMD ["--help"]
