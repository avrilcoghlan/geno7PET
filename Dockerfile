
FROM ghcr.io/astral-sh/uv:python3.14-bookworm-slim AS builder

WORKDIR /build

COPY src ./
COPY uv.lock pyproject.toml README.md ./

RUN uv build --wheel

FROM ghcr.io/astral-sh/uv:python3.14-bookworm-slim AS production

RUN apt-get update && \
    apt-get install -y --no-install-recommends ncbi-blast+ && \
    rm -rf /var/lib/apt/lists/*

# Copy the built wheel from the builder stage
COPY --from=builder /build/dist/*.whl /tmp/wheel/

RUN uv pip install --system --no-cache-dir /tmp/wheel/*.whl && \
    rm -rf /tmp/wheel

WORKDIR /app

ENTRYPOINT ["Geno7PET"]
