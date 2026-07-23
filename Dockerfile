# grapetree-rs — byte-identical MSTreeV2 with no `edmonds` binary.
#
# The default (bit-identical) mode delegates only the two non-portable numerics
# (harmonic weights' float32 SIMD sum, and branch_recraft's np.log) to NumPy via
# the tiny shims in shim/. Everything else — including the minimum spanning
# arborescence, a pure-Rust port of the `edmonds` binary — is Rust. So the runtime
# needs just the static binary + Python 3 + NumPy.
#
# For distance / MSTree / MSTreeV2 this image is a byte-identical drop-in. The NJ
# family additionally needs ete3 + FastME/RapidNJ/Ninja (see the note at the end);
# use `--native` for a fully self-contained (topology-only) result.

# ---- build: fully static musl binary ----
FROM rust:1-alpine AS build
RUN apk add --no-cache musl-dev
WORKDIR /src
COPY Cargo.toml Cargo.lock ./
COPY src ./src
COPY tests ./tests
RUN cargo build --release && strip target/release/grapetree

# ---- runtime: binary + NumPy shims (default = byte-identical) ----
FROM python:3.11-slim AS runtime
LABEL org.opencontainers.image.title="grapetree-rs" \
      org.opencontainers.image.description="Fast Rust port of GrapeTree's tree/distance engine — byte-identical MSTreeV2, no edmonds binary" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.source="https://github.com/genpat-it/grapetree-rs"

# NumPy pinned to the reference's version: the shims reproduce NumPy's own float32
# SIMD summation and np.log, so matching NumPy keeps the output byte-identical
# (subject to the usual CPU-SIMD reproducibility caveat — see DECISIONS.md).
RUN pip install --no-cache-dir "numpy==1.26.4"

# Binary and shims live side by side so the binary's resolve_bundled() finds shim/
# next to the executable at runtime.
COPY --from=build /src/target/release/grapetree /opt/grapetree-rs/grapetree
COPY shim /opt/grapetree-rs/shim

ENV PATH="/opt/grapetree-rs:${PATH}" \
    GT_PYTHON=python3

ENTRYPOINT ["grapetree"]
CMD ["--help"]

# ---- NJ family (optional) ----
# `NJ`/`RapidNJ`/`ninja` need ete3 + the FastME/RapidNJ/Ninja binaries (and a JRE
# for ninja) for byte-identity, exactly like upstream GrapeTree. To enable them,
# extend this stage with:
#   RUN pip install --no-cache-dir "ete3<3.1.4" \
#    && apt-get update && apt-get install -y --no-install-recommends default-jre-headless \
#    && rm -rf /var/lib/apt/lists/*
#   COPY binaries /opt/grapetree-rs/binaries
# (ete3 requires Python < 3.13.) Otherwise `--native` gives a self-contained,
# topology-identical NJ with no extra dependencies.
