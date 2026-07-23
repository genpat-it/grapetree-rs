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

# Runtime deps for all methods:
#  - NumPy (pinned to the reference's version): the MSTree/MSTreeV2 shims reproduce
#    NumPy's own float32 SIMD summation and np.log — matching NumPy keeps the output
#    byte-identical (subject to the usual CPU-SIMD caveat — see DECISIONS.md).
#  - ete3 (+ its undeclared `six` dep): the NJ-family post-processing shim
#    (midpoint root / unroot / write). Needs Python < 3.13 (this image is 3.11).
#  - a headless JRE: to run the bundled Ninja.jar for the `ninja` method.
RUN apt-get update \
 && apt-get install -y --no-install-recommends default-jre-headless \
 && rm -rf /var/lib/apt/lists/* \
 && pip install --no-cache-dir "numpy==1.26.4" "ete3==3.1.3" six

# Binary + shims + the bundled FastME/RapidNJ/Ninja binaries, side by side so the
# binary's resolve_bundled() finds shim/ and binaries/ next to the executable.
COPY --from=build /src/target/release/grapetree /opt/grapetree-rs/grapetree
COPY shim /opt/grapetree-rs/shim
COPY binaries /opt/grapetree-rs/binaries

ENV PATH="/opt/grapetree-rs:${PATH}" \
    GT_PYTHON=python3

ENTRYPOINT ["grapetree"]
CMD ["--help"]

# Coverage in this image:
#  - byte-identical: distance, MSTree, MSTreeV2 (NumPy shims + pure-Rust edmonds),
#    NJ, RapidNJ (FastME/RapidNJ binaries + ete3 shim).
#  - `ninja` RUNS and produces a valid NJ tree (Ninja.jar via a clean `java -jar`
#    invocation — grapetree-rs fixes upstream's broken caller: `-d64` removed in
#    Java ≥ 9, an `-Xmx` sized at 90% of total RAM, and an ete3 parse). It cannot be
#    *byte-validated* against upstream because upstream's `ninja` produces no output
#    at all on a modern JVM; we cross-check it as RF=0 to the byte-identical NJ.
