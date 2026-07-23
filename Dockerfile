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

# Runtime deps for byte-identity across ALL methods:
#  - NumPy (pinned to the reference's version): the MSTree/MSTreeV2 shims reproduce
#    NumPy's own float32 SIMD summation and np.log — matching NumPy keeps the output
#    byte-identical (subject to the usual CPU-SIMD caveat — see DECISIONS.md).
#  - ete3: the NJ-family post-processing shim (midpoint root / unroot / write),
#    exactly the reference toolchain. Needs Python < 3.13 (this image is 3.11).
RUN pip install --no-cache-dir "numpy==1.26.4" "ete3==3.1.3" six

# Binary + shims + the bundled FastME/RapidNJ/Ninja binaries, side by side so the
# binary's resolve_bundled() finds shim/ and binaries/ next to the executable.
COPY --from=build /src/target/release/grapetree /opt/grapetree-rs/grapetree
COPY shim /opt/grapetree-rs/shim
COPY binaries /opt/grapetree-rs/binaries

ENV PATH="/opt/grapetree-rs:${PATH}" \
    GT_PYTHON=python3

ENTRYPOINT ["grapetree"]
CMD ["--help"]

# Byte-identical coverage in this image: distance, MSTree, MSTreeV2 (NumPy shims,
# pure-Rust edmonds), and NJ / RapidNJ (FastME/RapidNJ binaries + ete3 shim).
# `ninja` is NOT usable — upstream's `ninja` calls `java -d64`, removed in Java ≥ 9,
# so even the reference can't run it; no JRE is bundled. `--native` gives a fully
# self-contained (topology-only, not byte-identical) run for the NJ family.
