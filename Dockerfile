# Multi-stage build → tiny runtime image with a static musl binary.
FROM rust:1-alpine AS build
RUN apk add --no-cache musl-dev
WORKDIR /src
COPY Cargo.toml ./
COPY src ./src
COPY tests ./tests
RUN cargo build --release && strip target/release/grapetree

FROM alpine:3.20
LABEL org.opencontainers.image.title="grapetree-rs" \
      org.opencontainers.image.description="Fast Rust port of GrapeTree's tree/distance engine (experimental)" \
      org.opencontainers.image.licenses="GPL-3.0-or-later" \
      org.opencontainers.image.source="https://github.com/genpat-it/grapetree-rs"
COPY --from=build /src/target/release/grapetree /usr/local/bin/grapetree
ENTRYPOINT ["grapetree"]
CMD ["--help"]
