#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p \
  "$ROOT_DIR/build" \
  "$ROOT_DIR/config" \
  "$ROOT_DIR/results" \
  "$ROOT_DIR/gas_files" \
  "$ROOT_DIR/maps/mylar" \
  "$ROOT_DIR/maps/mesh" \
  "$ROOT_DIR/plot/mylar" \
  "$ROOT_DIR/plot/mesh"

printf 'Directory di progetto pronte in %s\n' "$ROOT_DIR"