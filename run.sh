#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="${BUILD_DIR:-$ROOT_DIR/build}"

if (($# == 0)) || [[ ! "$1" =~ ^[0-9]+$ ]]; then
  echo "Usage: $0 <run_number> [simulation|analysis] [config.json [simulation.json]]" >&2
  exit 2
fi

RUN_NUMBER="$1"
shift
MODE="${1:-simulation}"

bash "$ROOT_DIR/setup.sh"

if (($# > 0)); then
  shift
fi

resolve_config_path() {
  if [[ "$1" = /* ]]; then
    printf '%s\n' "$1"
  else
    printf '%s/%s\n' "$ROOT_DIR" "$1"
  fi
}

case "$MODE" in
  simulation)
    if (($# > 3)); then
      echo "Usage: $0 <run_number> simulation [config.json [simulation.json]]" >&2
      exit 2
    fi
    if (($# > 0)); then shift; fi

    cmake -S "$ROOT_DIR" -B "$BUILD_DIR"
    cmake --build "$BUILD_DIR" --target drift
    mkdir -p "$ROOT_DIR/results"

    case $# in
      0)
        (cd "$BUILD_DIR" && ./drift "$RUN_NUMBER")
        ;;
      1)
        config_path="$(resolve_config_path "$1")"
        (cd "$BUILD_DIR" && ./drift "$RUN_NUMBER" "$config_path")
        ;;
      2)
        config_path="$(resolve_config_path "$1")"
        simulation_config_path="$(resolve_config_path "$2")"
        (cd "$BUILD_DIR" && ./drift "$RUN_NUMBER" "$config_path" "$simulation_config_path")
        ;;
    esac
    ;;
  analysis)
    if (($# > 1)); then
      echo "Usage: $0 <run_number> analysis" >&2
      exit 2
    fi
    if (($# > 0)); then shift; fi
    command -v root >/dev/null 2>&1 || {
      echo "ROOT executable not found in PATH." >&2
      exit 1
    }
    (cd "$ROOT_DIR" && root -l -q "analyse_results.C(\"$RUN_NUMBER\")")
    ;;
  *)
    echo "Usage: $0 <run_number> [simulation|analysis] [config.json [simulation.json]]" >&2
    exit 2
    ;;
esac