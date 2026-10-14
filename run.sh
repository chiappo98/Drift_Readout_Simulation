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

backup_previous_run() {
  local results_dir="$ROOT_DIR/results"
  local backup_root="$results_dir/backup"
  local backup_path
  local timestamp
  local suffix=1
  local -a matches=()
  local -a run_files=()

  shopt -s nullglob
  matches=("$results_dir/run${RUN_NUMBER}."* "$results_dir/run${RUN_NUMBER}_"*)
  shopt -u nullglob

  for backup_path in "${matches[@]}"; do
    if [[ -f "$backup_path" ]]; then
      run_files+=("$backup_path")
    fi
  done
  ((${#run_files[@]} > 0)) || return 0

  timestamp="$(date '+%Y%m%d_%H%M%S')"
  backup_path="$backup_root/run${RUN_NUMBER}_${timestamp}"
  while [[ -e "$backup_path" ]]; do
    backup_path="$backup_root/run${RUN_NUMBER}_${timestamp}_${suffix}"
    ((suffix += 1))
  done

  mkdir -p "$backup_path"
  mv -- "${run_files[@]}" "$backup_path/"
  printf 'Backup del run %s creato in %s\n' "$RUN_NUMBER" "$backup_path"
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
    backup_previous_run

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