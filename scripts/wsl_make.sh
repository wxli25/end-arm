#!/bin/bash
# helper: run make targets inside WSL from Windows:  wsl.exe -d Ubuntu-24.04 -u root -- bash scripts/wsl_make.sh <END_N> <BACKEND> <target...> [VAR=val ...]
cd "$(dirname "$0")/.." || exit 1
N=$1; B=$2; shift 2
exec make -s END_N="$N" BACKEND="$B" "$@"
