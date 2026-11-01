#!/usr/bin/env bash
set -euo pipefail
packet="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
node "$packet/verify.cjs" "$@" --pins-only
python3 "$packet/controls.py" | node "$packet/verify.cjs" "$@"
