#!/usr/bin/env bash
set -euo pipefail
packet_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
"${PYTHON:-python3.12}" "$packet_dir/controls.py" | node "$packet_dir/verify.cjs" "$@"
