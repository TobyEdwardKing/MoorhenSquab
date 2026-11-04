#!/usr/bin/env bash
# Build Toby's C++ core to npm/wasm/saxs.js (+ saxs.wasm if SINGLE_FILE is off).
# Needs an activated emsdk:   source ~/emsdk/emsdk_env.sh
# Extra arguments go to CMake, e.g.   npm run build:wasm -- -DSQUAB_DEBUG=ON
set -euo pipefail

here="$(cd "$(dirname "$0")/.." && pwd)"
build="$here/build-wasm"

if ! command -v emcmake >/dev/null 2>&1; then
    echo "emcmake not found on PATH. Activate emsdk first, e.g.:" >&2
    echo "    source ~/emsdk/emsdk_env.sh" >&2
    exit 1
fi

emcmake cmake -S "$here" -B "$build" -DCMAKE_BUILD_TYPE=Release "$@"
cmake --build "$build" -j

mkdir -p "$here/wasm"
cp "$build/saxs.js" "$here/wasm/saxs.js"
if [ -f "$build/saxs.wasm" ]; then
    cp "$build/saxs.wasm" "$here/wasm/saxs.wasm"
else
    rm -f "$here/wasm/saxs.wasm"
fi

{
    echo "emscripten: $(emcc --version | head -n 1)"
    echo "MoorhenSquab commit: $(git -C "$here/.." rev-parse --short HEAD 2>/dev/null || echo unknown)"
    echo "built: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$here/wasm/BUILD_INFO.txt"

echo "WASM module written to $here/wasm/"
ls -l "$here/wasm/"
