// Copies the Emscripten output (npm/wasm) and its hand-written typings into dist/wasm.
import { copyFileSync, existsSync, mkdirSync } from "node:fs";
import { dirname, join, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const npmDir = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const out = join(npmDir, "dist", "wasm");
mkdirSync(out, { recursive: true });

const js = join(npmDir, "wasm", "saxs.js");
if (!existsSync(js)) {
    console.error("npm/wasm/saxs.js not found - run `npm run build:wasm` first (needs emsdk).");
    process.exit(1);
}

copyFileSync(js, join(out, "saxs.js"));
copyFileSync(join(npmDir, "src", "wasm", "saxs.d.ts"), join(out, "saxs.d.ts"));
for (const optional of ["saxs.wasm", "BUILD_INFO.txt"]) {
    const p = join(npmDir, "wasm", optional);
    if (existsSync(p)) copyFileSync(p, join(out, optional));
}
console.log("copy-wasm: dist/wasm ready");
