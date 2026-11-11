// Smoke test for the built package: loads dist/, runs 1orc.cif through the
// C++ core and compares with intensities from a native build of the same
// sources. Run with `npm test` after `npm run build`.
import { readFileSync } from "node:fs";
import { dirname, join } from "node:path";
import { fileURLToPath } from "node:url";

const here = dirname(fileURLToPath(import.meta.url));
const { SquabCalculator, defaultQGrid, parseSAXSProfileText } = await import("../dist/index.js");

const reference = JSON.parse(readFileSync(join(here, "reference-1orc.json"), "utf8"));
const cif = readFileSync(join(here, "1orc.cif"), "utf8");

async function createCalculator() {
    try {
        return await SquabCalculator.create();
    } catch (e) {
        // The shipped build targets web + worker only. If this Node version
        // trips Emscripten's environment detection, pretend to be a worker.
        console.warn("create() failed under Node, retrying with worker shims:", e?.message ?? e);
        globalThis.self ??= globalThis;
        globalThis.location ??= { href: import.meta.url };
        return await SquabCalculator.create();
    }
}

let failures = 0;
const check = (ok, message) => {
    console.log(`${ok ? "ok  " : "FAIL"}  ${message}`);
    if (!ok) failures++;
};

const t0 = performance.now();
const calc = await createCalculator();
const t1 = performance.now();
console.log(`module created in ${(t1 - t0).toFixed(0)} ms`);

const { count, qMin, qMax } = reference.qGrid;
const prepared = calc.prepare(cif, { filename: "1orc.cif", qValues: defaultQGrid(count, qMin, qMax) });
const t2 = performance.now();
console.log(`prepare() in ${(t2 - t1).toFixed(0)} ms`);

check(prepared.atomCount === reference.atomCount, `atom count ${prepared.atomCount} (expected ${reference.atomCount})`);
check(prepared.multipoleL === reference.multipoleL, `multipole L ${prepared.multipoleL} (expected ${reference.multipoleL})`);

const curve = calc.calculateCurve({ solventDensity: 0.334, shellContrast: 0.0167, radiusAdjustmentFraction: 0 });
let worst = 0;
curve.intensities.forEach((v, i) => {
    const ref = reference.intensities[i];
    worst = Math.max(worst, Math.abs(v - ref) / Math.abs(ref));
});
check(curve.intensities.length === count, `curve has ${curve.intensities.length} points`);
check(worst < 1e-3, `max relative deviation from native reference: ${worst.toExponential(2)}`);

// Round trip through the fitting path: a synthetic "experimental" profile made
// from the model at known parameters should be recovered by the grid search.
const target = calc.calculateCurve({ shellContrast: 0.02, radiusAdjustmentFraction: 0.01 });
const datText = target.qValues
    .map((q, i) => `${q} ${3.5 * target.intensities[i]} ${0.01 * target.intensities[i]}`)
    .join("\n");
const profile = await parseSAXSProfileText(datText, "synthetic.dat");
const fit = await calc.fitToProfile(profile, {
    gridSize: 11,
    radiusAdjustmentRange: [-0.05, 0.05],
    shellContrastRange: [0, 0.04],
});
check(Math.abs(fit.radiusAdjustmentFraction - 0.01) < 1e-9, `fit radius adjustment ${fit.radiusAdjustmentFraction.toFixed(3)} (expected 0.010)`);
check(Math.abs(fit.shellContrast - 0.02) < 1e-9, `fit shell contrast ${fit.shellContrast.toFixed(4)} (expected 0.0200)`);
check(Math.abs(fit.fit.scale - 3.5) < 1e-3, `fit scale ${fit.fit.scale.toFixed(4)} (expected 3.5)`);

let threw = false;
try {
    calc.prepare("this is not a structure");
} catch (e) {
    threw = true;
    console.log("      bad input error:", String(e.message).split("\n")[0]);
}
check(threw && !calc.isPrepared, "garbage input raises a readable error");

console.log(failures ? `\n${failures} check(s) failed` : "\nall checks passed");
process.exit(failures ? 1 : 0);
