# moorhen-squab (NPM packaging of MoorhenSquab)

This folder turns Toby King's MoorhenSquab SAXS calculator into an ES-module
NPM package that a Moorhen plugin (e.g. MoorhenWebComponentServer) can import.

Nothing outside `npm/` is modified. The C++ in `../src` + `../include`, the
TypeScript helpers in `../src/utils`, `../src/types`, `../frontend/src/webgpu`,
and the tables/shaders in `../frontend/public` are all taken as they are at
build time, so the package always reflects Toby's current code.

## What's in the package

| Export | What it is |
| --- | --- |
| `SquabCalculator` | One WASM instance of the C++ core. `prepare(mmcif)` → `calculateCurve()` / `fitToProfile(profile)` |
| `SquabDebyeGPU`, `calculateDebyeCurveGPU` | Toby's WebGPU Debye path with the WGSL shaders bundled |
| `parseSAXSProfileText`, `parseSAXSProfile` | q / I(q) / σ profile parsing (ATSAS q-unit heuristic) |
| `fitSAXSScale`, `fitSAXSScaleWithOffset`, `compareSAXSProfiles`, `calculateDebyeReference` | Toby's fitting / reference helpers |
| `createSaxsModule` | The raw Emscripten factory, if you want the bindings directly |

The WASM is embedded in `dist/wasm/saxs.js` (`-sSINGLE_FILE`), and the JSON
tables and shaders are compiled in as strings, so consumers need no
`public/` asset copying and no special bundler config.

## Building locally (no NPM account needed)

Prerequisites: an activated emsdk, Node ≥ 18, CMake ≥ 3.18, git (CMake
fetches gemmi v0.7.0 and nlohmann/json v3.11.3 from GitHub unless
`../checkout/gemmi` is populated).

```bash
cd ~/MoorhenSquab/npm
source ~/emsdk/emsdk_env.sh
npm install          # typescript + @webgpu/types (dev only)
npm run build        # WASM (build-wasm/ → wasm/) then TypeScript (→ dist/)
npm test             # 1orc.cif vs. native reference + fit round-trip
npm pack             # → moorhen-squab-0.3.1.tgz
```

After a C++ change: `npm run build`. After a TS-only change: `npm run build:ts`.
Debug WASM build: `npm run build:wasm -- -DSQUAB_DEBUG=ON`.
Separate `.wasm` file instead of embedded: `npm run build:wasm -- -DSQUAB_SINGLE_FILE=OFF`
(then pass `locateFile` to `SquabCalculator.create`).

## Using it from another project

```bash
cd ~/MoorhenWebComponentServer
npm install ../MoorhenSquab/npm/moorhen-squab-0.3.1.tgz
```

Install the packed tarball rather than `file:../MoorhenSquab/npm`: a folder
install is a symlink, which pulls in the source tree and confuses Vite's
dependency pre-bundling. Re-run `npm pack` + `npm install <tgz>` after rebuilding.

```ts
import { SquabCalculator, parseSAXSProfileText } from "moorhen-squab";

const calc = await SquabCalculator.create();          // one per structure you want to keep prepared

// mmCIF from Moorhen, either via coot…
const mmcif = await molecule.getAtoms("mmcif");
// …or via gemmi:
// const mmcif = window.gemmiModule.get_mmcif_string_from_gemmi_struct(molecule.gemmiStructure);

// Theoretical curve on the default grid (500 points, q = 0.01–0.5 Å⁻¹)
calc.prepare(mmcif);
const curve = calc.calculateCurve();                  // { qValues, intensities }

// Fit to experimental data: prepare on the profile's own q grid first
const profile = await parseSAXSProfileText(datFileText, "data.dat");
calc.prepare(mmcif, { qValues: profile.q });
const best = await calc.fitToProfile(profile, { onProgress: f => console.log(f) });
// best.curve, best.fit.scale, best.fit.reducedChiSquared,
// best.radiusAdjustmentFraction, best.shellContrast
```

`prepare` accepts PDB text too: pass `{ filename: "model.pdb" }` (the
extension picks the parser; anything not `.pdb`/`.ent` is read as mmCIF).

C++ stdout (Toby's diagnostics) is captured rather than printed; read
`calc.log`, or pass `{ verbose: true }` to `create`. Errors thrown in C++
come back as JS `Error`s with the exception text and the last log lines.

## Running in a Web Worker (`SquabWorkerCalculator`)

The C++ runs synchronously, so on the main thread a large structure or a fit
freezes the page. `SquabWorkerCalculator` has the same calls as
`SquabCalculator`, returning promises, with the work done in a Web Worker:

```ts
import { SquabWorkerCalculator } from "moorhen-squab";

const calc = new SquabWorkerCalculator();          // worker starts on first call
calc.onBusyChange(busy => …);                      // true while anything is in flight
await calc.prepare(mmcif, { qValues, signal });    // aborting restarts the worker
const curve = await calc.calculateCurve();
const best = await calc.fitToProfile(profile, { onProgress, signal }); // abort stops between grid rows
calc.cancel();      // kill whatever is running (next call gets a fresh worker)
calc.terminate();   // done with it
```

The worker is shipped inside the package: `npm run build:ts` runs
`scripts/build-worker.mjs`, which inlines the worker entry (`src/worker/`),
`SquabCalculator`, the tables and the Emscripten module into one module source
(`dist/worker/workerSource.js`, checked to round-trip exactly). The client starts
it from a Blob URL, so hosts need no bundler configuration. Under a
Content-Security-Policy without `worker-src blob:`, pass `createWorker` to supply
your own `Worker`. `npm test` runs `test/worker-smoke.mjs`, which drives the real
bundle in a thread and checks results, progress, abort/cancel, errors, and that
the main thread keeps running during a fit.

## Moorhen plugin UI (`moorhen-squab/react`)

A ready-made SAXS plugin for Moorhen ≥ 1.0.2 (React ≥ 18 host). Put it inside the
web component:

```tsx
import { SaxsPlugin } from "moorhen-squab/react";

<moorhen-web-component id="my-moorhen">
    <SaxsPlugin />
</moorhen-web-component>
```

It finds the enclosing `<moorhen-web-component>`, waits for its instance, then:

* renders its UI as `slot="saxs"` light-DOM content (with its own styles, using
  Moorhen's `--moorhen-*` theme variables);
* registers a Moorhen side panel showing that slot (or a modal; the toolbar's
  *Open as window* / *Dock in side panel* button switches and reopens it);
* adds a **SAXS** main-menu entry and **Calculate → SAXS curves…**;
* runs every calculation in a `SquabWorkerCalculator`, so Moorhen stays usable,
  with progress and a **Cancel** button in the panel and a Moorhen snackbar
  ("SAXS: fitting … 40%", then done/failed) that shows even when the panel is closed.

The UI lists the loaded molecules (`instance.getMoleculeList()`, kept current via
the store subscription) with a checkbox each; new molecules start ticked. It covers
calculating curves for the ticked molecules (mmCIF from
`molecule.getAtoms("mmcif")`, i.e. `molecules_container_t`), loading a q / I(q) / σ
profile, the grid-search fit, a log I(q) plot, a results table and CSV export.
Because it is rendered by the host page's React and only projected into Moorhen's
frame, results survive closing the panel and switching placement.

Props (all optional): `id` (panel/modal id and slot name, default `"saxs"`),
`defaultPlacement` (`"sidePanel"` | `"modal"`), `addMenuEntries` (default `true`),
`instance` (pass a Moorhen instance instead of looking for the web component).
`useSaxsPlugin`, `SaxsPanel` and `SaxsPlot` are exported for hosts that want their
own layout. The package does not depend on `moorhen`: it uses the small structural
types in `src/react/moorhenTypes.ts`.

## Things worth knowing

* `prepare()` runs Toby's full pipeline in `main.cpp`, including his
  validation diagnostics (interpolation checks, several shell contrasts).
  That's extra work on every prepare; a lean entry point is an easy later win.
* `SquabCalculator` itself is synchronous and runs on the calling thread; use
  `SquabWorkerCalculator` in a page.
* The C++ keeps state in globals, so each `SquabCalculator` holds one
  prepared structure; separate calculators have separate WASM memory.
* The bundled gemmi is independent of Moorhen's own gemmi build.

## Licences

MoorhenSquab is Apache-2.0. The WASM includes gemmi (MPL-2.0,
https://github.com/project-gemmi/gemmi) and nlohmann/json (MIT).
`test/1orc.cif` is PDB entry 1ORC converted to mmCIF with gemmi.
