import createSaxsModule from "./wasm/saxs.js";
import type { FloatVector, SaxsModule, SaxsModuleOptions } from "./wasm/saxs.js";
import {
    elementsInfoJson,
    residueAtomTypesJson,
    saxsAtomTypesJson,
} from "./generated/assets.js";
import {
    fitSAXSScale,
    fitSAXSScaleWithOffset,
} from "./vendor/utils/SAXSScaling.js";
import type { SAXSCurve, SAXSFitResult } from "./vendor/utils/SAXSScaling.js";
import type { SAXSProfile } from "./vendor/types/SAXSProfile.js";

/** Defaults used by Toby's frontend (App.tsx). */
export const DEFAULT_SOLVENT_DENSITY = 0.334;
export const DEFAULT_SHELL_CONTRAST = 0.0167;

export interface SquabTables {
    elementsInfoJson: string;
    saxsAtomTypesJson: string;
    residueAtomTypesJson: string;
}

export interface SquabCreateOptions {
    /** Echo C++ stdout/stderr to the console. Default false (output is still kept in `log`). */
    verbose?: boolean;
    /** Lines of C++ output to keep for error messages / `log`. Default 500. */
    logLines?: number;
    /** Only needed for a build with SQUAB_SINGLE_FILE=OFF: where to find saxs.wasm. */
    locateFile?: (path: string, prefix: string) => string;
    /** Override the form-factor / atom-type tables bundled from frontend/public. */
    tables?: Partial<SquabTables>;
}

export interface PrepareOptions {
    /**
     * Used only to pick the parser: names ending .pdb/.ent are read as PDB,
     * anything else as mmCIF. Default "model.cif".
     */
    filename?: string;
    /** q grid in Å⁻¹. Default: `defaultQGrid()` (500 points, 0.01–0.5). */
    qValues?: ArrayLike<number>;
}

export interface PrepareResult {
    atomCount: number;
    multipoleL: number;
    qValues: number[];
}

export interface HydrationParameters {
    solventDensity?: number;
    shellContrast?: number;
    /** PEPSI-style Δr/r, must be within [-0.05, 0.05]. */
    radiusAdjustmentFraction?: number;
}

export interface FitOptions {
    /** Fit a constant background as well as the scale (fitSAXSScaleWithOffset). */
    fitBackgroundOffset?: boolean;
    solventDensity?: number;
    /** Grid points per parameter. Toby's app uses 100 (10,000 evaluations). */
    gridSize?: number;
    radiusAdjustmentRange?: [number, number];
    shellContrastRange?: [number, number];
    /** Called after each row of the grid with a fraction in [0, 1]. */
    onProgress?: (fraction: number) => void;
    /** Abort the search; the promise rejects with an AbortError. */
    signal?: AbortSignal;
}

export interface FitToProfileResult {
    curve: SAXSCurve;
    fit: SAXSFitResult;
    radiusAdjustmentFraction: number;
    shellContrast: number;
    solventDensity: number;
    evaluated: number;
}

/** Evenly spaced q grid, as used by the frontend when no experimental profile is loaded. */
export function defaultQGrid(count = 500, qMin = 0.01, qMax = 0.5): number[] {
    if (count < 2) return [qMin];
    return Array.from({ length: count }, (_, i) => qMin + (i * (qMax - qMin)) / (count - 1));
}

function linspace(lo: number, hi: number, n: number): number[] {
    if (n < 2) return [lo];
    return Array.from({ length: n }, (_, i) => lo + (i * (hi - lo)) / (n - 1));
}

const yieldToEventLoop = () => new Promise<void>((resolve) => setTimeout(resolve, 0));

/**
 * One Emscripten instance of Toby's SAXS core.
 *
 * The C++ side keeps a single prepared structure in globals, so a calculator
 * holds one structure at a time; create several calculators for several
 * structures (each has its own WASM memory).
 */
export class SquabCalculator {
    readonly module: SaxsModule;
    private readonly tables: SquabTables;
    private readonly logBuffer: string[] = [];
    private readonly logLimit: number;
    private prepared: PrepareResult | null = null;

    private constructor(module: SaxsModule, tables: SquabTables, logBuffer: string[], logLimit: number) {
        this.module = module;
        this.tables = tables;
        this.logBuffer = logBuffer;
        this.logLimit = logLimit;
    }

    static async create(options: SquabCreateOptions = {}): Promise<SquabCalculator> {
        const logLimit = options.logLines ?? 500;
        const logBuffer: string[] = [];
        const record = (echo: (s: string) => void) => (text: string) => {
            logBuffer.push(text);
            if (logBuffer.length > logLimit) logBuffer.splice(0, logBuffer.length - logLimit);
            if (options.verbose) echo(text);
        };

        const moduleOptions: SaxsModuleOptions = {
            print: record((s) => console.log("[squab]", s)),
            printErr: record((s) => console.warn("[squab]", s)),
        };
        if (options.locateFile) moduleOptions.locateFile = options.locateFile;

        const module = await createSaxsModule(moduleOptions);
        const tables: SquabTables = {
            elementsInfoJson,
            saxsAtomTypesJson,
            residueAtomTypesJson,
            ...options.tables,
        };
        return new SquabCalculator(module, tables, logBuffer, logLimit);
    }

    /** Recent C++ stdout/stderr lines (diagnostics printed by prepare). */
    get log(): readonly string[] {
        return this.logBuffer;
    }

    get isPrepared(): boolean {
        return this.prepared !== null;
    }

    get preparedState(): PrepareResult | null {
        return this.prepared;
    }

    /**
     * Load a structure (mmCIF by default) and precompute everything on the q grid.
     * Must be called again whenever the structure or the q grid changes.
     * To fit an experimental profile, pass `qValues: profile.q`.
     */
    prepare(structureText: string, options: PrepareOptions = {}): PrepareResult {
        if (!structureText || !structureText.trim()) {
            throw new Error("moorhen-squab: empty structure text");
        }
        const qValues = Array.from(options.qValues ?? defaultQGrid());
        if (qValues.length === 0) throw new Error("moorhen-squab: empty q grid");

        const filename = options.filename ?? "model.cif";
        const qVector = new this.module.FloatVector();
        this.prepared = null;
        this.logBuffer.length = 0;
        try {
            for (const q of qValues) qVector.push_back(q);
            const m = this.module;
            const t = this.tables;
            const args = [t.elementsInfoJson, t.saxsAtomTypesJson, t.residueAtomTypesJson, structureText, filename, qVector] as const;
            this.safeCall(
                m.squab_prepare ? () => m.squab_prepare!(...args) : null,
                () => m.prepare_debye_gpu(...args)
            );
        } finally {
            qVector.delete();
        }

        const atomCount = this.withVector(this.module.get_positions(), (v) => v.size() / 3);
        if (atomCount === 0) {
            throw this.error("no atoms were extracted from the structure");
        }
        this.prepared = { atomCount, multipoleL: this.module.get_multipole_L(), qValues };
        return this.prepared;
    }

    /** Production curve: atomic + excluded volume + hydration shell (Toby's `hydrationCurve`). */
    calculateCurve(params: HydrationParameters = {}): SAXSCurve {
        const state = this.requirePrepared();
        const m = this.module;
        const args = [
            params.solventDensity ?? DEFAULT_SOLVENT_DENSITY,
            params.shellContrast ?? DEFAULT_SHELL_CONTRAST,
            params.radiusAdjustmentFraction ?? 0,
        ] as const;
        const intensities = this.readVector(
            () => m.get_hydration_combined_multipole_intensities(...args),
            m.squab_hydration_combined_multipole_intensities
                ? () => m.squab_hydration_combined_multipole_intensities!(...args)
                : null
        );
        return { qValues: [...state.qValues], intensities };
    }

    /** In-vacuo multipole curve at order L (default: the L chosen by prepare). */
    atomicCurve(L?: number): SAXSCurve {
        const state = this.requirePrepared();
        const m = this.module;
        const order = L ?? state.multipoleL;
        const intensities = this.readVector(
            () => m.get_multipole_intensities(order),
            m.squab_multipole_intensities ? () => m.squab_multipole_intensities!(order) : null
        );
        return { qValues: [...state.qValues], intensities };
    }

    /** Atomic + excluded volume, no hydration shell. */
    excludedVolumeCurve(solventDensity = DEFAULT_SOLVENT_DENSITY): SAXSCurve {
        const state = this.requirePrepared();
        const m = this.module;
        const intensities = this.readVector(
            () => m.get_combined_multipole_intensities(solventDensity),
            m.squab_combined_multipole_intensities
                ? () => m.squab_combined_multipole_intensities!(solventDensity)
                : null
        );
        return { qValues: [...state.qValues], intensities };
    }

    /**
     * Toby's two-parameter grid search (radius adjustment × shell contrast),
     * scoring each candidate by reduced χ² against the experimental profile.
     * The calculator must have been prepared on the profile's own q grid.
     */
    async fitToProfile(profile: SAXSProfile, options: FitOptions = {}): Promise<FitToProfileResult> {
        const state = this.requirePrepared();
        if (profile.q.length !== state.qValues.length ||
            profile.q.some((q, i) => Math.abs(q - state.qValues[i]) > 1e-6 * Math.max(1, Math.abs(q)))) {
            throw new Error(
                "moorhen-squab: prepare() was not run on this profile's q grid; " +
                "call prepare(structure, { qValues: profile.q }) first"
            );
        }

        const n = options.gridSize ?? 100;
        const solventDensity = options.solventDensity ?? DEFAULT_SOLVENT_DENSITY;
        const radii = linspace(...(options.radiusAdjustmentRange ?? [-0.05, 0.05]), n);
        const contrasts = linspace(...(options.shellContrastRange ?? [0, 0.0334]), n);
        const fitCurve = options.fitBackgroundOffset ? fitSAXSScaleWithOffset : fitSAXSScale;

        let best: FitToProfileResult | null = null;
        let evaluated = 0;

        for (let row = 0; row < radii.length; row++) {
            if (options.signal?.aborted) {
                throw new DOMException("SAXS fit aborted", "AbortError");
            }
            const radiusAdjustmentFraction = radii[row];
            for (const shellContrast of contrasts) {
                const curve = this.calculateCurve({ solventDensity, shellContrast, radiusAdjustmentFraction });
                const fit = fitCurve(profile, curve);
                evaluated++;
                if (Number.isFinite(fit.reducedChiSquared) &&
                    (!best || fit.reducedChiSquared < best.fit.reducedChiSquared)) {
                    best = { curve, fit, radiusAdjustmentFraction, shellContrast, solventDensity, evaluated: 0 };
                }
            }
            options.onProgress?.((row + 1) / radii.length);
            await yieldToEventLoop();
        }

        if (!best) throw new Error("moorhen-squab: grid search produced no finite fit");
        best.evaluated = evaluated;
        return best;
    }

    /** Flattened xyz positions for the GPU Debye path. */
    gpuPositions(): Float32Array {
        this.requirePrepared();
        return Float32Array.from(this.readVector(() => this.module.get_positions()));
    }

    /** Amplitudes [q][atom] for the GPU Debye path. */
    gpuAmplitudes(): Float32Array {
        this.requirePrepared();
        return Float32Array.from(this.readVector(() => this.module.get_amplitudes()));
    }

    // ------------------------------------------------------------ internals

    private requirePrepared(): PrepareResult {
        if (!this.prepared) throw new Error("moorhen-squab: call prepare() first");
        return this.prepared;
    }

    private readVector(make: () => FloatVector, safeMake: (() => FloatVector) | null = null): number[] {
        const vector = this.safeCall(safeMake, make, (v) => v.delete());
        return this.withVector(vector, (v) => {
            const out = new Array<number>(v.size());
            for (let i = 0; i < out.length; i++) out[i] = v.get(i);
            return out;
        });
    }

    private withVector<T>(vector: FloatVector, use: (v: FloatVector) => T): T {
        try {
            return use(vector);
        } finally {
            vector.delete();
        }
    }

    /**
     * Prefer the squab_* binding (C++ exceptions caught, text in squab_last_error);
     * fall back to Toby's original binding for modules built without them.
     */
    private safeCall<T>(safe: (() => T) | null, unsafe: () => T, discard?: (result: T) => void): T {
        const lastError = this.module.squab_last_error;
        if (!safe || !lastError) return this.callWasm(unsafe);
        const result = this.callWasm(safe);
        const message = lastError.call(this.module);
        if (message) {
            try { discard?.(result); } catch { /* ignore */ }
            throw this.error(message);
        }
        return result;
    }

    private callWasm<T>(fn: () => T): T {
        try {
            return fn();
        } catch (e) {
            throw this.error(this.describeException(e), e);
        }
    }

    private describeException(e: unknown): string {
        if (e instanceof Error) return e.message;
        if (e && typeof e === "object" && "excPtr" in e) {
            // Module built before cpp/squab_npm_bindings.cpp existed: the text is not reachable.
            return "C++ exception (rebuild moorhen-squab to see its message; see the C++ output below)";
        }
        return String(e);
    }

    private error(message: string, cause?: unknown): Error {
        const tail = this.logBuffer.slice(-15).join("\n");
        const err = new Error(`moorhen-squab: ${message}` + (tail ? `\n--- C++ output (last lines) ---\n${tail}` : ""));
        if (cause !== undefined) (err as Error & { cause?: unknown }).cause = cause;
        return err;
    }
}
