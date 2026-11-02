// moorhen-squab: SAXS profiles from macromolecular structures (Toby E King).
// C++ core (Emscripten) + optional WebGPU Debye path + profile fitting helpers.

export {
    SquabCalculator,
    defaultQGrid,
    DEFAULT_SOLVENT_DENSITY,
    DEFAULT_SHELL_CONTRAST,
} from "./SquabCalculator.js";
export type {
    SquabCreateOptions,
    SquabTables,
    PrepareOptions,
    PrepareResult,
    HydrationParameters,
    FitOptions,
    FitToProfileResult,
} from "./SquabCalculator.js";

export { SquabWorkerCalculator } from "./SquabWorkerCalculator.js";
export type { SquabWorkerOptions } from "./SquabWorkerCalculator.js";

export { SquabDebyeGPU, calculateDebyeCurveGPU, isWebGPUAvailable } from "./gpu.js";
export type { SquabGPUOptions } from "./gpu.js";

export { fitSAXSScale, fitSAXSScaleWithOffset } from "./vendor/utils/SAXSScaling.js";
export type { SAXSCurve, SAXSFitResult } from "./vendor/utils/SAXSScaling.js";
export { parseSAXSProfile } from "./vendor/utils/SAXSProfileParser.js";
export { calculateDebyeReference } from "./vendor/utils/DebyeReference.js";
export type { DebyeReferenceResult } from "./vendor/utils/DebyeReference.js";
export { compareSAXSProfiles } from "./vendor/types/SAXSComparison.js";
export type { SAXSComparison } from "./vendor/types/SAXSComparison.js";
export type { SAXSProfile } from "./vendor/types/SAXSProfile.js";

export type { SaxsModule, FloatVector } from "./wasm/saxs.js";

import { parseSAXSProfile } from "./vendor/utils/SAXSProfileParser.js";
import type { SAXSProfile } from "./vendor/types/SAXSProfile.js";

/** Parse a q / I(q) / σ text profile (.dat etc.) that is already in memory. */
export function parseSAXSProfileText(text: string, name = "profile.dat"): Promise<SAXSProfile> {
    return parseSAXSProfile(new File([text], name));
}

/** Low-level Emscripten factory, for direct use of the C++ bindings. */
export { default as createSaxsModule } from "./wasm/saxs.js";
