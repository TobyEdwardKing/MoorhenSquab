
import { useEffect, useState } from "react";
import createSaxsModule from "./wasm/saxs.js";
import { DebyeGPU } from "./webgpu/DebyeGPU";
import { parseSAXSProfile } from "../../src/utils/SAXSProfileParser";
import type { SAXSProfile } from "../../src/types/SAXSProfile";
import {
    fitSAXSScale,
    type SAXSFitResult,
    type SAXSCurve
} from "../../src/utils/SAXSScaling";
import { SAXSPlot } from "./components/SAXSPlot";
import { calculateDebyeReference } from "../../src/utils/DebyeReference";

function App() {
    const [saxsProfile, setSaxsProfile] = useState<SAXSProfile | null>(null);
    const [calculatedCurve, setCalculatedCurve] = useState<SAXSCurve | null>(null);
    const [module, setModule] = useState<any>(null);
    const [filename, setFilename] = useState("");
    const [fileContents, setFileContents] = useState("");
    const [atomCount, setAtomCount] = useState<number | null>(null);
    const [summary, setSummary] = useState("");
    const [elementJSON, setElementJSON] = useState("");
    const [saxsAtomTypeJSON, setSaxsAtomTypeJSON] = useState("");
    const [residueAtomTypeJSON, setResidueAtomTypeJSON] = useState("");
    const [gpu, setGpu] = useState<DebyeGPU | null>(null);
    const [fitResult, setFitResult] = useState<SAXSFitResult | null>(null);
    const [isPreparing, setIsPreparing] = useState(false);
    const [isPrepared, setIsPrepared] = useState(false);

    const [structureText, setStructureText] = useState("");
    const [structureFilename, setStructureFilename] = useState("");
    const [bestParameters, setBestParameters] = useState<{
        radiusAdjustmentFraction: number;
        shellContrast: number;
    } | null>(null);
    // ---------------------------------------------------------
    // INITIALISE WASM + GPU
    // ---------------------------------------------------------

    useEffect(() => {
        async function loadModule() {
            try {
                const wasm = await createSaxsModule({
                    locateFile: (path: string) => {
                        if (path.endsWith(".wasm")) {
                            return "/saxs.wasm";
                        }
                        return path;
                    }
                });

                const elementResponse = await fetch("/elementsInfo.json");
                const saxsAtomTypeResponse = await fetch("/SAXSAtomTypes.json");
                const residueAtomTypeResponse = await fetch("/residueAtomTypes.json");

                const elementJson = await elementResponse.text();
                const saxsAtomTypeJson = await saxsAtomTypeResponse.text();
                const residueAtomTypeJson = await residueAtomTypeResponse.text();

                setElementJSON(elementJson);
                setSaxsAtomTypeJSON(saxsAtomTypeJson);
                setResidueAtomTypeJSON(residueAtomTypeJson);
                setModule(wasm);

                console.log("elementJSON:", elementJson.slice(0, 200));
                console.log("saxsAtomTypeJSON:", saxsAtomTypeJson.slice(0, 200));
                console.log("residueAtomTypeJSON:", residueAtomTypeJson.slice(0, 200));
                console.log(wasm);
                console.log(Object.keys(wasm));
                console.log("WASM loaded.");

                const gpuBackend = new DebyeGPU();
                await gpuBackend.initialise();
                setGpu(gpuBackend);
            } catch (error) {
                console.error("Failed to initialise SAXS calculator:", error);
            }
        }

        loadModule();
    }, []);

    // ---------------------------------------------------------
    // STRUCTURE FILE
    // ---------------------------------------------------------

    async function handleFile(event: React.ChangeEvent<HTMLInputElement>) {
        const file = event.target.files?.[0];
        if (!file) return;

        const text = await file.text();

        setFilename(file.name);
        setFileContents(text);
        setStructureFilename(file.name);
        setStructureText(text);

        // A new structure requires fresh buffers.
        setIsPrepared(false);
        setCalculatedCurve(null);
        setFitResult(null);
        setBestParameters(null);
        console.log("Structure loaded:", file.name);
        console.log(text.substring(0, 300));
    }

    // ---------------------------------------------------------
    // SAXS PROFILE FILE
    // ---------------------------------------------------------

    async function handleSaxsFile(event: React.ChangeEvent<HTMLInputElement>) {
        const file = event.target.files?.[0];
        if (!file) return;

        try {
            const profile = await parseSAXSProfile(file);

            setSaxsProfile(profile);

            // A new q grid requires fresh buffers.
            setIsPrepared(false);
            setCalculatedCurve(null);
            setFitResult(null);
            setBestParameters(null);
            console.log("SAXS profile loaded:", file.name);
            console.log(
                "q range (A^-1):",
                profile.q[0],
                "to",
                profile.q[profile.q.length - 1]
            );
            console.log("Data points:", profile.q.length);
        } catch (error) {
            console.error("Failed to load SAXS profile:", error);
        }
    }

    // ---------------------------------------------------------
    // PREPARE WASM + GPU BUFFERS
    //
    // Called only when the user clicks the Calculate button.
    // Returns true if preparation succeeds.
    // ---------------------------------------------------------

    async function prepareBuffers(): Promise<boolean> {
        if (
            !module ||
            !gpu ||
            !structureText ||
            !structureFilename ||
            !saxsProfile ||
            !elementJSON ||
            !saxsAtomTypeJSON ||
            !residueAtomTypeJSON
        ) {
            console.error("Cannot prepare: required data is missing.");
            return false;
        }

        if (isPreparing) {
            return false;
        }

        if (isPrepared) {
            return true;
        }

        setIsPreparing(true);

        const qValues = new module.FloatVector();

        try {
            console.log("===== PREPARING SAXS GPU BUFFERS =====");

            for (const q of saxsProfile.q) {
                qValues.push_back(q);
            }

            console.log("Preparing amplitudes for", qValues.size(), "q values");

            module.prepare_debye_gpu(
                elementJSON,
                saxsAtomTypeJSON,
                residueAtomTypeJSON,
                structureText,
                structureFilename,
                qValues
            );

            const posVec = module.get_positions();
            const amplitudeVec = module.get_amplitudes();

            try {
                const positions = new Float32Array(posVec.size());

                for (let i = 0; i < posVec.size(); i++) {
                    positions[i] = posVec.get(i);
                }

                const amplitudes = new Float32Array(amplitudeVec.size());

                for (let i = 0; i < amplitudeVec.size(); i++) {
                    amplitudes[i] = amplitudeVec.get(i);
                }

                const expectedAmplitudeCount =
                    saxsProfile.q.length * (positions.length / 3);

                if (amplitudes.length !== expectedAmplitudeCount) {
                    throw new Error(
                        `Amplitude buffer has wrong size. Expected ` +
                        `${expectedAmplitudeCount}, got ${amplitudes.length}.`
                    );
                }

                gpu.uploadBuffers(positions, amplitudes);

                setAtomCount(positions.length / 3);
                setIsPrepared(true);

                console.log("===== GPU BUFFERS READY =====");
                console.log("Atoms:", positions.length / 3);
                console.log("q values:", saxsProfile.q.length);
                console.log("Amplitudes:", amplitudes.length);

                return true;
            } finally {
                posVec.delete();
                amplitudeVec.delete();
            }
        } catch (error) {
            console.error("Failed to prepare GPU buffers:", error);
            return false;
        } finally {
            qValues.delete();
            setIsPreparing(false);
        }
    }

    // ---------------------------------------------------------
    // MULTIPOLE L VALIDATION
    // ---------------------------------------------------------

    function validateMultipoleL(debyeIntensities: number[]) {
        if (!module) {
            console.error("WASM module is not ready.");
            return;
        }

        if (!saxsProfile) {
            console.error("No SAXS profile loaded.");
            return;
        }

        const testLValues = [2, 4, 6, 8, 10, 12];

        console.log("========================================");
        console.log("MULTIPOLE L VALIDATION");
        console.log("Reference: GPU Debye calculation");
        console.log("q points:", debyeIntensities.length);
        console.log("========================================");
        console.log(
            "L | coefficients | max rel. error | RMS rel. error | q at max error"
        );
        console.log(
            "---------------------------------------------------------------------"
        );

        for (const L of testLValues) {
            const vector = module.get_multipole_intensities(L);

            try {
                const multipoleIntensities = new Float32Array(vector.size());

                for (let i = 0; i < vector.size(); i++) {
                    multipoleIntensities[i] = vector.get(i);
                }

                if (multipoleIntensities.length !== debyeIntensities.length) {
                    console.error(
                        `L=${L}: intensity array size mismatch. ` +
                        `Multipole=${multipoleIntensities.length}, ` +
                        `Debye=${debyeIntensities.length}`
                    );
                    continue;
                }

                let maxRelativeError = 0;
                let maxErrorIndex = 0;
                let sumSquaredRelativeError = 0;
                let validPoints = 0;

                for (let i = 0; i < debyeIntensities.length; i++) {
                    const debye = debyeIntensities[i];
                    const multipole = multipoleIntensities[i];

                    if (Math.abs(debye) < 1e-20) continue;

                    const relativeError =
                        Math.abs(multipole - debye) / Math.abs(debye);

                    if (relativeError > maxRelativeError) {
                        maxRelativeError = relativeError;
                        maxErrorIndex = i;
                    }

                    sumSquaredRelativeError += relativeError * relativeError;
                    validPoints++;
                }

                const rmsRelativeError =
                    validPoints > 0
                        ? Math.sqrt(sumSquaredRelativeError / validPoints)
                        : 0;

                const qAtMaxError = saxsProfile.q[maxErrorIndex];
                const coefficientCount = (L + 1) * (L + 1);

                console.log(
                    `${L} | ${coefficientCount} | ${maxRelativeError} | ` +
                    `${rmsRelativeError} | ${qAtMaxError} Å^-1`
                );
            } finally {
                vector.delete();
            }
        }

        console.log("========================================");
    }

    // ---------------------------------------------------------
    // CALCULATE SAXS CURVE
    // ---------------------------------------------------------

    async function calculateSAXSCurve() {
        if (!gpu) {
            console.error("GPU is not ready.");
            return;
        }

        if (!module) {
            console.error("WASM module is not ready.");
            return;
        }

        if (!saxsProfile) {
            console.error("No SAXS profile loaded.");
            return;
        }

        const MULTIPOLE_L = module.get_multipole_L();        
        const SOLVENT_DENSITY = 0.334;
        const SHELL_CONTRAST = 0.0167;

        try {
            const qValues = new Float32Array(saxsProfile.q);

            console.log("Calculating curve on experimental q grid:");
            console.log(qValues.length, "q values");
            console.log("q range:", qValues[0], "to", qValues[qValues.length - 1]);

            const batchSize = 100;
            const allQValues: number[] = [];
            const allIntensities: number[] = [];

            const numberOfBatches = Math.ceil(qValues.length / batchSize);

            console.log("GPU batch size:", batchSize);
            console.log("Number of GPU batches:", numberOfBatches);

            for (let start = 0; start < qValues.length; start += batchSize) {
                const end = Math.min(start + batchSize, qValues.length);
                const batch = qValues.slice(start, end);
                const batchNumber = Math.floor(start / batchSize) + 1;

                console.log(
                    `Calculating GPU batch ${batchNumber}/${numberOfBatches}: ` +
                    `q values ${start}–${end - 1}`
                );

                const batchCurve = await gpu.calculateCurve(batch, start);

                allQValues.push(...batchCurve.qValues);
                allIntensities.push(...batchCurve.intensities);

                console.log(`GPU batch ${batchNumber}/${numberOfBatches} complete`);
            }

            const gpuCurve: SAXSCurve = {
                qValues: allQValues,
                intensities: allIntensities
            };

            // -------------------------------------------------
            // COMPARE ATOMIC, EXCLUDED-VOLUME AND HYDRATION MODELS
            // -------------------------------------------------

            const atomicVector =
                module.get_multipole_intensities(MULTIPOLE_L);

            const combinedZeroVector =
                module.get_combined_multipole_intensities(0.0);

            const combinedSolventVector =
                module.get_combined_multipole_intensities(SOLVENT_DENSITY);

            const hydrationCombinedVector =
                module.get_hydration_combined_multipole_intensities(
                    SOLVENT_DENSITY,
                    SHELL_CONTRAST,
                    0.0
                );

            let combinedHydrationCurve: SAXSCurve | null = null;
            let bestHydrationCurve: SAXSCurve | null = null;
            let bestHydrationFit: SAXSFitResult | null = null;

            let bestRadiusAdjustmentFraction = 0;
            let bestShellContrast = 0;
            try {
                const atomic = Array.from(
                    { length: atomicVector.size() },
                    (_, i) => atomicVector.get(i)
                );

                const combinedZero = Array.from(
                    { length: combinedZeroVector.size() },
                    (_, i) => combinedZeroVector.get(i)
                );

                const combinedSolvent = Array.from(
                    { length: combinedSolventVector.size() },
                    (_, i) => combinedSolventVector.get(i)
                );

                const hydrationCombined = Array.from(
                    { length: hydrationCombinedVector.size() },
                    (_, i) => hydrationCombinedVector.get(i)
                );

                const expectedLength = saxsProfile.q.length;

                const lengths = [
                    atomic.length,
                    combinedZero.length,
                    combinedSolvent.length,
                    hydrationCombined.length,
                    gpuCurve.intensities.length
                ];

                if (lengths.some(length => length !== expectedLength)) {
                    throw new Error(
                        `Intensity array size mismatch. Expected ${expectedLength}; ` +
                        `got atomic=${atomic.length}, ` +
                        `zero-density=${combinedZero.length}, ` +
                        `excluded-volume=${combinedSolvent.length}, ` +
                        `hydration=${hydrationCombined.length}, ` +
                        `GPU=${gpuCurve.intensities.length}.`
                    );
                }

                const atomicCurve: SAXSCurve = {
                    qValues: [...saxsProfile.q],
                    intensities: atomic
                };

                const excludedVolumeCurve: SAXSCurve = {
                    qValues: [...saxsProfile.q],
                    intensities: combinedSolvent
                };

                const hydrationCurve: SAXSCurve = {
                    qValues: [...saxsProfile.q],
                    intensities: hydrationCombined
                };

                combinedHydrationCurve = hydrationCurve;

                const atomicFit = fitSAXSScale(saxsProfile, atomicCurve);
                const excludedVolumeFit = fitSAXSScale(
                    saxsProfile,
                    excludedVolumeCurve
                );
                const hydrationFit = fitSAXSScale(saxsProfile, hydrationCurve);

                console.log("========== SAXS MODEL FIT COMPARISON ==========");

                console.log("GPU Debye reference:");
                console.log("  Scale:", fitResult?.scale);
                console.log("  Reduced chi-squared:", fitResult?.reducedChiSquared);

                console.log(`Atomic multipoles, L=${MULTIPOLE_L}:`);
                console.log("  Scale:", atomicFit.scale);
                console.log("  Reduced chi-squared:", atomicFit.reducedChiSquared);

                console.log(
                    `Atomic + excluded volume, L=${MULTIPOLE_L}, ` +
                    `rho=${SOLVENT_DENSITY}:`
                );
                console.log("  Scale:", excludedVolumeFit.scale);
                console.log(
                    "  Reduced chi-squared:",
                    excludedVolumeFit.reducedChiSquared
                );

                console.log(
                    `Atomic + excluded volume + hydration, L=${MULTIPOLE_L}, ` +
                    `rho=${SOLVENT_DENSITY}, shell contrast=${SHELL_CONTRAST}:`
                );
                console.log("  Scale:", hydrationFit.scale);
                console.log("  Reduced chi-squared:", hydrationFit.reducedChiSquared);

                console.log("First five intensity comparisons:");

                for (let i = 0; i < Math.min(5, atomic.length); i++) {
                    console.log({
                        q: saxsProfile.q[i],
                        atomic: atomic[i],
                        excludedVolume: combinedSolvent[i],
                        hydration: hydrationCombined[i],
                        gpuDebye: gpuCurve.intensities[i]
                    });
                }

                console.log("================================================");

 
                // -------------------------------------------------
                // EXHAUSTIVE GRID SEARCH: RADIUS ADJUSTMENT + SHELL CONTRAST
                // -------------------------------------------------

                const GRID_SIZE = 100;

                const radiusAdjustmentValues = Array.from(
                    { length: GRID_SIZE },
                    (_, i) => -0.05 + i * (0.10 / (GRID_SIZE - 1))
                );

                const shellContrastValues = Array.from(
                    { length: GRID_SIZE },
                    (_, i) => i * (0.0334 / (GRID_SIZE - 1))
                );

                let bestReducedChiSquared = Infinity;
                let evaluatedPairs = 0;

                console.log("========== TWO-PARAMETER GRID SEARCH ==========");
                console.log("Radius adjustment fractions:", radiusAdjustmentValues);
                console.log("Shell contrast range: 0 to 0.0334 e/A^3");
                console.log("Radius candidates:", radiusAdjustmentValues.length);
                console.log("Contrast candidates:", shellContrastValues.length);
                console.log(
                    "Total candidate pairs:",
                    radiusAdjustmentValues.length * shellContrastValues.length
                );

                for (const radiusAdjustmentFraction of radiusAdjustmentValues) {
                    for (const shellContrast of shellContrastValues) {
                        const candidateVector =
                            module.get_hydration_combined_multipole_intensities(
                                SOLVENT_DENSITY,
                                shellContrast,
                                radiusAdjustmentFraction
                            );

                        try {
                            const intensities = Array.from(
                                { length: candidateVector.size() },
                                (_, i) => candidateVector.get(i)
                            );

                            if (intensities.length !== saxsProfile.q.length) {
                                throw new Error(
                                    `Grid candidate returned ${intensities.length} ` +
                                    `intensities; expected ${saxsProfile.q.length}.`
                                );
                            }

                            const candidateCurve: SAXSCurve = {
                                qValues: [...saxsProfile.q],
                                intensities
                            };

                            // fitSAXSScale optimises the multiplicative scale
                            // for this candidate before comparing its reduced chi-squared.
                            const candidateFit = fitSAXSScale(
                                saxsProfile,
                                candidateCurve
                            );

                            evaluatedPairs++;

                            if (
                                Number.isFinite(candidateFit.reducedChiSquared) &&
                                candidateFit.reducedChiSquared < bestReducedChiSquared
                            ) {
                                bestReducedChiSquared =
                                    candidateFit.reducedChiSquared;

                                bestRadiusAdjustmentFraction =
                                    radiusAdjustmentFraction;

                                bestShellContrast = shellContrast;

                                bestHydrationCurve = candidateCurve;
                                bestHydrationFit = candidateFit;
                            }
                        } finally {
                            candidateVector.delete();
                        }
                    }
                }

                if (!bestHydrationCurve || !bestHydrationFit) {
                    throw new Error(
                        "Grid search failed to produce a valid candidate fit."
                    );
                }

                combinedHydrationCurve = bestHydrationCurve;

                setBestParameters({
                    radiusAdjustmentFraction: bestRadiusAdjustmentFraction,
                    shellContrast: bestShellContrast
                });

                console.log("Grid search complete.");
                console.log("Evaluated parameter pairs:", evaluatedPairs);
                console.log("Best radius adjustment fraction:",
                    bestRadiusAdjustmentFraction);
                console.log("Best shell contrast (e/A^3):", bestShellContrast);
                console.log("Best fitted scale:", bestHydrationFit.scale);
                console.log(
                    "Best chi-squared:",
                    bestHydrationFit.chiSquared
                );
                console.log(
                    "Best reduced chi-squared:",
                    bestHydrationFit.reducedChiSquared
                );
                console.log("===============================================");


                // Verify rho=0 removes the excluded-volume term.
                let maxZeroDifference = 0;

                for (let i = 0; i < atomic.length; i++) {
                    const difference =
                        Math.abs(combinedZero[i] - atomic[i]) /
                        Math.max(Math.abs(atomic[i]), 1e-30);

                    maxZeroDifference = Math.max(maxZeroDifference, difference);
                }

                console.log("Excluded-volume sanity check:");
                console.log(
                    "  rho=0 maximum relative difference from atomic-only:",
                    maxZeroDifference
                );
            } finally {
                atomicVector.delete();
                combinedZeroVector.delete();
                combinedSolventVector.delete();
                hydrationCombinedVector.delete();
            }

            validateMultipoleL(gpuCurve.intensities);

            console.log("===== COMPLETE GPU CURVE =====");
            console.log("Calculated points:", gpuCurve.qValues.length);
            console.log("Expected points:", saxsProfile.q.length);
            console.log("==============================");

            if (!combinedHydrationCurve) {
                throw new Error("Combined hydration curve was not generated.");
            }

            const productionCurve = combinedHydrationCurve;
            setCalculatedCurve(productionCurve);

            const result = fitSAXSScale(saxsProfile, productionCurve);
            setFitResult(result);

            console.log("===== SAXS SCALE FIT =====");
            console.log("Scale:", result.scale);
            console.log("Chi-squared:", result.chiSquared);
            console.log("Reduced chi-squared:", result.reducedChiSquared);
            console.log("==========================");

            for (
                let i = 0;
                i < Math.min(10, productionCurve.qValues.length);
                i++
            ) {
                console.log(
                    "Production q =",
                    productionCurve.qValues[i],
                    "I =",
                    productionCurve.intensities[i]
                );
            }

            console.log("Production SAXS curve calculation finished.");
            console.log("Production model: atomic + excluded volume + hydration");
            console.log("Shell contrast:", SHELL_CONTRAST);
        } catch (error) {
            console.error("Failed to calculate SAXS curve:", error);
        }
    }

    // ---------------------------------------------------------
    // EXPORT
    // ---------------------------------------------------------

    function exportSAXSData() {
        if (!saxsProfile || !calculatedCurve || !fitResult) {
            console.error("No SAXS data available to export.");
            return;
        }

        const rows = [
            [
                "q_A^-1",
                "I_exp",
                "sigma",
                "I_calc_raw",
                "I_calc_scaled"
            ].join(",")
        ];

        for (let i = 0; i < saxsProfile.q.length; i += 50) {
            rows.push([
                saxsProfile.q[i],
                saxsProfile.intensity[i],
                saxsProfile.sigma[i],
                calculatedCurve.intensities[i],
                fitResult.scaledIntensities[i]
            ].join(","));
        }

        const blob = new Blob([rows.join("\n")], {
            type: "text/csv;charset=utf-8;"
        });

        const url = URL.createObjectURL(blob);
        const link = document.createElement("a");

        link.href = url;
        link.download = "saxs_comparison.csv";
        link.click();

        URL.revokeObjectURL(url);
    }

    // ---------------------------------------------------------
    // TEST DEBYE GPU
    // ---------------------------------------------------------

    async function testDebyeGPU() {
        if (!gpu) {
            console.error("GPU is not ready.");
            return;
        }

        console.log("===== DEBYE GPU TEST =====");

        const positions = new Float32Array([
            0, 0, 0,
            1, 0, 0,
            0, 1, 0
        ]);

        const atomCount = 3;

        const qValues = new Float32Array([
            0.0,
            0.1,
            1.0
        ]);

        const amplitudes = new Float32Array([
            1, 2, 3,
            2, 4, 6,
            3, 6, 9
        ]);

        gpu.uploadBuffers(positions, amplitudes);

        const gpuResult = await gpu.calculateCurve(qValues, 0);

        const cpuResult = calculateDebyeReference(
            positions,
            amplitudes,
            qValues,
            atomCount
        );

        console.log("===== CPU vs GPU =====");

        for (let i = 0; i < qValues.length; i++) {
            const cpu = cpuResult.intensities[i];
            const gpuValue = gpuResult.intensities[i];
            const difference = gpuValue - cpu;
            const relativeError = cpu !== 0 ? difference / cpu : 0;

            console.log(
                `q=${qValues[i]} CPU=${cpu} GPU=${gpuValue} ` +
                `difference=${difference} relative error=${relativeError}`
            );
        }

        console.log("======================");
        console.log("===== qOffset TEST =====");
    }

    // ---------------------------------------------------------
    // UI
    // ---------------------------------------------------------

    return (
        <div style={{ padding: "2rem", fontFamily: "sans-serif" }}>
            <h1>SAXS Profile Calculator</h1>

            <p>
                WASM status: {module ? "Loaded ✓" : "Loading..."}
            </p>

            <h2>Structure</h2>

            <input
                type="file"
                accept=".pdb,.cif,.mmcif,.ent"
                onChange={handleFile}
            />

            <h2>Experimental SAXS profile</h2>

            <input
                type="file"
                accept=".dat,.txt,.csv"
                onChange={handleSaxsFile}
            />

            {saxsProfile && (
                <div>
                    <h3>SAXS profile loaded</h3>

                    <p>Data points: {saxsProfile.q.length}</p>

                    <p>
                        q range: {saxsProfile.q[0]} to{" "}
                        {saxsProfile.q[saxsProfile.q.length - 1]} Å⁻¹
                    </p>
                </div>
            )}

            {gpu && saxsProfile && structureText && (
                <button
                    onClick={async () => {
                        if (isPreparing) return;

                        if (!isPrepared) {
                            const success = await prepareBuffers();

                            if (!success) return;
                        }

                        await calculateSAXSCurve();
                    }}
                    disabled={isPreparing}
                    style={{
                        marginTop: "1rem",
                        padding: "0.75rem 1.5rem",
                        fontSize: "1rem"
                    }}
                >
                    {isPreparing
                        ? "Preparing model..."
                        : isPrepared
                            ? "Calculate SAXS curve"
                            : "Prepare model and calculate SAXS curve"}
                </button>
            )}

            {calculatedCurve && fitResult && (
                <button onClick={exportSAXSData}>
                    Export SAXS data
                </button>
            )}

            {calculatedCurve && (
                <div>
                    <h3>Calculated SAXS curve</h3>

                    <p>Data points: {calculatedCurve.qValues.length}</p>

                    <p>
                        q range: {calculatedCurve.qValues[0]} to{" "}
                        {calculatedCurve.qValues[
                            calculatedCurve.qValues.length - 1
                        ]} Å⁻¹
                    </p>
                </div>
            )}

            {fitResult && (
                <div>
                    <h3>SAXS fit</h3>
                    <p>Scale factor: {fitResult.scale}</p>
                    <p>χ²: {fitResult.chiSquared}</p>
                    <p>Reduced χ²: {fitResult.reducedChiSquared}</p>

                    {bestParameters && (
                        <>
                            <h4>Best-fit model parameters</h4>
                            <p>
                                Radius adjustment fraction (δr/rₘ):{" "}
                                {bestParameters.radiusAdjustmentFraction.toFixed(4)}
                            </p>
                            <p>
                                Shell contrast (e/Å³):{" "}
                                {bestParameters.shellContrast.toFixed(6)}
                            </p>
                        </>
                    )}
                </div>
            )}

            {calculatedCurve && saxsProfile && (
                <SAXSPlot
                    experimentalQ={saxsProfile.q}
                    experimentalIntensity={saxsProfile.intensity}
                    calculatedQ={calculatedCurve.qValues}
                    calculatedIntensity={calculatedCurve.intensities}
                    scaledIntensity={fitResult?.scaledIntensities}
                />
            )}

            <button
                onClick={testDebyeGPU}
                style={{
                    marginTop: "1rem",
                    padding: "0.75rem 1.5rem",
                    fontSize: "1rem"
                }}
            >
                Test Debye GPU
            </button>

            {filename && (
                <div>
                    <h2>{filename}</h2>

                    <p>File size: {fileContents.length} characters</p>

                    <p>Atom count: {atomCount ?? "Unknown"}</p>

                    <h3>Atom summary</h3>

                    <pre
                        style={{
                            background: "#eee",
                            padding: "1rem",
                            overflowX: "auto"
                        }}
                    >
                        {summary}
                    </pre>

                    <textarea
                        value={fileContents}
                        readOnly
                        rows={20}
                        cols={100}
                    />
                </div>
            )}
        </div>
    );
}

export default App;
