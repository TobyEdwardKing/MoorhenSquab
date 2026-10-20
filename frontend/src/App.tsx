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

    const [saxsProfile, setSaxsProfile] =
        useState<SAXSProfile | null>(null);

    const [calculatedCurve, setCalculatedCurve] =
        useState<SAXSCurve | null>(null);

    const [module, setModule] =
        useState<any>(null);

    const [filename, setFilename] =
        useState("");

    const [fileContents, setFileContents] =
        useState("");

    const [atomCount, setAtomCount] =
        useState<number | null>(null);

    const [summary, setSummary] =
        useState("");

    const [elementJSON, setElementJSON] =
        useState("");

    const [saxsAtomTypeJSON, setSaxsAtomTypeJSON] =
        useState("");

    const [residueAtomTypeJSON, setResidueAtomTypeJSON] =
        useState("");

    const [gpu, setGpu] =
        useState<DebyeGPU | null>(null);

    const [fitResult, setFitResult] =
        useState<SAXSFitResult | null>(null);


    // ---------------------------------------------------------
    // STRUCTURE DATA
    // ---------------------------------------------------------

    const [structureText, setStructureText] =
        useState("");

    const [structureFilename, setStructureFilename] =
        useState("");


    // ---------------------------------------------------------
    // INITIALISE WASM + GPU
    // ---------------------------------------------------------

    useEffect(() => {

        async function loadModule() {

            try {

                const wasm =
                    await createSaxsModule({

                        locateFile: (path: string) => {

                            if (path.endsWith(".wasm")) {
                                return "/saxs.wasm";
                            }

                            return path;
                        }

                    });


                const elementResponse =
                    await fetch("/elementsInfo.json");

                const saxsAtomTypeResponse =
                    await fetch("/SAXSAtomTypes.json");

                const residueAtomTypeResponse =
                    await fetch("/residueAtomTypes.json");


                const elementJson =
                    await elementResponse.text();

                const saxsAtomTypeJson =
                    await saxsAtomTypeResponse.text();

                const residueAtomTypeJson =
                    await residueAtomTypeResponse.text();


                setElementJSON(elementJson);
                setSaxsAtomTypeJSON(saxsAtomTypeJson);
                setResidueAtomTypeJSON(residueAtomTypeJson);

                setModule(wasm);


                console.log(
                    "elementJSON:",
                    elementJson.slice(0, 200)
                );

                console.log(
                    "saxsAtomTypeJSON:",
                    saxsAtomTypeJson.slice(0, 200)
                );

                console.log(
                    "residueAtomTypeJSON:",
                    residueAtomTypeJson.slice(0, 200)
                );

                console.log(wasm);
                console.log(Object.keys(wasm));
                console.log("WASM loaded.");


                const gpuBackend =
                    new DebyeGPU();

                await gpuBackend.initialise();

                setGpu(gpuBackend);

            } catch (error) {

                console.error(
                    "Failed to initialise SAXS calculator:",
                    error
                );
            }
        }


        loadModule();

    }, []);


    // ---------------------------------------------------------
    // STRUCTURE FILE
    // ---------------------------------------------------------

    async function handleFile(
        event: React.ChangeEvent<HTMLInputElement>
    ) {

        const file =
            event.target.files?.[0];

        if (!file) {
            return;
        }


        const text =
            await file.text();


        setFilename(file.name);
        setFileContents(text);

        setStructureFilename(file.name);
        setStructureText(text);


        console.log(
            "Structure loaded:",
            file.name
        );

        console.log(
            text.substring(0, 300)
        );
    }


    // ---------------------------------------------------------
    // SAXS PROFILE FILE
    // ---------------------------------------------------------

    async function handleSaxsFile(
        event: React.ChangeEvent<HTMLInputElement>
    ) {

        const file =
            event.target.files?.[0];

        if (!file) {
            return;
        }


        try {

            const profile =
                await parseSAXSProfile(file);


            setSaxsProfile(profile);


            console.log(
                "SAXS profile loaded:",
                file.name
            );


            console.log(
                "q range (A^-1):",
                profile.q[0],
                "to",
                profile.q[
                    profile.q.length - 1
                ]
            );


            console.log(
                "Data points:",
                profile.q.length
            );


        } catch (error) {

            console.error(
                "Failed to load SAXS profile:",
                error
            );
        }
    }


    // ---------------------------------------------------------
    // PREPARE WASM + GPU BUFFERS
    //
    // This happens once both the structure and SAXS profile
    // are available.
    // ---------------------------------------------------------

    useEffect(() => {

        if (!module) {
            return;
        }

        if (!gpu) {
            return;
        }

        if (!structureText) {
            return;
        }

        if (!structureFilename) {
            return;
        }

        if (!saxsProfile) {
            return;
        }

        if (!elementJSON) {
            return;
        }

        if (!saxsAtomTypeJSON) {
            return;
        }

        if (!residueAtomTypeJSON) {
            return;
        }


        async function prepareBuffers() {

            const qValues = new module.FloatVector();

            try {

                console.log(
                    "===== PREPARING SAXS GPU BUFFERS ====="
                );


                for (const q of saxsProfile.q) {
                    qValues.push_back(q);
                }


                const qCount = qValues.size();

                console.log(
                    "Preparing amplitudes for",
                    qCount,
                    "q values"
                );

                console.log(
                    "q range:",
                    qCount > 0 ? qValues.get(0) : undefined,
                    "to",
                    qCount > 0 ? qValues.get(qCount - 1) : undefined
                );


                // -------------------------------------------------
                // Generate positions + q-dependent amplitudes
                // in WASM.
                // -------------------------------------------------

                module.prepare_debye_gpu(
                    elementJSON,
                    saxsAtomTypeJSON,
                    residueAtomTypeJSON,
                    structureText,
                    structureFilename,
                    qValues
                );


                const posVec =
                    module.get_positions();

                const amplitudeVec =
                    module.get_amplitudes();

                try {

                    console.log(
                        "positions size:",
                        posVec.size()
                    );

                    console.log(
                        "amplitudes size:",
                        amplitudeVec.size()
                    );


                    // -------------------------------------------------
                    // Copy Embind vectors into typed arrays.
                    // -------------------------------------------------

                    const positions =
                        new Float32Array(
                            posVec.size()
                        );


                    for (
                        let i = 0;
                        i < posVec.size();
                        i++
                    ) {

                        positions[i] =
                            posVec.get(i);
                    }


                    const amplitudes =
                        new Float32Array(
                            amplitudeVec.size()
                        );


                    for (
                        let i = 0;
                        i < amplitudeVec.size();
                        i++
                    ) {

                        amplitudes[i] =
                            amplitudeVec.get(i);
                    }


                    console.log(
                        "positions:",
                        positions.length,
                        "floats",
                        positions.byteLength,
                        "bytes"
                    );


                    console.log(
                        "amplitudes:",
                        amplitudes.length,
                        "floats",
                        amplitudes.byteLength,
                        "bytes"
                    );


                    // -------------------------------------------------
                    // Sanity check.
                    // -------------------------------------------------

                    const expectedAmplitudeCount =
                        saxsProfile.q.length *
                        (positions.length / 3);


                    console.log(
                        "Expected amplitude count:",
                        expectedAmplitudeCount
                    );

                    console.log(
                        "Actual amplitude count:",
                        amplitudes.length
                    );


                    if (
                        amplitudes.length !==
                        expectedAmplitudeCount
                    ) {

                        throw new Error(
                            `Amplitude buffer has wrong size. ` +
                            `Expected ${expectedAmplitudeCount}, ` +
                            `got ${amplitudes.length}.`
                        );
                    }


                    // -------------------------------------------------
                    // Upload to GPU.
                    // -------------------------------------------------

                    gpu.uploadBuffers(
                        positions,
                        amplitudes
                    );


                    setAtomCount(
                        positions.length / 3
                    );


                    console.log(
                        "===== GPU BUFFERS READY ====="
                    );

                    console.log(
                        "Atoms:",
                        positions.length / 3
                    );

                    console.log(
                        "q values:",
                        saxsProfile.q.length
                    );

                    console.log(
                        "Amplitudes:",
                        amplitudes.length
                    );

                } finally {
                    posVec.delete();
                    amplitudeVec.delete();
                }

            } catch (error) {

                console.error(
                    "Failed to prepare GPU buffers:",
                    error
                );

            } finally {
                qValues.delete();
            }
        }


        prepareBuffers();

    }, [
        module,
        gpu,
        structureText,
        structureFilename,
        saxsProfile,
        elementJSON,
        saxsAtomTypeJSON,
        residueAtomTypeJSON
    ]);

    // ---------------------------------------------------------
    // MULTIPOLE L VALIDATION
    //
    // Compare multipole reconstruction against the existing
    // GPU Debye reference for L = 2,4,6,8,10,12.
    // ---------------------------------------------------------

    function validateMultipoleL(
        debyeIntensities: number[]
    ) {

        if (!module) {
            console.error(
                "WASM module is not ready."
            );
            return;
        }

        if (!saxsProfile) {
            console.error(
                "No SAXS profile loaded."
            );
            return;
        }

        const testLValues = [
            2,
            4,
            6,
            8,
            10,
            12
        ];

        console.log(
            "========================================"
        );

        console.log(
            "MULTIPOLE L VALIDATION"
        );

        console.log(
            "Reference: GPU Debye calculation"
        );

        console.log(
            "q points:",
            debyeIntensities.length
        );

        console.log(
            "========================================"
        );

        console.log(
            "L | coefficients | max rel. error | RMS rel. error | q at max error"
        );

        console.log(
            "---------------------------------------------------------------------"
        );


        for (const L of testLValues) {

            const vector =
                module.get_multipole_intensities(
                    L
                );

            try {

                const multipoleIntensities =
                    new Float32Array(
                        vector.size()
                    );


                for (
                    let i = 0;
                    i < vector.size();
                    i++
                ) {

                    multipoleIntensities[i] =
                        vector.get(i);
                }


                if (
                    multipoleIntensities.length !==
                    debyeIntensities.length
                ) {

                    console.error(
                        `L=${L}: intensity array size mismatch. ` +
                        `Multipole=${multipoleIntensities.length}, ` +
                        `Debye=${debyeIntensities.length}`
                    );

                    continue;
                }


                let maxRelativeError = 0.0;
                let maxErrorIndex = 0;

                let sumSquaredRelativeError = 0.0;
                let validPoints = 0;


                for (
                    let i = 0;
                    i < debyeIntensities.length;
                    i++
                ) {

                    const debye =
                        debyeIntensities[i];

                    const multipole =
                        multipoleIntensities[i];


                    /*
                     * The Debye intensity should be positive.
                     *
                     * Ignore exactly-zero reference points rather
                     * than dividing by zero.
                     */

                    if (Math.abs(debye) < 1e-20) {
                        continue;
                    }


                    const relativeError =
                        Math.abs(
                            multipole - debye
                        ) /
                        Math.abs(debye);


                    if (
                        relativeError >
                        maxRelativeError
                    ) {

                        maxRelativeError =
                            relativeError;

                        maxErrorIndex =
                            i;
                    }


                    sumSquaredRelativeError +=
                        relativeError *
                        relativeError;

                    validPoints++;
                }


                const rmsRelativeError =
                    validPoints > 0
                        ? Math.sqrt(
                            sumSquaredRelativeError /
                            validPoints
                        )
                        : 0.0;


                const qAtMaxError =
                    saxsProfile.q[
                        maxErrorIndex
                    ];


                const coefficientCount =
                    (L + 1) * (L + 1);


                console.log(
                    `${L} | ` +
                    `${coefficientCount} | ` +
                    `${maxRelativeError} | ` +
                    `${rmsRelativeError} | ` +
                    `${qAtMaxError} Å^-1`
                );

            } finally {
                vector.delete();
            }
        }


        console.log(
            "========================================"
        );
    }

    // ---------------------------------------------------------
    // CALCULATE SAXS CURVE
    // ---------------------------------------------------------

    async function calculateSAXSCurve() {

        if (!gpu) {

            console.error(
                "GPU is not ready."
            );

            return;
        }


        if (!saxsProfile) {

            console.error(
                "No SAXS profile loaded."
            );

            return;
        }


        try {

            /*
             * Use the experimental q grid directly.
             *
             * SAXSProfile.q is already converted to
             * Angstrom^-1 by the parser.
             */

            const qValues =
                new Float32Array(
                    saxsProfile.q
                );


            console.log(
                "Calculating curve on experimental q grid:"
            );


            console.log(
                qValues.length,
                "q values"
            );


            console.log(
                "q range:",
                qValues[0],
                "to",
                qValues[
                    qValues.length - 1
                ]
            );


            /*
             * Calculate the curve in batches.
             *
             * A single dispatch containing all 1949 q values
             * caused the D3D12 device to hang. Batching keeps
             * each GPU dispatch below that workload.
             */

            const batchSize = 100;

            const allQValues: number[] = [];
            const allIntensities: number[] = [];

            const numberOfBatches =
                Math.ceil(
                    qValues.length / batchSize
                );


            console.log(
                "GPU batch size:",
                batchSize
            );

            console.log(
                "Number of GPU batches:",
                numberOfBatches
            );


            for (
                let start = 0;
                start < qValues.length;
                start += batchSize
            ) {

                const end =
                    Math.min(
                        start + batchSize,
                        qValues.length
                    );


                const batch =
                    qValues.slice(
                        start,
                        end
                    );


                const batchNumber =
                    Math.floor(
                        start / batchSize
                    ) + 1;


                console.log(
                    `Calculating GPU batch ${batchNumber}/${numberOfBatches}: q values ${start}–${end - 1}`
                );

                const batchCurve =
                    await gpu.calculateCurve(
                        batch,
                        start
                    );


                allQValues.push(
                    ...batchCurve.qValues
                );


                allIntensities.push(
                    ...batchCurve.intensities
                );


                console.log(
                    `GPU batch ${batchNumber}/${numberOfBatches} complete`
                );
            }


            /*
             * Reconstruct the complete calculated curve.
             */

            const gpuCurve: SAXSCurve = {

                qValues:
                    allQValues,

                intensities:
                    allIntensities

            };


            // -------------------------------------------------
            // TEST COMBINED ATOMIC + EXCLUDED-VOLUME INTENSITY
            // -------------------------------------------------

            const atomicVector =
                module.get_multipole_intensities(12);

            const combinedZeroVector =
                module.get_combined_multipole_intensities(0.0);

            const combinedSolventVector =
                module.get_combined_multipole_intensities(0.334);

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

                if (
                    atomic.length !== combinedZero.length ||
                    atomic.length !== combinedSolvent.length ||
                    atomic.length !== gpuCurve.intensities.length
                ) {
                    throw new Error(
                        "Atomic, combined, and GPU intensity arrays must have matching lengths"
                    );
                }

                let maxZeroDifference = 0;
                let maxSolventDifference = 0;

                let sumZeroSquaredDifference = 0;
                let sumSolventSquaredDifference = 0;

                for (let i = 0; i < atomic.length; i++) {
                    const zeroDifference =
                        Math.abs(combinedZero[i] - atomic[i]) /
                        Math.max(Math.abs(atomic[i]), 1e-30);

                    const solventDifference =
                        Math.abs(combinedSolvent[i] - atomic[i]) /
                        Math.max(Math.abs(atomic[i]), 1e-30);

                    maxZeroDifference = Math.max(
                        maxZeroDifference,
                        zeroDifference
                    );

                    maxSolventDifference = Math.max(
                        maxSolventDifference,
                        solventDifference
                    );

                    sumZeroSquaredDifference +=
                        zeroDifference * zeroDifference;

                    sumSolventSquaredDifference +=
                        solventDifference * solventDifference;
                }
                // Compare fits without changing the displayed GPU curve.

                const atomicCurve: SAXSCurve = {
                    qValues: [...saxsProfile.q],
                    intensities: atomic
                };

                const combinedCurve: SAXSCurve = {
                    qValues: [...saxsProfile.q],
                    intensities: combinedSolvent
                };

                const atomicFit = fitSAXSScale(
                    saxsProfile,
                    atomicCurve
                );

                const combinedFit = fitSAXSScale(
                    saxsProfile,
                    combinedCurve
                );

                console.log("===== FIT COMPARISON =====");

                console.log("GPU Debye:");
                console.log("  Scale:", fitResult?.scale);
                console.log("  Reduced chi-squared:", fitResult?.reducedChiSquared);

                console.log("Atomic multipoles, L=12:");
                console.log("  Scale:", atomicFit.scale);
                console.log("  Reduced chi-squared:", atomicFit.reducedChiSquared);

                console.log("Atomic + excluded volume, rho=0.334:");
                console.log("  Scale:", combinedFit.scale);
                console.log("  Reduced chi-squared:", combinedFit.reducedChiSquared);

                console.log("First five atomic/combined intensity ratios:");

                for (let i = 0; i < Math.min(5, atomic.length); i++) {
                    console.log({
                        q: saxsProfile.q[i],
                        atomic: atomic[i],
                        combined: combinedSolvent[i],
                        combinedToAtomicRatio:
                            combinedSolvent[i] / Math.max(atomic[i], 1e-30)
                    });
                }

                console.log("==========================");
                console.log("===== EXCLUDED-VOLUME TEST =====");

                console.log(
                    "rho = 0: maximum relative difference from atomic-only:",
                    maxZeroDifference
                );

                console.log(
                    "rho = 0: RMS relative difference from atomic-only:",
                    Math.sqrt(
                        sumZeroSquaredDifference / atomic.length
                    )
                );

                console.log(
                    "rho = 0.334: maximum relative difference from atomic-only:",
                    maxSolventDifference
                );

                console.log(
                    "rho = 0.334: RMS relative difference from atomic-only:",
                    Math.sqrt(
                        sumSolventSquaredDifference / atomic.length
                    )
                );

                console.log(
                    "First five values:",
                    Array.from(
                        { length: Math.min(5, atomic.length) },
                        (_, i) => ({
                            q: saxsProfile.q[i],
                            atomic: atomic[i],
                            excludedVolumeAtZeroDensity: combinedZero[i],
                            atomicPlusExcludedVolume: combinedSolvent[i],
                            gpuDebye: gpuCurve.intensities[i]
                        })
                    )
                );

                console.log("================================");

            } finally {
                atomicVector.delete();
                combinedZeroVector.delete();
                combinedSolventVector.delete();
            }


            /*
             * -------------------------------------------------
             * Validate multipole expansion against the Debye
             * reference.
             *
             * The GPU Debye curve is already calculated, so
             * use it directly as the reference rather than
             * performing another Debye calculation.
             * -------------------------------------------------
             */

            validateMultipoleL(
                gpuCurve.intensities
            );


            console.log(
                "===== COMPLETE GPU CURVE ====="
            );

            console.log(
                "Calculated points:",
                gpuCurve.qValues.length
            );

            console.log(
                "Expected points:",
                saxsProfile.q.length
            );

            console.log(
                "=============================="
            );


            setCalculatedCurve(
                gpuCurve
            );


            /*
             * Fit calculated curve to experimental
             * intensities using a multiplicative scale.
             */

            const result =
                fitSAXSScale(
                    saxsProfile,
                    gpuCurve
                );


            setFitResult(
                result
            );


            console.log(
                "===== SAXS SCALE FIT ====="
            );

            console.log(
                "Scale:",
                result.scale
            );

            console.log(
                "Chi-squared:",
                result.chiSquared
            );

            console.log(
                "Reduced chi-squared:",
                result.reducedChiSquared
            );

            console.log(
                "=========================="
            );


            /*
             * Print first few calculated values.
             */

            for (
                let i = 0;
                i < Math.min(
                    10,
                    gpuCurve.qValues.length
                );
                i++
            ) {

                console.log(
                    "q =",
                    gpuCurve.qValues[i],
                    "I =",
                    gpuCurve.intensities[i]
                );
            }


            console.log(
                "GPU curve calculation finished."
            );


        } catch (error) {

            console.error(
                "Failed to calculate SAXS curve:",
                error
            );
        }
    }


    // ---------------------------------------------------------
    // EXPORT
    // ---------------------------------------------------------

    function exportSAXSData() {

        if (
            !saxsProfile ||
            !calculatedCurve ||
            !fitResult
        ) {

            console.error(
                "No SAXS data available to export."
            );

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


        for (
            let i = 0;
            i < saxsProfile.q.length;
            i += 50
        ) {

            rows.push([

                saxsProfile.q[i],
                saxsProfile.intensity[i],
                saxsProfile.sigma[i],
                calculatedCurve.intensities[i],
                fitResult.scaledIntensities[i]

            ].join(","));
        }


        const blob =
            new Blob(
                [rows.join("\n")],
                {
                    type:
                        "text/csv;charset=utf-8;"
                }
            );


        const url =
            URL.createObjectURL(blob);

        const link =
            document.createElement("a");


        link.href = url;

        link.download =
            "saxs_comparison.csv";

        link.click();


        URL.revokeObjectURL(url);
    }

    async function testDebyeGPU() {

        if (!gpu) {
            console.error("GPU is not ready.");
            return;
        }

        console.log("===== DEBYE GPU TEST =====");

        // Three atoms:
        //
        // atom 0 = (0, 0, 0)
        // atom 1 = (1, 0, 0)
        // atom 2 = (0, 1, 0)

        const positions = new Float32Array([
            0, 0, 0,
            1, 0, 0,
            0, 1, 0
        ]);

        const atomCount = 3;

        // Three q values
        const qValues = new Float32Array([
            0.0,
            0.1,
            1.0
        ]);

        // q-dependent amplitudes.
        //
        // q = 0.0:
        //   [1, 2, 3]
        //
        // q = 0.1:
        //   [1, 2, 3]
        //
        // q = 1.0:
        //   [1, 2, 3]
        //
        // Flattened as q-major:
        //
        // [q0 atom0, q0 atom1, q0 atom2,
        //  q1 atom0, q1 atom1, q1 atom2,
        //  q2 atom0, q2 atom1, q2 atom2]

        const amplitudes = new Float32Array([
            1, 2, 3,
            2, 4, 6,
            3, 6, 9
        ]);

        // Upload our tiny artificial system.
        gpu.uploadBuffers(
            positions,
            amplitudes
        );

        // Calculate with GPU.
        const gpuResult =
            await gpu.calculateCurve(
                qValues,
                0
            );

        // Calculate independently on CPU.
        const cpuResult =
            calculateDebyeReference(
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

            const relativeError =
                cpu !== 0
                    ? difference / cpu
                    : 0;

            console.log(
                `q=${qValues[i]} ` +
                `CPU=${cpu} ` +
                `GPU=${gpuValue} ` +
                `difference=${difference} ` +
                `relative error=${relativeError}`
            );
        }

        console.log("======================");
        console.log("===== qOffset TEST =====");


    }

    // ---------------------------------------------------------
    // UI
    // ---------------------------------------------------------

    return (

        <div
            style={{
                padding: "2rem",
                fontFamily: "sans-serif"
            }}
        >

            <h1>
                SAXS Profile Calculator
            </h1>


            <p>
                WASM status:
                {" "}
                {module
                    ? "Loaded ✓"
                    : "Loading..."}
            </p>


            <h2>
                Structure
            </h2>


            <input
                type="file"
                accept=".pdb,.cif,.mmcif,.ent"
                onChange={handleFile}
            />


            <h2>
                Experimental SAXS profile
            </h2>


            <input
                type="file"
                accept=".dat,.txt,.csv"
                onChange={handleSaxsFile}
            />


            {saxsProfile && (

                <div>

                    <h3>
                        SAXS profile loaded
                    </h3>


                    <p>
                        Data points:
                        {" "}
                        {saxsProfile.q.length}
                    </p>


                    <p>
                        q range:
                        {" "}
                        {saxsProfile.q[0]}
                        {" "}
                        to
                        {" "}
                        {saxsProfile.q[
                            saxsProfile.q.length - 1
                        ]}
                        {" "}
                        Å⁻¹
                    </p>

                </div>

            )}


            {gpu && saxsProfile && atomCount && (

                <button
                    onClick={calculateSAXSCurve}
                    style={{
                        marginTop: "1rem",
                        padding: "0.75rem 1.5rem",
                        fontSize: "1rem"
                    }}
                >
                    Calculate SAXS curve
                </button>

            )}


            {calculatedCurve && fitResult && (

                <button
                    onClick={exportSAXSData}
                >
                    Export SAXS data
                </button>
            )}


            {calculatedCurve && (

                <div>

                    <h3>
                        Calculated SAXS curve
                    </h3>


                    <p>
                        Data points:
                        {" "}
                        {calculatedCurve.qValues.length}
                    </p>


                    <p>
                        q range:
                        {" "}
                        {calculatedCurve.qValues[0]}
                        {" "}
                        to
                        {" "}
                        {calculatedCurve.qValues[
                            calculatedCurve.qValues.length - 1
                        ]}
                        {" "}
                        Å⁻¹
                    </p>

                </div>

            )}


            {fitResult && (

                <div>

                    <h3>
                        SAXS fit
                    </h3>


                    <p>
                        Scale factor:
                        {" "}
                        {fitResult.scale}
                    </p>

                    <p>
                        χ²:
                        {" "}
                        {fitResult.chiSquared}
                    </p>


                    <p>
                        Reduced χ²:
                        {" "}
                        {fitResult.reducedChiSquared}
                    </p>

                </div>

            )}


            {calculatedCurve && saxsProfile && (

                <SAXSPlot
                    experimentalQ={saxsProfile.q}
                    experimentalIntensity={saxsProfile.intensity}
                    calculatedQ={calculatedCurve.qValues}
                    calculatedIntensity={calculatedCurve.intensities}
                    scaledIntensity={
                        fitResult?.scaledIntensities
                    }
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

                    <h2>
                        {filename}
                    </h2>


                    <p>
                        File size:
                        {" "}
                        {fileContents.length}
                        {" "}
                        characters
                    </p>


                    <p>
                        Atom count:
                        {" "}
                        {atomCount ?? "Unknown"}
                    </p>


                    <h3>
                        Atom summary
                    </h3>


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