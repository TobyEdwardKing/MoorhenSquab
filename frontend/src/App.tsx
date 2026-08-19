import { useEffect, useState } from "react";
import createSaxsModule from "./wasm/saxs.js";
import { DebyeGPU } from "./webgpu/DebyeGPU";

function App() {

    const [module, setModule] = useState<any>(null);

    const [filename, setFilename] = useState("");
    const [fileContents, setFileContents] = useState("");
    const [atomCount, setAtomCount] = useState<number | null>(null);
    const [summary, setSummary] = useState("");
    // const [formFactorJSON, setFormFactorJSON] = useState("");
    // const [atomTypeJSON, setAtomTypeJSON] = useState("");
    const [elementJSON, setElementJSON] = useState("");
    const [saxsAtomTypeJSON, setSaxsAtomTypeJSON] = useState("");
    const [residueAtomTypeJSON, setResidueAtomTypeJSON] = useState("");
    const [gpu, setGpu] =
        useState<DebyeGPU|null>(null);

    useEffect(() => {

    async function loadModule() {

        const wasm = await createSaxsModule({

            locateFile: (path: string) => {

                if (path.endsWith(".wasm")) {
                    return "/saxs.wasm";
                }

                return path;
            }

        });


        // const elementResponse = await fetch("/elementsInfo.json");
        // const atomTypeResponse = await fetch("/residueAtomTypes.json");

        // const elementJson = await elementResponse.text();
        // const atomTypeJson = await atomTypeResponse.text();

        // setFormFactorJSON(elementJson);
        // setAtomTypeJSON(atomTypeJson);

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
        console.log("elementJSON:", elementJSON.slice(0, 200));
        console.log("saxsAtomTypeJSON:", saxsAtomTypeJSON.slice(0, 200));
        console.log("residueAtomTypeJSON:", residueAtomTypeJSON.slice(0, 200));
        console.log(wasm);
        console.log(Object.keys(wasm));
        console.log("WASM loaded.");
        const gpuBackend =
            new DebyeGPU();

        await gpuBackend.initialise();

        setGpu(gpuBackend);

    }

        loadModule();

    }, []);

    async function handleFile(event: React.ChangeEvent<HTMLInputElement>) {

        const file = event.target.files?.[0];

        if (!file || module === null)
            return;

        setFilename(file.name);

        const text = await file.text();

        setFileContents(text);
        console.log(file.name);
        console.log(fileContents.substring(0, 300));
        console.log(text.substring(0,300));
        // const nAtoms = module.atom_count(
        //     text,
        //     file.name
        // );
        // setAtomCount(nAtoms);

        // const summary =
        //     module.atom_summary(
        //         text,
        //         file.name);

        // console.log(summary);

        // setSummary(summary);
        
        // const test =
        //     module.form_factor_test(
        //         formFactorJSON,
        //         text,
        //         file.name
        //     );

        // console.log(test);
        // module.compute_debye_curve(
        //     formFactorJSON,
        //     text,
        //     file.name,
        //     0.0,   // qMin
        //     0.5,   // qMax
        //     0.005  // qStep
        // );
        console.log("ABOUT TO CALL prepare_debye_gpu");
        console.log("elementJSON first chars:", JSON.stringify(elementJSON.slice(0, 100)));
        console.log("saxsAtomTypeJSON first chars:", JSON.stringify(saxsAtomTypeJSON.slice(0, 100)));
        console.log("residueAtomTypeJSON first chars:", JSON.stringify(residueAtomTypeJSON.slice(0, 100)));

        module.prepare_debye_gpu(
            elementJSON,
            saxsAtomTypeJSON,
            residueAtomTypeJSON,
            text,
            file.name,
            0.1
        );
        const posVec = module.get_positions();
        const ffVec = module.get_amplitudes();
        console.log("positions size", posVec.size());
        console.log("form factors size", ffVec.size());
        const positions = new Float32Array(posVec.size());

        for (let i = 0; i < posVec.size(); i++)
        {
            positions[i] = posVec.get(i);
        }

        const formFactors = new Float32Array(ffVec.size());

        for (let i = 0; i < ffVec.size(); i++)
        {
            formFactors[i] = ffVec.get(i);
        }
        console.log(
        positions.length,
        positions.byteLength,
        formFactors.length,
        formFactors.byteLength
        );
        if (gpu)
            {
                gpu.uploadBuffers(
                    positions,
                    formFactors
                );

                console.log("Buffers uploaded to GPU");

                const gpuCurve =
                    await gpu.calculateCurve(
                        0.0,
                        0.5,
                        0.005
                    );

                console.log(
                    "GPU curve:",
                    gpuCurve
                );

                console.log(
                    "Curve contains",
                    gpuCurve.qValues.length,
                    "points"
                );

                for (let i = 0; i < gpuCurve.qValues.length; i++)
                {
                    console.log(
                        gpuCurve.qValues[i],
                        gpuCurve.intensities[i]
                    );

                    if (Math.abs(gpuCurve.qValues[i] - 0.1) < 1e-6)
                    {
                        console.log(
                            "GPU intensity at q = 0.1:",
                            gpuCurve.intensities[i]
                        );
                    }
                }

                console.log("GPU dispatch finished");
            }
}

    return (

        <div style={{ padding: "2rem", fontFamily: "sans-serif" }}>

            <h1>SAXS Profile Calculator</h1>

            <p>
                WASM status:
                {" "}
                {module ? "Loaded ✓" : "Loading..."}
            </p>

            <input
                type="file"
                accept=".pdb,.cif,.mmcif,.ent"
                onChange={handleFile}
            />

            {filename &&

                <div>

                    <h2>{filename}</h2>

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

            }

        </div>

    );

}

export default App;