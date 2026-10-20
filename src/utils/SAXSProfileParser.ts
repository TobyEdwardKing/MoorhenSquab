import { SAXSProfile } from "../types/SAXSProfile";

export async function parseSAXSProfile(
    file: File
): Promise<SAXSProfile> {

    const text = await file.text();

    const q: number[] = [];
    const intensity: number[] = [];
    const sigma: number[] = [];

    const lines = text.split(/\r?\n/);

    for (const rawLine of lines) {

        const line = rawLine.trim();

        // Ignore empty lines
        if (!line) {
            continue;
        }

        // Ignore common SAXS comment/header lines
        if (
            line.startsWith("#") ||
            line.startsWith(";") ||
            line.startsWith("//")
        ) {
            continue;
        }

        // Support whitespace, tabs, or commas
        const columns = line
            .split(/[\s,]+/)
            .filter(value => value.length > 0);

        // Need at least q and I(q)
        if (columns.length < 2) {
            continue;
        }

        const qValue = Number(columns[0]);
        const intensityValue = Number(columns[1]);

        // Skip headers or malformed rows
        if (
            !Number.isFinite(qValue) ||
            !Number.isFinite(intensityValue)
        ) {
            continue;
        }

        // Third column is uncertainty, if present
        let sigmaValue = 1.0;

        if (columns.length >= 3) {

            const parsedSigma = Number(columns[2]);

            if (
                Number.isFinite(parsedSigma) &&
                parsedSigma > 0
            ) {
                sigmaValue = parsedSigma;
            }
        }

        q.push(qValue);
        intensity.push(intensityValue);
        sigma.push(sigmaValue);
    }

    if (q.length === 0) {
        throw new Error(
            `Could not find valid SAXS data in ${file.name}. ` +
            `Expected at least two numeric columns: q and I(q).`
        );
    }

    // =========================================================
    // ATSAS q/s unit heuristic
    //
    // If the file does not specify units, ATSAS uses the
    // maximum scattering-vector value as a heuristic:
    //
    //     qMax < 1.0  -> assume Å^-1
    //     qMax >= 1.0 -> assume nm^-1
    //
    // Our internal convention is always Å^-1.
    //
    // Therefore, if nm^-1 is inferred, convert:
    //
    //     nm^-1 -> Å^-1
    //     q / 10
    // =========================================================

    const qMax = Math.max(...q);

    let inferredUnit: "A^-1" | "nm^-1";

    if (qMax >= 1.0) {

        inferredUnit = "nm^-1";

        for (let i = 0; i < q.length; i++) {
            q[i] *= 0.1;
        }

    } else {

        inferredUnit = "A^-1";
    }

    console.log("===== Parsed SAXS profile =====");

    console.log(
        "File:",
        file.name
    );

    console.log(
        "Data points:",
        q.length
    );

    console.log(
        "Original q range:",
        qMax >= 1.0
            ? `${q[0] * 10} to ${q[q.length - 1] * 10} nm^-1`
            : `${q[0]} to ${q[q.length - 1]} A^-1`
    );

    console.log(
        "Inferred q unit:",
        inferredUnit
    );

    console.log(
        "Internal q range (A^-1):",
        q[0],
        "to",
        q[q.length - 1]
    );

    console.log(
        "First point:",
        {
            q: q[0],
            intensity: intensity[0],
            sigma: sigma[0]
        }
    );

    console.log("==============================");

    return {
        q,
        intensity,
        sigma
    };
}
