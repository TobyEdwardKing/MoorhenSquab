    export interface DebyeReferenceResult {
    qValues: number[];
    intensities: number[];
}

export function calculateDebyeReference(
    positions: Float32Array,
    amplitudes: Float32Array,
    qValues: Float32Array,
    atomCount: number
): DebyeReferenceResult {

    const intensities: number[] = [];

    for (let qIndex = 0; qIndex < qValues.length; qIndex++) {

        const q = qValues[qIndex];

        // Amplitudes are flattened as:
        // [q0 atom0, q0 atom1, ..., q0 atomN,
        //  q1 atom0, q1 atom1, ..., q1 atomN, ...]

        const amplitudeOffset = qIndex * atomCount;

        let total = 0.0;

        for (let i = 0; i < atomCount; i++) {

            const ix = i * 3;

            const xi = positions[ix];
            const yi = positions[ix + 1];
            const zi = positions[ix + 2];

            const fi = amplitudes[amplitudeOffset + i];

            for (let j = 0; j < atomCount; j++) {

                const jx = j * 3;

                const xj = positions[jx];
                const yj = positions[jx + 1];
                const zj = positions[jx + 2];

                const fj = amplitudes[amplitudeOffset + j];

                const dx = xi - xj;
                const dy = yi - yj;
                const dz = zi - zj;

                const r = Math.sqrt(
                    dx * dx +
                    dy * dy +
                    dz * dz
                );

                const qr = q * r;

                const sinc =
                    Math.abs(qr) < 1e-6
                        ? 1.0
                        : Math.sin(qr) / qr;

                total += fi * fj * sinc;
            }
        }

        intensities.push(total);
    }

    return {
        qValues: Array.from(qValues),
        intensities
    };
}