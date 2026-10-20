import type { SAXSProfile } from "../types/SAXSProfile";

export interface SAXSComparison {

    scaleFactor: number;

    calculatedIntensity: number[];

    scaledIntensity: number[];

    residuals: number[];

    chiSquared: number;

    reducedChiSquared: number;
}


export function compareSAXSProfiles(
    experimental: SAXSProfile,
    calculatedIntensity: number[]
): SAXSComparison {

    if (
        experimental.q.length !==
        calculatedIntensity.length
    ) {
        throw new Error(
            "Experimental and calculated profiles " +
            "must contain the same number of points."
        );
    }


    let numerator = 0;
    let denominator = 0;


    for (
        let i = 0;
        i < experimental.q.length;
        i++
    ) {

        const Iexp =
            experimental.intensity[i];

        const Icalc =
            calculatedIntensity[i];

        const sigma =
            experimental.sigma[i];


        numerator +=
            Iexp *
            Icalc /
            (sigma * sigma);

        denominator +=
            Icalc *
            Icalc /
            (sigma * sigma);
    }


    if (denominator === 0) {
        throw new Error(
            "Cannot determine intensity scale."
        );
    }


    const scaleFactor =
        numerator /
        denominator;


    const scaledIntensity: number[] = [];
    const residuals: number[] = [];


    let chiSquared = 0;


    for (
        let i = 0;
        i < experimental.q.length;
        i++
    ) {

        const scaled =
            scaleFactor *
            calculatedIntensity[i];


        const residual =
            (
                experimental.intensity[i] -
                scaled
            ) /
            experimental.sigma[i];


        scaledIntensity.push(
            scaled
        );

        residuals.push(
            residual
        );

        chiSquared +=
            residual * residual;
    }


    const degreesOfFreedom =
        experimental.q.length - 1;


    const reducedChiSquared =
        chiSquared /
        degreesOfFreedom;


    return {

        scaleFactor,

        calculatedIntensity,

        scaledIntensity,

        residuals,

        chiSquared,

        reducedChiSquared

    };
}