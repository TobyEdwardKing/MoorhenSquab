import type { SAXSProfile } from "../types/SAXSProfile";

export interface SAXSCurve {
    qValues: number[];
    intensities: number[];
}

export interface SAXSFitResult {
    scale: number;

    chiSquared: number;

    reducedChiSquared: number;

    residuals: number[];

    scaledIntensities: number[];

    qValues: number[];

    difference: number[];
    
    offset?: number;
}


/**
 * Fit a multiplicative scale factor between
 * an experimental SAXS profile and a calculated
 * SAXS curve.
 *
 * The scale factor is chosen to minimise:
 *
 *   chi2 = sum[
 *       (I_exp - scale * I_calc)^2 / sigma^2
 *   ]
 *
 * The analytical solution is:
 *
 *   scale =
 *       sum(I_exp * I_calc / sigma^2)
 *       --------------------------------
 *       sum(I_calc^2 / sigma^2)
 *
 * No normalisation of either curve is performed.
 */
export function fitSAXSScale(
    experimental: SAXSProfile,
    calculated: SAXSCurve
): SAXSFitResult {

    if (
        experimental.q.length !==
        calculated.qValues.length
    ) {
        throw new Error(
            "Experimental and calculated curves have different numbers of points."
        );
    }

    const n =
        experimental.q.length;

    if (n === 0) {
        throw new Error(
            "Cannot fit scale to an empty curve."
        );
    }

    let numerator = 0.0;
    let denominator = 0.0;

    /*
     * Calculate the optimal multiplicative scale.
     *
     * Every point is weighted by 1 / sigma^2.
     */
    for (let i = 0; i < n; i++) {

        const Iexp =
            experimental.intensity[i];

        const Icalc =
            calculated.intensities[i];

        const sigma =
            experimental.sigma[i];

        if (!Number.isFinite(Iexp)) {
            continue;
        }

        if (!Number.isFinite(Icalc)) {
            continue;
        }

        if (
            !Number.isFinite(sigma) ||
            sigma <= 0
        ) {
            continue;
        }

        const weight =
            1.0 / (sigma * sigma);

        numerator +=
            Iexp *
            Icalc *
            weight;

        denominator +=
            Icalc *
            Icalc *
            weight;
    }

    if (denominator === 0) {
        throw new Error(
            "Cannot determine scale: calculated intensities have zero weighted magnitude."
        );
    }

    const scale =
        numerator / denominator;

    /*
     * Calculate scaled curve, differences and residuals.
     */
    const scaledIntensities: number[] = [];
    const residuals: number[] = [];
    const difference = new Array<number>(n);

    let chiSquaredSum = 0.0;
    let validPoints = 0;

    for (let i = 0; i < n; i++) {

        const Iexp =
            experimental.intensity[i];

        const Icalc =
            calculated.intensities[i];

        const sigma =
            experimental.sigma[i];

        /*
         * Check that all values needed for the fit
         * are valid before calculating anything.
         */
        if (
            !Number.isFinite(Iexp) ||
            !Number.isFinite(Icalc) ||
            !Number.isFinite(sigma) ||
            sigma <= 0
        ) {

            scaledIntensities.push(
                NaN
            );

            difference[i] =
                NaN;

            residuals.push(
                NaN
            );

            continue;
        }

        /*
         * Apply the fitted scale factor.
         */
        const scaled =
            scale * Icalc;

        scaledIntensities.push(
            scaled
        );

        /*
         * Raw intensity difference:
         *
         *   I_exp - I_calc_scaled
         */
        difference[i] =
            Iexp - scaled;

        /*
         * Normalised residual:
         *
         *   (I_exp - I_calc_scaled) / sigma
         */
        const residual =
            (Iexp - scaled) / sigma;

        residuals.push(
            residual
        );

        chiSquaredSum +=
            residual * residual;

        validPoints++;
    }

    /*
     * This is the conventional reduced chi-squared:
     *
     *   chi2 / degrees_of_freedom
     *
     * We fitted one parameter (the scale), therefore:
     *
     *   DOF = N - 1
     */
    const chiSquared =
        chiSquaredSum;

    const degreesOfFreedom =
        validPoints - 1;

    const reducedChiSquared =
        degreesOfFreedom > 0
            ? chiSquared / degreesOfFreedom
            : NaN;

    /*
     * Print the first 10 points for debugging.
     */
    // console.table(
    //     experimental.q.slice(0, 10).map((q, i) => ({
    //         q,
    //         Iexp: experimental.intensity[i],
    //         scaledIcalc: scaledIntensities[i],
    //         difference: difference[i],
    //     }))
    // );

    return {

        scale,

        chiSquared,

        reducedChiSquared,

        residuals,

        scaledIntensities,

        qValues:
            [...experimental.q],

        difference
    };
}


export function fitSAXSScaleWithOffset(
    experimental: SAXSProfile,
    calculated: SAXSCurve
): SAXSFitResult {

    if (
        experimental.q.length !== calculated.qValues.length
    ) {
        throw new Error(
            "Experimental and calculated curves have different numbers of points."
        );
    }

    const n = experimental.q.length;

    if (n < 2) {
        throw new Error(
            "At least two points are required to fit scale and offset."
        );
    }

    let sumTheorOverSigma2 = 0;
    let sumTheorSquaredOverSigma2 = 0;

    let validPoints = 0;

    // Calculate b from SI equation (5).
    for (let i = 0; i < n; i++) {
        const Iexp = experimental.intensity[i];
        const Itheor = calculated.intensities[i];
        const sigma = experimental.sigma[i];

        if (
            !Number.isFinite(Iexp) ||
            !Number.isFinite(Itheor) ||
            !Number.isFinite(sigma) ||
            sigma <= 0
        ) {
            continue;
        }

        const weight = 1 / (sigma * sigma);

        sumTheorOverSigma2 += Itheor * weight;
        sumTheorSquaredOverSigma2 +=
            Itheor * Itheor * weight;

        validPoints++;
    }

    if (
        validPoints < 2 ||
        sumTheorSquaredOverSigma2 === 0
    ) {
        throw new Error(
            "Cannot determine scale and offset from the available data."
        );
    }

    const b =
        sumTheorOverSigma2 /
        sumTheorSquaredOverSigma2;

    // Calculate the ordinary scale factor c.
    const scaleResult = fitSAXSScale(
        experimental,
        calculated
    );

    const c = scaleResult.scale;

    // Calculate the optimal offset kappa using the SI expression.
    let offsetNumerator = 0;
    let offsetDenominator = 0;

    for (let i = 0; i < n; i++) {
        const Iexp = experimental.intensity[i];
        const Itheor = calculated.intensities[i];
        const sigma = experimental.sigma[i];

        if (
            !Number.isFinite(Iexp) ||
            !Number.isFinite(Itheor) ||
            !Number.isFinite(sigma) ||
            sigma <= 0
        ) {
            continue;
        }

        const weight = 1 / (sigma * sigma);

        const residualWithoutOffset =
            Iexp - c * Itheor;

        const offsetTerm =
            1 - b * Itheor;

        offsetNumerator +=
            residualWithoutOffset *
            offsetTerm *
            weight;

        offsetDenominator +=
            offsetTerm *
            offsetTerm *
            weight;
    }

    if (
        !Number.isFinite(offsetDenominator) ||
        offsetDenominator <= 0
    ) {
        throw new Error(
            "Cannot determine the optimal background offset."
        );
    }

    const offset =
        offsetNumerator / offsetDenominator;

    // Calculate the offset-adjusted scale c(kappa).
    const adjustedScale =
        c + offset * b;

    // Calculate the fitted curve and goodness of fit.
    const scaledIntensities: number[] = [];
    const residuals: number[] = [];
    const difference = new Array<number>(n);

    let chiSquared = 0;
    let fittedPoints = 0;

    for (let i = 0; i < n; i++) {
        const Iexp = experimental.intensity[i];
        const Itheor = calculated.intensities[i];
        const sigma = experimental.sigma[i];

        if (
            !Number.isFinite(Iexp) ||
            !Number.isFinite(Itheor) ||
            !Number.isFinite(sigma) ||
            sigma <= 0
        ) {
            scaledIntensities.push(NaN);
            residuals.push(NaN);
            difference[i] = NaN;
            continue;
        }

        const fitted =
            adjustedScale * Itheor - offset;

        // SI residual: (Iexp + kappa - c(kappa)*Itheor) / sigma.
        const residual =
            (Iexp + offset - adjustedScale * Itheor) /
            sigma;

        scaledIntensities.push(fitted);
        difference[i] = Iexp - fitted;
        residuals.push(residual);

        chiSquared += residual * residual;
        fittedPoints++;
    }

    // Match the SI's stated normalisation: N - 1.
    const reducedChiSquared =
        fittedPoints > 1
            ? chiSquared / (fittedPoints - 1)
            : NaN;

    return {
        scale: adjustedScale,
        offset,
        chiSquared,
        reducedChiSquared,
        residuals,
        scaledIntensities,
        qValues: [...experimental.q],
        difference
    };
}


