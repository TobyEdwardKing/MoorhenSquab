import { DebyeGPU } from "./vendor/webgpu/DebyeGPU.js";
import { debyeWgsl, reduceWgsl } from "./generated/assets.js";
import type { SAXSCurve } from "./vendor/utils/SAXSScaling.js";
import type { SquabCalculator } from "./SquabCalculator.js";

export interface SquabGPUOptions {
    /** Reuse an existing device instead of requesting a new adapter/device. */
    device?: GPUDevice;
    powerPreference?: GPUPowerPreference;
}

/** True when this browser exposes WebGPU at all (an adapter may still be unavailable). */
export function isWebGPUAvailable(): boolean {
    return typeof navigator !== "undefined" && "gpu" in navigator && !!navigator.gpu;
}

/**
 * Toby's DebyeGPU with the WGSL shaders bundled, instead of fetched from
 * "/debye.wgsl" and "/reduce.wgsl" on the host page.
 */
export class SquabDebyeGPU extends DebyeGPU {
    declare device: GPUDevice;
    declare pipeline: GPUComputePipeline;
    declare reducePipeline: GPUComputePipeline;

    async initialise(options: SquabGPUOptions = {}): Promise<void> {
        let device = options.device;
        if (!device) {
            if (!isWebGPUAvailable()) throw new Error("WebGPU not supported");
            const adapter = await navigator.gpu.requestAdapter({ powerPreference: options.powerPreference });
            if (!adapter) throw new Error("No GPU adapter");
            device = await adapter.requestDevice();
        }
        this.device = device;
        this.pipeline = device.createComputePipeline({
            layout: "auto",
            compute: { module: device.createShaderModule({ code: debyeWgsl }), entryPoint: "main" },
        });
        this.reducePipeline = device.createComputePipeline({
            layout: "auto",
            compute: { module: device.createShaderModule({ code: reduceWgsl }), entryPoint: "main" },
        });
    }
}

/**
 * Atom–atom Debye curve on the GPU for the structure currently prepared in
 * `calculator` (Toby's "GPU Debye reference"; no solvent terms yet).
 */
export async function calculateDebyeCurveGPU(
    calculator: SquabCalculator,
    gpu: SquabDebyeGPU,
    batchSize = 100
): Promise<SAXSCurve> {
    const state = calculator.preparedState;
    if (!state) throw new Error("moorhen-squab: call prepare() first");

    gpu.uploadBuffers(calculator.gpuPositions(), calculator.gpuAmplitudes());

    const q = new Float32Array(state.qValues);
    const qValues: number[] = [];
    const intensities: number[] = [];
    for (let start = 0; start < q.length; start += batchSize) {
        const batch = await gpu.calculateCurve(q.slice(start, start + batchSize), start);
        qValues.push(...batch.qValues);
        intensities.push(...batch.intensities);
    }
    return { qValues, intensities };
}
