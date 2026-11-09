/**
 * WebGPU implementation of the Debye scattering calculation.
 *
 * Responsibilities:
 * - Upload atomic coordinates
 * - Keep scattering amplitudes in CPU memory
 * - Upload amplitudes in q batches to avoid GPU storage-buffer limits
 * - Dispatch compute shader
 * - Calculate Debye intensities
 * - Return the SAXS intensity curve
 */

export class DebyeGPU {

    device!: GPUDevice;

    pipeline!: GPUComputePipeline;
    reductionPipeline!: GPUComputePipeline;

    positionBuffer!: GPUBuffer;

    amplitudeBuffer!: GPUBuffer;

    qBuffer!: GPUBuffer;

    outputBuffer!: GPUBuffer;

    paramBuffer!: GPUBuffer;

    bindGroup!: GPUBindGroup;

    atomCount = 0;
    qCount = 0;
    finalOutputBuffer!: GPUBuffer;
    reductionBindGroup!: GPUBindGroup;
    private amplitudes!: Float32Array;


    async initialise() {

        if (!navigator.gpu) {
            throw new Error("WebGPU not supported");
        }

        const adapter =
            await navigator.gpu.requestAdapter();

        if (!adapter) {
            throw new Error("No GPU adapter");
        }

        this.device =
            await adapter.requestDevice();

        const shader =
            await fetch("/debye.wgsl")
                .then(r => {
                    if (!r.ok) {
                        throw new Error("Could not load debye.wgsl");
                    }
                    return r.text();
                });

        const module =
            this.device.createShaderModule({
                code: shader
            });

        this.pipeline =
            this.device.createComputePipeline({
                layout: "auto",
                compute: {
                    module,
                    entryPoint: "main"
                }
            });
        const reductionShader = await fetch("/debye-reduce.wgsl")
            .then(r => {
                if (!r.ok) {
                    throw new Error("Could not load debye-reduce.wgsl");
                }
                return r.text();
            });

        const reductionModule = this.device.createShaderModule({
            code: reductionShader
        });

        this.reductionPipeline = this.device.createComputePipeline({
            layout: "auto",
            compute: {
                module: reductionModule,
                entryPoint: "main"
            }
        });
    }


    uploadBuffers(
        positions: Float32Array,
        amplitudes: Float32Array
    ) {
        if (positions.length === 0 || positions.length % 3 !== 0) {
            throw new Error("Invalid position array.");
        }

        this.atomCount = positions.length / 3;

        if (amplitudes.length % this.atomCount !== 0) {
            throw new Error(
                "Amplitude array size is inconsistent with atom count."
            );
        }

        if (
            positions.byteLength >
            this.device.limits.maxStorageBufferBindingSize
        ) {
            throw new Error(
                "Position buffer exceeds the GPU storage-buffer binding limit."
            );
        }

        // Keep amplitudes in CPU memory.
        // They will be uploaded in smaller q batches during calculation.
        this.amplitudes = amplitudes;

        // Upload atomic positions once.
        this.positionBuffer =
            this.device.createBuffer({
                size: positions.byteLength,
                usage:
                    GPUBufferUsage.STORAGE |
                    GPUBufferUsage.COPY_DST
            });

        this.device.queue.writeBuffer(
            this.positionBuffer,
            0,
            positions
        );

        console.log(
            "positionBuffer",
            positions.byteLength
        );

        console.log(
            "amplitudeBuffer (CPU)",
            amplitudes.byteLength
        );
    }


    /**
     * Determine how many q values can fit into one amplitude buffer.
     *
     * A 64 MiB cap is used even if the GPU supports larger buffers.
     * This leaves room for other GPU resources.
     */
    private getMaxQBatchSize(): number {

        const bytesPerQ =
            this.atomCount * Float32Array.BYTES_PER_ELEMENT;

        const maxBufferBytes = Math.min(
            this.device.limits.maxStorageBufferBindingSize,
            64 * 1024 * 1024
        );

        const maxQ =
            Math.floor(maxBufferBytes / bytesPerQ);

        if (maxQ < 1) {
            throw new Error(
                "One q-value amplitude array exceeds the GPU storage-buffer limit."
            );
        }

        return maxQ;
    }


    /**
     * Prepare GPU buffers for one batch of q values.
     */
    private prepareQBuffers(
        qValues: Float32Array,
        amplitudes: Float32Array,
        qOffset: number
    ) {
        this.qCount = qValues.length;

        const iTileCount = Math.ceil(this.atomCount / 64);

        // q values.
        this.qBuffer = this.device.createBuffer({
            size: qValues.byteLength,
            usage:
                GPUBufferUsage.STORAGE |
                GPUBufferUsage.COPY_DST
        });

        this.device.queue.writeBuffer(
            this.qBuffer,
            0,
            qValues
        );

        // Amplitudes for this q batch only.
        this.amplitudeBuffer = this.device.createBuffer({
            size: amplitudes.byteLength,
            usage:
                GPUBufferUsage.STORAGE |
                GPUBufferUsage.COPY_DST
        });

        this.device.queue.writeBuffer(
            this.amplitudeBuffer,
            0,
            amplitudes
        );

        // One partial intensity for every q / iTile combination.
        const partialOutputBytes =
            qValues.length * iTileCount * Float32Array.BYTES_PER_ELEMENT;
        if (
            partialOutputBytes >
            this.device.limits.maxStorageBufferBindingSize
        ) {
            throw new Error(
                "Partial output buffer exceeds the GPU storage-buffer binding limit."
            );
        }
        this.outputBuffer = this.device.createBuffer({
            size: partialOutputBytes,
            usage:
                GPUBufferUsage.STORAGE |
                GPUBufferUsage.COPY_SRC
        });

        // One final intensity per q value.
        this.finalOutputBuffer = this.device.createBuffer({
            size: qValues.length * Float32Array.BYTES_PER_ELEMENT,
            usage:
                GPUBufferUsage.STORAGE |
                GPUBufferUsage.COPY_SRC
        });

        // Shader parameters.
        const paramData = new ArrayBuffer(16);
        const params = new Uint32Array(paramData);

        params[0] = this.atomCount;
        params[1] = this.qCount;
        params[2] = qOffset;
        params[3] = iTileCount;

        this.paramBuffer = this.device.createBuffer({
            size: 16,
            usage:
                GPUBufferUsage.UNIFORM |
                GPUBufferUsage.COPY_DST
        });

        this.device.queue.writeBuffer(
            this.paramBuffer,
            0,
            paramData
        );

        // Main Debye shader bind group.
        this.bindGroup = this.device.createBindGroup({
            layout: this.pipeline.getBindGroupLayout(0),
            entries: [
                {
                    binding: 0,
                    resource: { buffer: this.positionBuffer }
                },
                {
                    binding: 1,
                    resource: { buffer: this.amplitudeBuffer }
                },
                {
                    binding: 2,
                    resource: { buffer: this.outputBuffer }
                },
                {
                    binding: 3,
                    resource: { buffer: this.qBuffer }
                },
                {
                    binding: 4,
                    resource: { buffer: this.paramBuffer }
                }
            ]
        });

        // Reduction shader bind group.
        this.reductionBindGroup = this.device.createBindGroup({
            layout: this.reductionPipeline.getBindGroupLayout(0),
            entries: [
                {
                    binding: 0,
                    resource: { buffer: this.outputBuffer }
                },
                {
                    binding: 1,
                    resource: { buffer: this.finalOutputBuffer }
                },
                {
                    binding: 2,
                    resource: { buffer: this.paramBuffer }
                }
            ]
        });
    }


    /**
     * Calculate Debye intensities.
     *
     * The q grid is automatically divided into batches so that
     * the amplitude buffer does not exceed the GPU storage limit.
     */
    async calculateCurve(
        qValues: Float32Array,
        qOffset = 0
    ) {
        if (this.atomCount === 0) {
            throw new Error("No atom buffers have been uploaded.");
        }

        if (qValues.length === 0) {
            throw new Error("No q values supplied.");
        }

        const totalQCount =
            this.amplitudes.length / this.atomCount;

        if (
            qOffset < 0 ||
            qOffset + qValues.length > totalQCount
        ) {
            throw new Error(
                "Requested q range exceeds the uploaded amplitude data."
            );
        }

        // Temporary conservative batch limit while testing.
        const maxBatchSize = Math.min(
            this.getMaxQBatchSize(),
            this.atomCount > 5000 ? 1 : 100
        );

        const iTileCount = Math.ceil(this.atomCount / 64);

        const allIntensities: number[] = [];

        const curveStart = performance.now();

        console.log(`[GPU] Using q batch size: ${maxBatchSize}`);
        console.log(`[GPU] Atom blocks per q: ${iTileCount}`);

        for (
            let start = 0;
            start < qValues.length;
            start += maxBatchSize
        ) {
            const end = Math.min(
                start + maxBatchSize,
                qValues.length
            );

            const qBatch = qValues.slice(start, end);

            const amplitudeStart =
                (qOffset + start) * this.atomCount;

            const amplitudeEnd =
                (qOffset + end) * this.atomCount;

            const amplitudeBatch = this.amplitudes.subarray(
                amplitudeStart,
                amplitudeEnd
            );

            this.prepareQBuffers(
                qBatch,
                amplitudeBatch,
                qOffset + start
            );

            const readBuffer = this.device.createBuffer({
                size: qBatch.byteLength,
                usage:
                    GPUBufferUsage.COPY_DST |
                    GPUBufferUsage.MAP_READ
            });

            try {
                const encoder = this.device.createCommandEncoder();

                // ----------------------------------------------------
                // First pass: calculate partial Debye intensities.
                // ----------------------------------------------------

                const pass = encoder.beginComputePass();

                pass.setPipeline(this.pipeline);
                pass.setBindGroup(0, this.bindGroup);

                pass.dispatchWorkgroups(
                    qBatch.length,
                    iTileCount
                );

                pass.end();

                // ----------------------------------------------------
                // Second pass: combine the partial intensities.
                // ----------------------------------------------------

                const reductionPass = encoder.beginComputePass();

                reductionPass.setPipeline(this.reductionPipeline);
                reductionPass.setBindGroup(0, this.reductionBindGroup);

                reductionPass.dispatchWorkgroups(qBatch.length);

                reductionPass.end();

                // Copy final intensities back to CPU.
                encoder.copyBufferToBuffer(
                    this.finalOutputBuffer,
                    0,
                    readBuffer,
                    0,
                    qBatch.byteLength
                );

                this.device.queue.submit([
                    encoder.finish()
                ]);

                await readBuffer.mapAsync(GPUMapMode.READ);

                const mapped = readBuffer.getMappedRange();

                const batchIntensities = new Float32Array(mapped);

                for (let i = 0; i < batchIntensities.length; i++) {
                    allIntensities.push(batchIntensities[i]);
                }

                readBuffer.unmap();

            } finally {
                readBuffer.destroy();

                this.qBuffer.destroy();
                this.amplitudeBuffer.destroy();
                this.outputBuffer.destroy();
                this.finalOutputBuffer.destroy();
                this.paramBuffer.destroy();
            }
        }

        const curveEnd = performance.now();

        console.log(
            `[GPU] COMPLETE: ${qValues.length} q values in ` +
            `${(curveEnd - curveStart).toFixed(2)} ms`
        );

        return {
            qValues: Array.from(qValues),
            intensities: allIntensities
        };
    }
}