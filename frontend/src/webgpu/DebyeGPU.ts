/**
 * WebGPU implementation of the Debye scattering calculation.
 *
 * Responsibilities:
 * - Upload atomic coordinates
 * - Upload scattering amplitudes
 * - Dispatch compute shader
 * - Calculate Debye intensities for batches of q values
 * - Return the SAXS intensity curve
 */

export class DebyeGPU {

    device!: GPUDevice;

    pipeline!: GPUComputePipeline;

    positionBuffer!: GPUBuffer;

    amplitudeBuffer!: GPUBuffer;

    qBuffer!: GPUBuffer;

    outputBuffer!: GPUBuffer;

    paramBuffer!: GPUBuffer;

    bindGroup!: GPUBindGroup;

    atomCount = 0;
    qCount = 0;

    async initialise() {


        if (!navigator.gpu) {
            throw new Error(
                "WebGPU not supported"
            );
        }


        const adapter =
            await navigator.gpu.requestAdapter();


        if (!adapter) {
            throw new Error(
                "No GPU adapter"
            );
        }


        this.device =
            await adapter.requestDevice();



        const shader =
            await fetch(
                "/debye.wgsl"
            )
            .then(r => r.text());



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
        
    }

    uploadBuffers(
        positions: Float32Array,
        amplitudes: Float32Array
    )
    {
        this.atomCount =
            positions.length / 3;

        console.log(
            "positionBuffer",
            positions.byteLength
        );

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
            "amplitudeBuffer",
            amplitudes.byteLength
        );

        this.amplitudeBuffer =
            this.device.createBuffer({

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
    }

    private prepareQBuffers(
        qValues: Float32Array,
        qOffset: number
    )
    {
        this.qCount = qValues.length;


        this.qBuffer =
            this.device.createBuffer({

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


        this.outputBuffer =
            this.device.createBuffer({

                size:
                    qValues.length * 4,

                usage:
                    GPUBufferUsage.STORAGE |
                    GPUBufferUsage.COPY_SRC

            });


        /*
        * Params:
        *
        * uint atomCount
        * uint qCount
        *
        * 16 bytes because of uniform-buffer alignment.
        */

        this.paramBuffer =
            this.device.createBuffer({

                size: 16,

                usage:
                    GPUBufferUsage.UNIFORM |
                    GPUBufferUsage.COPY_DST

            });


        const paramData =
            new ArrayBuffer(16);

        const params =
            new Uint32Array(paramData);

        params[0] = this.atomCount;
        params[1] = this.qCount;
        params[2] = qOffset;
        
        this.device.queue.writeBuffer(
            this.paramBuffer,
            0,
            paramData
        );


        /*
        * Bind everything required by debye.wgsl.
        */

        this.bindGroup =
            this.device.createBindGroup({

                layout:
                    this.pipeline
                        .getBindGroupLayout(0),

                entries: [

                    {
                        binding: 0,
                        resource: {
                            buffer:
                                this.positionBuffer
                        }
                    },

                    {
                        binding: 1,
                        resource: {
                            buffer:
                                this.amplitudeBuffer
                        }
                    },

                    {
                        binding: 2,
                        resource: {
                            buffer:
                                this.outputBuffer
                        }
                    },

                    {
                        binding: 3,
                        resource: {
                            buffer:
                                this.qBuffer
                        }
                    },

                    {
                        binding: 4,
                        resource: {
                            buffer:
                                this.paramBuffer
                        }
                    }

                ]

            });
    }
    async readOutput()
    {
        const size = 4;


        const staging =
            this.device.createBuffer({

                size,

                usage:
                    GPUBufferUsage.COPY_DST |
                    GPUBufferUsage.MAP_READ

            });


        const encoder =
            this.device.createCommandEncoder();


        encoder.copyBufferToBuffer(
            this.outputBuffer,
            0,
            staging,
            0,
            size
        );


        this.device.queue.submit([
            encoder.finish()
        ]);


        await staging.mapAsync(
            GPUMapMode.READ
        );


        const data =
            new Float32Array(
                staging
                .getMappedRange()
            );


        return data;
    }

    

async calculateCurve(
    qValues: Float32Array,
    qOffset = 0
)
{
    if (this.atomCount === 0) {
        throw new Error(
            "No atom buffers have been uploaded."
        );
    }

    if (qValues.length === 0) {
        throw new Error(
            "No q values supplied."
        );
    }


    const curveStart =
        performance.now();


    this.prepareQBuffers(
        qValues,
        qOffset
    );

    /*
     * Readback buffer.
     */

    const readBuffer =
        this.device.createBuffer({

            size:
                qValues.byteLength,

            usage:
                GPUBufferUsage.COPY_DST |
                GPUBufferUsage.MAP_READ

        });


    /*
     * Encode GPU work.
     */

    const encoder =
        this.device.createCommandEncoder();


    const pass =
        encoder.beginComputePass();


    pass.setPipeline(
        this.pipeline
    );

    pass.setBindGroup(
        0,
        this.bindGroup
    );


    /*
     * One workgroup per q value.
     */

    pass.dispatchWorkgroups(
        qValues.length
    );


    pass.end();


    /*
     * Copy all calculated intensities
     * back to the CPU in one operation.
     */

    encoder.copyBufferToBuffer(
        this.outputBuffer,
        0,
        readBuffer,
        0,
        qValues.byteLength
    );


    this.device.queue.submit([
        encoder.finish()
    ]);


    /*
     * Wait for GPU completion.
     */

    await readBuffer.mapAsync(
        GPUMapMode.READ
    );


    const intensities =
        new Float32Array(
            readBuffer.getMappedRange()
        );


    /*
     * Copy because the mapped buffer
     * becomes invalid after unmap().
     */

    const intensityArray =
        new Float32Array(
            intensities
        );


    readBuffer.unmap();


    const curveEnd =
        performance.now();


    console.log(
        `[GPU] COMPLETE: ` +
        `${qValues.length} q values in ` +
        `${(curveEnd - curveStart).toFixed(2)} ms`
    );


    return {

        qValues:
            Array.from(qValues),

        intensities:
            Array.from(intensityArray)

    };
}
}