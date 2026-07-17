/**
 * WebGPU implementation of the Debye scattering calculation.
 *
 * Responsibilities:
 * - Upload atomic coordinates
 * - Upload scattering amplitudes
 * - Dispatch compute shader
 * - Reduce partial intensities
 * - Return SAXS intensity curve
 */

export class DebyeGPU {

    device!: GPUDevice;

    pipeline!: GPUComputePipeline;

    positionBuffer!: GPUBuffer;

    amplitudeBuffer!: GPUBuffer;

    outputBuffer!: GPUBuffer;

    paramBuffer!: GPUBuffer;

    bindGroup!: GPUBindGroup;

    atomCount = 0;
    reducePipeline!: GPUComputePipeline;

    reduceBindGroup!: GPUBindGroup;

    reduceBuffer!: GPUBuffer;

    reduceParamBuffer!: GPUBuffer;

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
        const reduceShader =
            await fetch("/reduce.wgsl")
                .then(r => r.text());

        const reduceModule =
            this.device.createShaderModule({
                code: reduceShader
            });

        this.reducePipeline =
            this.device.createComputePipeline({

                layout: "auto",

                compute: {
                    module: reduceModule,
                    entryPoint: "main"
                }

            });
    }

    uploadBuffers(
        positions: Float32Array,
        amplitudes: Float32Array
    )
    {
        this.atomCount = amplitudes.length;
        console.log("positionBuffer", positions.byteLength);
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
        console.log("amplitudeBuffer", amplitudes.byteLength);

        this.amplitudeBuffer  =
            this.device.createBuffer({

                size: amplitudes.byteLength,

                usage:
                    GPUBufferUsage.STORAGE |
                    GPUBufferUsage.COPY_DST

            });

        this.device.queue.writeBuffer(
            this.amplitudeBuffer ,
            0,
            amplitudes
        );

        this.outputBuffer =
            this.device.createBuffer({

                size: this.atomCount * 4,
                usage:
                    GPUBufferUsage.STORAGE |
                    GPUBufferUsage.COPY_SRC |
                    GPUBufferUsage.COPY_DST

            });
        
        this.paramBuffer =
            this.device.createBuffer({

                size:16,

                usage:
                    GPUBufferUsage.UNIFORM |
                    GPUBufferUsage.COPY_DST

            });


        this.bindGroup =
            this.device.createBindGroup({

                layout:
                    this.pipeline.getBindGroupLayout(0),

                entries:[

                {
                    binding:0,
                    resource:{
                        buffer:this.positionBuffer
                    }
                },

                {
                    binding:1,
                    resource:{
                        buffer:this.amplitudeBuffer,

                    }
                },

                {
                    binding:2,
                    resource:{
                        buffer:this.outputBuffer
                    }
                },

                {
                    binding:3,
                    resource:{
                        buffer:this.paramBuffer
                    }
                }

                ]

            });
this.reduceBuffer =
    this.device.createBuffer({

        size: 4,

        usage:
            GPUBufferUsage.STORAGE |
            GPUBufferUsage.COPY_SRC

    });

this.reduceParamBuffer =
    this.device.createBuffer({

        size: 16,

        usage:
            GPUBufferUsage.UNIFORM |
            GPUBufferUsage.COPY_DST

    });

this.reduceBindGroup =
    this.device.createBindGroup({

        layout:
            this.reducePipeline.getBindGroupLayout(0),

        entries: [

            {
                binding: 0,
                resource: {
                    buffer: this.outputBuffer
                }
            },

            {
                binding: 1,
                resource: {
                    buffer: this.reduceBuffer
                }
            },

            {
                binding: 2,
                resource: {
                    buffer: this.reduceParamBuffer
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

    

    async calculate(q:number)
    {
        const buffer = new ArrayBuffer(16);

        const u32 = new Uint32Array(buffer);
        const f32 = new Float32Array(buffer);

        u32[0] = this.atomCount;
        f32[1] = q;

        this.device.queue.writeBuffer(
            this.paramBuffer,
            0,
            buffer
        );
        this.device.queue.writeBuffer(
            this.reduceParamBuffer,
            0,
            buffer
        );

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


        pass.dispatchWorkgroups(
            Math.ceil(this.atomCount / 64)
        );
        


        pass.end();
        const reducePass =
            encoder.beginComputePass();

        reducePass.setPipeline(
            this.reducePipeline
        );

        reducePass.setBindGroup(
            0,
            this.reduceBindGroup
        );

        reducePass.dispatchWorkgroups(1);

        reducePass.end();


        this.device.queue.submit([
            encoder.finish()
        ]);

        const readBuffer =
        this.device.createBuffer({

            size: 4,
            usage:
                GPUBufferUsage.COPY_DST |
                GPUBufferUsage.MAP_READ

        });

    const copyEncoder =
        this.device.createCommandEncoder();

    copyEncoder.copyBufferToBuffer(
        this.reduceBuffer,
        0,
        readBuffer,
        0,
        4
    );

    this.device.queue.submit([
        copyEncoder.finish()
    ]);

    await readBuffer.mapAsync(
        GPUMapMode.READ
    );

    const result =
        new Float32Array(
            readBuffer.getMappedRange()
        );

    const intensity = result[0];

    console.log(
        "GPU intensity:",
        intensity
    );

    readBuffer.unmap();

    return intensity;

    // for (let i = 0; i < result.length; i++)
    // {
    //     intensity += result[i];
    // }

    // console.log(
    //     "GPU intensity:",
    //     intensity
    // );

    // readBuffer.unmap();

    // return intensity;

        // return this.outputBuffer;
    }
    async calculateCurve(
        qMin: number,
        qMax: number,
        qStep: number
    )
    {
        const qValues: number[] = [];
        const intensities: number[] = [];

        for (
            let q = qMin;
            q <= qMax + 1e-8;
            q += qStep
        )
        {
            const intensity =
                await this.calculate(q);

            qValues.push(q);
            intensities.push(intensity);
        }

        return {
            qValues,
            intensities
        };
    }
}