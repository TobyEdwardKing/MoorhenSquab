#ifndef GPUCONTEXT_HPP
#define GPUCONTEXT_HPP

#include <webgpu/webgpu_cpp.h>
// TEMP STUB (no WebGPU yet)

// class GPUContext
// {
// public:
//     GPUContext();

//     bool initialise();
//     bool is_valid() const;

// private:
//     bool initialised_ = false;
// };

// #endif


class GPUContext
{
public:

    GPUContext();

    bool initialise();

    bool is_valid() const;

    const wgpu::Device& device() const;
    const wgpu::Queue& queue() const;

private:

    wgpu::Instance instance_;

    wgpu::Adapter adapter_;

    wgpu::Device device_;

    wgpu::Queue queue_;

    bool initialised_ = false;
};

#endif