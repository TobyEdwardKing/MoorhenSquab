#ifndef DEBYEPIPELINE_HPP
#define DEBYEPIPELINE_HPP

#include <vector>
#include <webgpu/webgpu_cpp.h>

#include "GPUContext.hpp"

struct AtomData
{
    float x;
    float y;
    float z;
    uint32_t Z;
};

class DebyePipeline
{
public:
    explicit DebyePipeline(GPUContext& context);

    bool initialise();

    bool set_data(const std::vector<AtomData>& atoms,
                  const std::vector<float>& formFactors,
                  const std::vector<float>& qValues);

    bool compute();

    const std::vector<float>& result() const;

private:
    GPUContext& context_;

    // Pipeline
    wgpu::ShaderModule shader_;
    wgpu::ComputePipeline pipeline_;
    wgpu::BindGroup bind_group_;

    // Buffers
    wgpu::Buffer atom_buffer_;
    wgpu::Buffer form_factor_buffer_;
    wgpu::Buffer q_buffer_;
    wgpu::Buffer intensity_buffer_;
    wgpu::Buffer staging_buffer_;
    wgpu::Buffer uniform_buffer_;

    // CPU copy
    std::vector<float> host_result_;

    uint32_t atom_count_ = 0;
    uint32_t q_count_ = 0;

    void create_buffers();
    void create_pipeline();
};

#endif

// #ifndef DEBYEPIPELINE_HPP
// #define DEBYEPIPELINE_HPP

// #include <vector>

// #include "atomExtractor.hpp"
// #include "formFactor.hpp"

// class GPUContext; // forward declare

// class DebyePipeline
// {
// public:
//     DebyePipeline(GPUContext& gpu);

//     void initialise();

//     void set_data(
//         const std::vector<AtomData>& atoms,
//         const FormFactorTable& table,
//         const std::vector<float>& qValues);

//     void compute();

//     std::vector<float> result() const;

// private:
//     GPUContext& gpu_;

//     std::vector<float> output_;
// };

// #endif