#include "DebyePipeline.hpp"

#include <fstream>
#include <sstream>
#include <cstring>
#include <iostream>

// ---------------- WGSL loader ----------------
static std::string loadWGSL(const std::string& path)
{
    std::ifstream file(path);
    if (!file.is_open()) {
        std::cerr << "Failed to open WGSL file: " << path << std::endl;
        return "";
    }

    std::stringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

// ---------------- buffer helper ----------------
static wgpu::Buffer makeBuffer(const wgpu::Device& device,
                               size_t size,
                               wgpu::BufferUsage usage)
{
    wgpu::BufferDescriptor desc = {};
    desc.size = size;
    desc.usage = usage;
    desc.mappedAtCreation = false;

    return device.CreateBuffer(&desc);
}

// ---------------- constructor ----------------
DebyePipeline::DebyePipeline(GPUContext& context)
    : context_(context)
{}

// ---------------- initialise ----------------
bool DebyePipeline::initialise()
{
    if (!context_.is_valid())
        return false;

    std::string code = loadWGSL("shaders/Debye.wgsl");
    if (code.empty())
        return false;

    wgpu::ShaderModuleWGSLDescriptor wgsl = {};
    wgsl.code = code.c_str();

    wgpu::ShaderModuleDescriptor desc = {};
    desc.nextInChain = &wgsl;

    shader_ = context_.device().CreateShaderModule(&desc);

    // ---------------- bind group layout ----------------
    wgpu::BindGroupLayoutEntry bgl[5] = {};

    // atoms
    bgl[0].binding = 0;
    bgl[0].visibility = wgpu::ShaderStage::Compute;
    bgl[0].buffer.type = wgpu::BufferBindingType::ReadOnlyStorage;

    // form factors
    bgl[1] = bgl[0];
    bgl[1].binding = 1;

    // q values
    bgl[2] = bgl[0];
    bgl[2].binding = 2;

    // intensity output
    bgl[3] = bgl[0];
    bgl[3].binding = 3;
    bgl[3].buffer.type = wgpu::BufferBindingType::Storage;

    // uniforms
    bgl[4] = bgl[0];
    bgl[4].binding = 4;
    bgl[4].buffer.type = wgpu::BufferBindingType::Uniform;

    wgpu::BindGroupLayoutDescriptor bglDesc = {};
    bglDesc.entryCount = 5;
    bglDesc.entries = bgl;

    auto layout = context_.device().CreateBindGroupLayout(&bglDesc);

    wgpu::PipelineLayoutDescriptor plDesc = {};
    plDesc.bindGroupLayoutCount = 1;
    plDesc.bindGroupLayouts = &layout;

    wgpu::ComputePipelineDescriptor cpDesc = {};
    cpDesc.layout = context_.device().CreatePipelineLayout(&plDesc);
    cpDesc.compute.module = shader_;
    cpDesc.compute.entryPoint = "main";

    pipeline_ = context_.device().CreateComputePipeline(&cpDesc);

    return true;
}

// ---------------- set data ----------------
bool DebyePipeline::set_data(const std::vector<AtomData>& atoms,
                             const std::vector<float>& formFactors,
                             const std::vector<float>& qValues)
{
    auto& device = context_.device();
    auto& queue  = context_.queue();

    atom_count_ = atoms.size();
    q_count_ = qValues.size();
    FormFactorTable fft;
    fft.load("elementsInfo.json"); // or stored filename

    std::vector<float> gpuFormFactors =
        fft.build_gpu_table(qValues);
        
    // atoms
    atom_buffer_ = makeBuffer(
        device,
        atoms.size() * sizeof(AtomData),
        wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst
    );

    queue.WriteBuffer(atom_buffer_, 0,
        atoms.data(),
        atoms.size() * sizeof(AtomData)
    );

    // form factors

    form_factor_buffer_ = makeBuffer(
        device,
        gpuFormFactors.size() * sizeof(float),
        wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst
    );

    queue.WriteBuffer(
        form_factor_buffer_,
        0,
        gpuFormFactors.data(),
        gpuFormFactors.size() * sizeof(float)
    );

    // q values
    q_buffer_ = makeBuffer(
        device,
        qValues.size() * sizeof(float),
        wgpu::BufferUsage::Storage | wgpu::BufferUsage::CopyDst
    );

    queue.WriteBuffer(q_buffer_, 0,
        qValues.data(),
        qValues.size() * sizeof(float)
    );

    // intensity output
    intensity_buffer_ = makeBuffer(
        device,
        qValues.size() * sizeof(float),
        wgpu::BufferUsage::Storage |
        wgpu::BufferUsage::CopySrc |
        wgpu::BufferUsage::CopyDst
    );

    // staging buffer
    staging_buffer_ = makeBuffer(
        device,
        qValues.size() * sizeof(float),
        wgpu::BufferUsage::MapRead | wgpu::BufferUsage::CopyDst
    );

    // uniforms
    struct Uniforms
    {
        uint32_t atomCount;
        uint32_t qCount;
        uint32_t pad0;
        uint32_t pad1;
    };

    Uniforms u = {
        atom_count_,
        q_count_,
        0,
        0
    };

    uniform_buffer_ = makeBuffer(
        device,
        sizeof(Uniforms),
        wgpu::BufferUsage::Uniform | wgpu::BufferUsage::CopyDst
    );

    queue.WriteBuffer(uniform_buffer_, 0, &u, sizeof(u));

    return true;
}

// ---------------- compute ----------------
bool DebyePipeline::compute()
{
    auto& device = context_.device();
    auto& queue  = context_.queue();

    // bind group
    wgpu::BindGroupEntry e[5] = {};

    e[0].binding = 0;
    e[0].buffer = atom_buffer_;

    e[1].binding = 1;
    e[1].buffer = form_factor_buffer_;

    e[2].binding = 2;
    e[2].buffer = q_buffer_;

    e[3].binding = 3;
    e[3].buffer = intensity_buffer_;

    e[4].binding = 4;
    e[4].buffer = uniform_buffer_;

    wgpu::BindGroupDescriptor bgDesc = {};
    bgDesc.layout = pipeline_.GetBindGroupLayout(0);
    bgDesc.entryCount = 5;
    bgDesc.entries = e;

    bind_group_ = device.CreateBindGroup(&bgDesc);

    // command encoding
    wgpu::CommandEncoder encoder = device.CreateCommandEncoder();

    {
        wgpu::ComputePassEncoder pass = encoder.BeginComputePass();
        pass.SetPipeline(pipeline_);
        pass.SetBindGroup(0, bind_group_);

        uint32_t workgroups = (q_count_ + 63) / 64;
        pass.DispatchWorkgroups(workgroups);

        pass.End();
    }

    encoder.CopyBufferToBuffer(
        intensity_buffer_, 0,
        staging_buffer_, 0,
        q_count_ * sizeof(float)
    );

    wgpu::CommandBuffer cmd = encoder.Finish();
    queue.Submit(1, &cmd);

    return true;
}

// ---------------- readback ----------------
const std::vector<float>& DebyePipeline::result() const
{
    host_result_.resize(q_count_);

    staging_buffer_.MapAsync(
        wgpu::MapMode::Read,
        0,
        q_count_ * sizeof(float),
        [](wgpu::MapAsyncStatus status, void*) {},
        nullptr
    );

    std::memcpy(
        host_result_.data(),
        staging_buffer_.GetMappedRange(),
        q_count_ * sizeof(float)
    );

    staging_buffer_.Unmap();

    return host_result_;
}