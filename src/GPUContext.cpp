#include "GPUContext.hpp"

GPUContext::GPUContext()
{
}

bool GPUContext::initialise()
{
    initialised_ = true;
    return true;
}

bool GPUContext::is_valid() const
{
    return initialised_;
}
// #include "GPUContext.hpp"

// #include <iostream>

// GPUContext::GPUContext()
// {
//     instance_ = wgpu::CreateInstance();
// }

// bool GPUContext::initialise()
// {
//     if (!instance_) {
//         std::cerr << "Failed to create WebGPU instance\n";
//         return false;
//     }

//     // 1. Request adapter (prefer high performance GPU if available)
//     wgpu::RequestAdapterOptions options = {};
//     options.powerPreference = wgpu::PowerPreference::HighPerformance;

//     bool adapterReady = false;

//     instance_.RequestAdapter(
//         &options,
//         wgpu::CallbackMode::WaitAnyOnly,
//         [&adapterReady, this](wgpu::RequestAdapterStatus status,
//                               wgpu::Adapter adapter,
//                               const char* message)
//         {
//             if (status != wgpu::RequestAdapterStatus::Success) {
//                 std::cerr << "RequestAdapter failed: "
//                           << (message ? message : "unknown error") << "\n";
//                 adapterReady = false;
//                 return;
//             }

//             adapter_ = std::move(adapter);
//             adapterReady = true;
//         }
//     );

//     if (!adapterReady || !adapter_) {
//         std::cerr << "No suitable WebGPU adapter found\n";
//         return false;
//     }

//     // 2. Create device
//     wgpu::DeviceDescriptor deviceDesc = {};

//     // Optional: add debug/error logging
//     deviceDesc.SetDefault();

//     bool deviceReady = false;

//     adapter_.RequestDevice(
//         &deviceDesc,
//         wgpu::CallbackMode::WaitAnyOnly,
//         [&deviceReady, this](wgpu::RequestDeviceStatus status,
//                              wgpu::Device device,
//                              const char* message)
//         {
//             if (status != wgpu::RequestDeviceStatus::Success) {
//                 std::cerr << "RequestDevice failed: "
//                           << (message ? message : "unknown error") << "\n";
//                 deviceReady = false;
//                 return;
//             }

//             device_ = std::move(device);
//             deviceReady = true;
//         }
//     );

//     if (!deviceReady || !device_) {
//         std::cerr << "Failed to create WebGPU device\n";
//         return false;
//     }

//     // 3. Get queue
//     queue_ = device_.GetQueue();

//     if (!queue_) {
//         std::cerr << "Failed to get WebGPU queue\n";
//         return false;
//     }

//     initialised_ = true;

//     std::cout << "WebGPU initialised successfully\n";
//     return true;
// }

// bool GPUContext::is_valid() const
// {
//     return initialised_ && device_ && queue_;
// }

// const wgpu::Device& GPUContext::device() const
// {
//     return device_;
// }

// const wgpu::Queue& GPUContext::queue() const
// {
//     return queue_;
// }