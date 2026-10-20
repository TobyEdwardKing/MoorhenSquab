#pragma once

#include <vector>

#include "atomExtractor.hpp"
#include "formFactor.hpp"
#include "scatteringModel.hpp"

class DebyeCalculator
{
public:

    // Generate q-dependent scattering amplitudes for GPU use.
    //
    // positions:
    // [
    //   x0,y0,z0,
    //   x1,y1,z1,
    //   ...
    // ]
    //
    // amplitudes:
    // [
    //   q0_atom0,
    //   q0_atom1,
    //   ...
    //   q0_atomN,
    //   q1_atom0,
    //   q1_atom1,
    //   ...
    // ]
    //
    // In other words:
    //
    // amplitudes[qIndex * atomCount + atomIndex]
    //
    void prepare_buffers(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& table,
        const std::vector<float>& qValues);


    // CPU curve access
    size_t curve_size() const;

    double q_value(size_t i) const;

    double intensity_value(size_t i) const;


    // WebGPU data export.
    std::vector<float> get_positions() const;

    std::vector<float> get_amplitudes() const;


    // Number of atoms represented in the GPU buffers
    size_t atom_count() const;


private:

    // CPU-generated SAXS curve
    std::vector<double> qValues;
    std::vector<double> intensities;


    // GPU upload buffers

    // positions:
    // [
    //   x0,y0,z0,
    //   x1,y1,z1,
    //   ...
    // ]
    std::vector<float> positions;


    // q-dependent amplitudes:
    //
    // [
    //   q0_atom0,
    //   q0_atom1,
    //   ...
    //   q1_atom0,
    //   q1_atom1,
    //   ...
    // ]
    //
    // Indexed as:
    //
    // amplitudes[qIndex * atomCount + atomIndex]
    std::vector<float> amplitudes;

    std::vector<ShellScatterer> shellScatterers;
};