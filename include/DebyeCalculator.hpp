#pragma once

#include <vector>

#include "atomExtractor.hpp"
#include "formFactor.hpp"

#include "scatteringModel.hpp"

class DebyeCalculator
{
public:

    // Generate CPU SAXS curve
    void compute_curve(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& formFactors,
        double qMin,
        double qMax,
        double qStep);


    // Prepare contiguous arrays for GPU upload.
    //
    // positions:
    // [
    //   x0,y0,z0,
    //   x1,y1,z1,
    //   ...
    // ]
    //
    // formFactors:
    // [
    //   f0,
    //   f1,
    //   ...
    // ]
    void prepare_buffers(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& table,
        double q);


    // CPU curve access
    size_t curve_size() const;

    double q_value(size_t i) const;

    double intensity_value(size_t i) const;


    // WebGPU data export.
    //
    // These return copies that Emscripten can transfer
    // safely into JavaScript.
    std::vector<float> get_positions() const;

    std::vector<float> get_form_factors() const;


    // Number of atoms represented in the GPU buffers
    size_t atom_count() const;


private:

    // CPU-generated SAXS curve
    std::vector<double> qValues;

    std::vector<double> intensities;


    // GPU upload buffers
    //
    // positions:
    // [
    //   x0,y0,z0,
    //   x1,y1,z1,
    //   ...
    // ]
    std::vector<float> positions;


    // form factors:
    // [
    //   f0,
    //   f1,
    //   ...
    // ]
    std::vector<float> formFactors;

    std::vector<ShellScatterer> shellScatterers;
    
    // Single-q Debye calculation
    double calculate_point(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& formFactors,
        double q) const;

    double calculate_shell_shell(
        const std::vector<ShellScatterer>& shell,
        double q
    ) const;

    double calculate_atom_shell(
    const std::vector<AtomData>& atoms,
    const std::vector<ShellScatterer>& shell,
    const std::vector<double>& amplitudes,
    double q
    ) const;
    
};