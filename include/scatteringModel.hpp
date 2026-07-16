#ifndef SCATTERINGMODEL_HPP
#define SCATTERINGMODEL_HPP

#include <vector>

#include "atomExtractor.hpp"
#include "formFactor.hpp"
#include "hydrationShell.hpp"

struct ShellScatterer
{
    double x;
    double y;
    double z;

    double amplitude;
};

class ScatteringModel
{
public:

    void compute_amplitudes(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& table,
        double q
    );

    const std::vector<double>& amplitudes() const;
    const std::vector<ShellScatterer>& shell_scatterers() const;
private:

    // Effective atomic amplitudes
    std::vector<double> amplitudes_;

    // Hydration shell scatterers
    std::vector<ShellScatterer> shellScatterers_;

    // Bulk solvent electron density
    double solventDensity_ = 0.334;

    // Excess shell density
    double shellContrast_ = 0.03;
};

#endif