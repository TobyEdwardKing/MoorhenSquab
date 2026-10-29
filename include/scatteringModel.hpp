#ifndef SCATTERINGMODEL_HPP
#define SCATTERINGMODEL_HPP

#include <complex>
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

    // Calculate ordinary atomic scattering amplitudes
    void compute_amplitudes(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& table,
        double q
    );

    const std::vector<double>& amplitudes() const;

    const std::vector<ShellScatterer>& shell_scatterers() const;

    // Calculate spherical-harmonic multipole coefficients:
    //
    // A_lm(q) =
    //   4 pi i^l sum_i f_i(q) j_l(q r_i) Y_lm*(theta_i, phi_i)
    //
    // Stored as:
    // multipoles_[qIndex][lm_index]
    void compute_multipoles(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& table,
        const std::vector<double>& qGrid,
        int L
    );

    // Calculate I_aa(q) using only terms with l <= maxL.
    //
    // The multipoles themselves may have been calculated at a
    // higher L. This allows several truncated curves to be
    // reconstructed from one high-L calculation.
    double intensity_from_multipoles(
        size_t qIndex,
        int maxL
    ) const;

    // Access the calculated multipole coefficients.
    const std::vector<
        std::vector<std::complex<double>>
    >& multipoles() const;

    double combined_intensity_from_multipoles(
        size_t qIndex,
        int maxL,
        const std::vector<
            std::vector<std::complex<double>>
        >& excludedVolumeCoefficients,
        double solventDensity,
        const std::vector<
            std::vector<std::complex<double>>
        >& hydrationCoefficients,
        double shellContrast
    ) const;
private:

    // Effective atomic amplitudes
    std::vector<double> amplitudes_;

    // Hydration shell scatterers
    std::vector<ShellScatterer> shellScatterers_;

    // Bulk solvent electron density
    double solventDensity_ = 0.334;

    // Excess shell density
    double shellContrast_ = 0.03;

    // Maximum spherical-harmonic angular momentum
    int L_ = 0;

    // Multipole coefficients:
    //
    // multipoles_[qIndex][lmIndex]
    //
    // Number of coefficients per q:
    // (L + 1)^2
    std::vector<
        std::vector<std::complex<double>>
    > multipoles_;
};

#endif