#ifndef HYDRATION_SHELL_HPP
#define HYDRATION_SHELL_HPP

#include <cstddef>
#include <vector>
#include <complex>

#include "atomExtractor.hpp"
#include "formFactor.hpp"

// One grid cell representing a small volume of the hydration shell.
// Coordinates are in Angstroms and volume is in Angstroms cubed.
struct HydrationPoint
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double volume = 0.0;
};

class HydrationShell
{
public:
    // Generate the paper-style coarse hydration-shell grid.
    // The shell width is calculated from the coordinate-derived molecular Rg.
    // gridSpacing is the grid-cell width in Angstroms.
    std::vector<HydrationPoint> generate(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& table,
        double gridSpacing = 3.5
    ) const;

    static double total_volume(const std::vector<HydrationPoint>& shell);

    // Radius of gyration of the shell-point distribution about its centroid.
    static double radius_of_gyration(const std::vector<HydrationPoint>& points);

    double shell_thickness() const { return lastShellThickness_; }
    double grid_spacing() const { return lastGridSpacing_; }

    void compute_multipoles(
        const std::vector<HydrationPoint>& shellPoints,
        const FormFactorTable& table,
        const std::vector<double>& qGrid,
        const std::vector<double>& centre,
        int L,
        bool useInterpolation = true
    );

    // Compare interpolated coefficients with direct shell sums at
    // experimental-grid points near the sparse-grid midpoints.
    void validate_interpolation(
        const std::vector<HydrationPoint>& shellPoints,
        const std::vector<double>& qGrid,
        const std::vector<double>& centre,
        int maxTests = 25
    ) const;

    const std::vector<std::vector<std::complex<double>>>& coefficients() const
    {
        return multipoles_;
    }

    // Compare combined intensity using interpolated hydration coefficients
    // against direct hydration sums. Apply the same excluded-volume radius
    // correction as the production calculation: G = 1 + 3*fraction.
    void validate_combined_intensity_interpolation(
        const std::vector<HydrationPoint>& shellPoints,
        const std::vector<double>& qGrid,
        const std::vector<double>& centre,
        const std::vector<std::vector<std::complex<double>>>& atomic,
        const std::vector<std::vector<std::complex<double>>>& excluded,
        double solventDensity,
        double shellContrast,
        double radiusAdjustmentFraction,
        int maxTests = 25
    ) const;

    double water_form_factor(std::size_t qIndex) const
    {
        return waterFormFactors_.at(qIndex);
    }
private:
    std::vector<std::vector<std::complex<double>>> multipoles_;
    int L_ = 0;
    mutable double lastShellThickness_ = 5.0;
    mutable double lastGridSpacing_ = 3.5;
    std::vector<double> waterFormFactors_;
};

void run_spherical_bessel_test();

#endif
