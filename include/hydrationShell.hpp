
#ifndef HYDRATION_SHELL_HPP
#define HYDRATION_SHELL_HPP

#include <cstddef>
#include <vector>
#include <complex>

#include "atomExtractor.hpp"
#include "formFactor.hpp"


// One grid cell representing a small volume of the
// hydration shell. Coordinates are in Angstroms and
// volume is in Angstroms cubed.
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
    // Generate a grid representation of the hydration shell.
    //
    // shellThickness: thickness outside the molecular surface (A).
    // gridSpacing:    distance between grid-cell centres (A).
    //
    // Defaults are starting values, not fitted parameters.
    std::vector<HydrationPoint> generate(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& table,
        double shellThickness = 5.0,
        double gridSpacing = 1.0
    ) const;

    // Sum the represented volume of the generated shell (A^3).
    static double total_volume(
        const std::vector<HydrationPoint>& shell
    );
    
    void compute_multipoles(
        const std::vector<HydrationPoint>& shellPoints,
        const std::vector<double>& qGrid,
        const std::vector<double>& centre,
        int L
    );

    const std::vector<std::vector<std::complex<double>>>&
    coefficients() const
    {
        return multipoles_;
    }
private:
    struct AtomSphere
    {
        double x;
        double y;
        double z;
        double radius;
    };
    std::vector<std::vector<std::complex<double>>> multipoles_;
    int L_ = 0;
};

#endif 