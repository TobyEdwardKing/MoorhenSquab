#ifndef HYDRATIONSHELL_HPP
#define HYDRATIONSHELL_HPP

#include <vector>

#include "atomExtractor.hpp"
#include "formFactor.hpp"

struct HydrationPoint
{
    double x;
    double y;
    double z;
};

class HydrationShell
{
public:

    // Generate solvent-accessible hydration shell
    std::vector<HydrationPoint> generate(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& table
    ) const;

private:

    static constexpr double WATER_RADIUS = 1.4;

    // Candidate directions on a unit sphere
    static std::vector<std::array<double,3>> directions;


};

#endif