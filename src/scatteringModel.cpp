#include "scatteringModel.hpp"

#include <cmath>
#include <iostream>

void ScatteringModel::compute_amplitudes(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table,
    double q
)
{
    amplitudes_.clear();
    amplitudes_.reserve(atoms.size());

    for (const auto& atom : atoms)
    {
        if (!table.contains(atom.Z))
        {
            amplitudes_.push_back(0.0);
            continue;
        }

        const auto& ff =
            table[atom.Z];

        const double atomic =
            ff.evaluate(q);

        const double excluded =
            ff.excluded_amplitude(q);

        amplitudes_.push_back(
            atomic -
            solventDensity_ * excluded
        );
    }

    shellScatterers_.clear();

    HydrationShell shell;

    auto shellPoints =
        shell.generate(
            atoms,
            table
        );

    shellScatterers_.reserve(
        shellPoints.size()
    );
    for (const auto& point : shellPoints)
        {
            ShellScatterer scatterer;

            scatterer.x = point.x;
            scatterer.y = point.y;
            scatterer.z = point.z;

            //
            // Temporary constant amplitude.
            // We'll replace this with a proper water form factor.
            //

            scatterer.amplitude = shellContrast_;

            shellScatterers_.push_back(scatterer);
        }
    static bool printed = false;

    if (!printed)
    {
        std::cout
            << "Hydration scatterers = "
            << shellScatterers_.size()
            << "\n";

        printed = true;
    }
    std::cout << "\n";
    std::cout << "Scattering model @ q = "
            << q
            << "\n";

    std::cout
        << "Atom\tZ\tAtomic\tExcluded\tEffective\n";

    double totalAtomic = 0.0;
    double totalExcluded = 0.0;
    double totalEffective = 0.0;

    for (size_t i = 0;
        i < std::min<size_t>(5, atoms.size());
        ++i)
    {
        const auto& ff =
            table[atoms[i].Z];

        const double atomic =
            ff.evaluate(q);

        const double excluded =
            solventDensity_ *
            ff.excluded_amplitude(q);

        const double effective =
            amplitudes_[i];

        std::cout
            << i << "\t"
            << atoms[i].Z << "\t"
            << atomic << "\t"
            << excluded << "\t"
            << effective
            << "\n";
    }

    for (size_t i = 0;
        i < atoms.size();
        ++i)
    {
        if (!table.contains(atoms[i].Z))
            continue;

        const auto& ff =
            table[atoms[i].Z];

        totalAtomic +=
            ff.evaluate(q);

        totalExcluded +=
            solventDensity_ *
            ff.excluded_amplitude(q);

        totalEffective +=
            amplitudes_[i];
    }

    std::cout << "\n";
    std::cout
        << "Atomic total    : "
        << totalAtomic
        << "\n";

    std::cout
        << "Excluded total  : "
        << totalExcluded
        << "\n";

    std::cout
        << "Effective total : "
        << totalEffective
        << "\n\n";
}

const std::vector<double>&
ScatteringModel::amplitudes() const
{
    return amplitudes_;
}

const std::vector<ShellScatterer>&
ScatteringModel::shell_scatterers() const
{
    return shellScatterers_;
}