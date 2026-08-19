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

    // for (const auto& atom : atoms)
    //     {
    //         const FormFactor* ff = nullptr;

    //         if (!atom.atomType.empty() &&
    //             table.contains_atom_type(atom.atomType))
    //         {
    //             ff = &table.get_atom_type(atom.atomType);
    //         }
    //         else if (table.contains_element(atom.Z))
    //         {
    //             ff = &table.get_element(atom.Z);
    //         }
    //         else
    //         {
    //             amplitudes_.push_back(0.0);
    //             continue;
    //         }

    //         const double atomic =
    //             ff->evaluate(q);

    //         // Excluded volume still comes from the element table
    //         if (!table.contains_element(atom.Z))
    //         {
    //             amplitudes_.push_back(0.0);
    //             continue;
    //         }

    //         const auto& elementFF =
    //             table.get_element(atom.Z);

    //         const double excluded =
    //             solventDensity_ *
    //             elementFF.excluded_amplitude(q);

    //         amplitudes_.push_back(
    //             atomic -
    //             excluded
    //         );
    //     }
    for (const auto& atom : atoms)
        {
            const FormFactor* ff = nullptr;

            // Protein atoms with a recognised atom type use
            // the protein/SAXS atom-type form factor.
            if (!atom.atomType.empty() &&
                table.contains_atom_type(atom.atomType))
            {
                ff = &table.get_atom_type(atom.atomType);
            }
            // Everything else falls back to the elemental form factor.
            // This includes ligand atoms.
            else if (table.contains_element(atom.Z))
            {
                ff = &table.get_element(atom.Z);
            }
            else
            {
                amplitudes_.push_back(0.0);
                continue;
            }

            const double atomic =
                ff->evaluate(q);

            // Elemental form factor is required for excluded volume.
            if (!table.contains_element(atom.Z))
            {
                amplitudes_.push_back(0.0);
                continue;
            }

            const auto& elementFF =
                table.get_element(atom.Z);

            const double excluded =
                solventDensity_ *
                elementFF.excluded_amplitude(q);

            amplitudes_.push_back(
                atomic - excluded
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

    const auto& hydrogen = table.get_element(0);
    const auto& oxygen   = table.get_element(7);
    const double waterAmplitude =
        oxygen.evaluate(q)
        + 2.0 * hydrogen.evaluate(q);

    for (const auto& point : shellPoints)
    {
        ShellScatterer scatterer;

        scatterer.x = point.x;
        scatterer.y = point.y;
        scatterer.z = point.z;

        scatterer.amplitude =
            shellContrast_ * waterAmplitude;

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
        const FormFactor* ff = nullptr;

        if (!atoms[i].atomType.empty() &&
            table.contains_atom_type(atoms[i].atomType))
        {
            ff = &table.get_atom_type(atoms[i].atomType);
        }
        else if (table.contains_element(atoms[i].Z))
        {
            ff = &table.get_element(atoms[i].Z);
        }
        else
        {
            continue;
        }
        const double atomic =
            ff->evaluate(q);

        const auto& elementFF =
            table.get_element(atoms[i].Z);

        const double excluded =
            solventDensity_ *
            elementFF.excluded_amplitude(q);

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

    // for (size_t i = 0;
    //     i < atoms.size();
    //     ++i)
    // {
    //     if (!table.contains_element(atoms[i].Z))
    //         continue;

    //     const FormFactor* ff = nullptr;

    //     if (!atoms[i].atomType.empty() &&
    //         table.contains_atom_type(atoms[i].atomType))
    //     {
    //         ff = &table.get_atom_type(atoms[i].atomType);
    //     }
    //     else if (table.contains_element(atoms[i].Z))
    //     {
    //         ff = &table.get_element(atoms[i].Z);
    //     }
    //     else
    //     {
    //         continue;
    //     }

    //     totalAtomic +=
    //         ff->evaluate(q);
    //     if (!table.contains_element(atoms[i].Z))
    //     {
    //         amplitudes_.push_back(0.0);
    //         continue;
    //     }

    //     const auto& elementFF =
    //         table.get_element(atoms[i].Z);

    //     totalExcluded +=
    //         solventDensity_ *
    //         elementFF.excluded_amplitude(q);

    //     totalEffective +=
    //         amplitudes_[i];
    // }

    for (size_t i = 0;
            i < atoms.size();
            ++i)
        {
            if (!table.contains_element(atoms[i].Z))
                continue;

            const FormFactor* ff = nullptr;

            if (!atoms[i].atomType.empty() &&
                table.contains_atom_type(atoms[i].atomType))
            {
                ff = &table.get_atom_type(atoms[i].atomType);
            }
            else
            {
                // Ligands and other atoms without a recognised
                // protein atom type use the elemental form factor.
                ff = &table.get_element(atoms[i].Z);
            }

            totalAtomic +=
                ff->evaluate(q);

            const auto& elementFF =
                table.get_element(atoms[i].Z);

            totalExcluded +=
                solventDensity_ *
                elementFF.excluded_amplitude(q);

            totalEffective +=
                amplitudes_[i];
        }
        
    double totalShellAmplitude = 0.0;

    for (const auto& s : shellScatterers_)
    {
        totalShellAmplitude += s.amplitude;
    }

    const double totalGPUAmplitude =
        totalEffective + totalShellAmplitude;

    const double expectedI0 =
        totalGPUAmplitude * totalGPUAmplitude;
    std::cout << "\n";

    std::cout
        << "===== Scatterer preparation =====\n";

    std::cout
        << "Atoms              : "
        << atoms.size()
        << "\n";

    std::cout
        << "Hydration points   : "
        << shellScatterers_.size()
        << "\n\n";

    std::cout
        << "Atomic amplitude sum : "
        << totalAtomic
        << "\n";

    std::cout
        << "Solvent exclusion   : "
        << totalExcluded
        << "\n";

    std::cout
        << "Effective atom sum  : "
        << totalEffective
        << "\n";

    std::cout
        << "Hydration sum       : "
        << totalShellAmplitude
        << "\n";

    std::cout
        << "Total GPU amplitude : "
        << totalGPUAmplitude
        << "\n";

    std::cout
        << "Expected I(0)       : "
        << expectedI0
        << "\n";

    std::cout
        << "===============================\n\n";
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