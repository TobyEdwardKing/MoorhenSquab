#include "DebyeCalculator.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include "scatteringModel.hpp"

void DebyeCalculator::prepare_buffers(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table,
    double q)
{
    positions.clear();
    formFactors.clear();

    positions.reserve(atoms.size() * 3);
    formFactors.reserve(atoms.size());

    //
    // Compute amplitudes ONCE
    //

    ScatteringModel model;

    model.compute_amplitudes(
        atoms,
        table,
        q
    );


    const auto& amplitudes =
        model.amplitudes();


    shellScatterers =
        model.shell_scatterers();

    //
    // Fill buffers
    //

    for (size_t i = 0; i < atoms.size(); ++i)
    {
        positions.push_back(atoms[i].x);
        positions.push_back(atoms[i].y);
        positions.push_back(atoms[i].z);

        formFactors.push_back(
            static_cast<float>(
                amplitudes[i]
            )
        );
    }
}


// These two functions are the WASM -> JavaScript boundary.
// They return copies that Emscripten can safely expose.

std::vector<float> DebyeCalculator::get_positions() const
{
    return positions;
}


std::vector<float> DebyeCalculator::get_form_factors() const
{
    return formFactors;
}


size_t DebyeCalculator::atom_count() const
{
    return formFactors.size();
}


double DebyeCalculator::calculate_point(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table,
    double q) const
{
    std::cout
        << "DebyeCalculator called\n";


    std::cout
        << "Number of atoms: "
        << atoms.size()
        << "\n";


    std::cout
        << "q value: "
        << q
        << "\n";


    double intensity = 0.0;


    const size_t n =
        atoms.size();



    // Precompute atomic amplitudes

    std::vector<double> f(n);
    std::vector<double> g(n);

    for (size_t i = 0; i < n; i++)
    {
        if (table.contains(atoms[i].Z))
        {
            const auto& ff = table[atoms[i].Z];

            f[i] = ff.evaluate(q);
            g[i] = ff.excluded_amplitude(q);
        }
        else
        {
            f[i] = 0.0;
            g[i] = 0.0;
        }
    }

    // Symmetry:
    // calculate only j >= i
    // double off-diagonal terms

    for (size_t i = 0; i < n; i++)
    {
        const auto& atom_i =
            atoms[i];


        const double fi =
            f[i];


        if (fi == 0.0)
            continue;



        for (size_t j = i; j < n; j++)
        {
            const auto& atom_j =
                atoms[j];


            const double fj =
                f[j];


            if (fj == 0.0)
                continue;



            const double dx =
                atom_i.x - atom_j.x;


            const double dy =
                atom_i.y - atom_j.y;


            const double dz =
                atom_i.z - atom_j.z;



            const double r =
                std::sqrt(
                    dx * dx +
                    dy * dy +
                    dz * dz);



            const double qr =
                q * r;



            const double sinc =
                std::abs(qr) < 1e-12
                    ? 1.0
                    : std::sin(qr) / qr;



            const double gi = g[i];
            const double gj = g[j];

            const double effective_i = fi - gi;
            const double effective_j = fj - gj;

            const double contribution =
                effective_i *
                effective_j *
                sinc;



            if (i == j)
            {
                intensity += contribution;
            }
            else
            {
                intensity +=
                    2.0 * contribution;
            }
        }
    }


    std::vector<double> effectiveAmplitudes(
        atoms.size()
    );


    for(size_t i=0;i<atoms.size();i++)
    {
        if(table.contains(atoms[i].Z))
        {
            const auto& ff =
                table[atoms[i].Z];

            effectiveAmplitudes[i] =
                ff.evaluate(q)
                -
                0.334 *
                ff.excluded_amplitude(q);
        }
        else
        {
            effectiveAmplitudes[i]=0.0;
        }
    }


    double shellShell =
        calculate_shell_shell(
            shellScatterers,
            q
        );


    double atomShell =
        calculate_atom_shell(
            atoms,
            shellScatterers,
            effectiveAmplitudes,
            q
        );





    static bool printed = false;

    if (!printed)
    {
        std::cout
            << "\n===== Hydration diagnostics =====\n";

        std::cout
            << "Atoms: "
            << atoms.size()
            << "\n";

        std::cout
            << "Shell points: "
            << shellScatterers.size()
            << "\n";

        if (!shellScatterers.empty())
        {
            std::cout
                << "First shell amplitude: "
                << shellScatterers[0].amplitude
                << "\n";
        }

        std::cout
            << "Atom-atom: "
            << intensity
            << "\n";

        std::cout
            << "Atom-shell: "
            << atomShell
            << "\n";

        std::cout
            << "Shell-shell: "
            << shellShell
            << "\n";

        std::cout
            << "Total: "
            << intensity + atomShell + shellShell
            << "\n";

        std::cout
            << "================================\n";

        printed = true;
    }

    return intensity + shellShell + atomShell;
}



void DebyeCalculator::compute_curve(

    const std::vector<AtomData>& atoms,
    const FormFactorTable& formFactors,
    double qMin,
    double qMax,
    double qStep)
{
    std::cout << "compute_curve()\n";
    std::cout << "qMin = " << qMin << "\n";
    std::cout << "qMax = " << qMax << "\n";
    std::cout << "qStep = " << qStep << "\n";
    const size_t nPoints =
        static_cast<size_t>(
            (qMax - qMin) / qStep)
        + 1;

    std::cout << "nPoints = "
            << nPoints
            << "\n";

    qValues.clear();
    intensities.clear();



    qValues.reserve(nPoints);
    intensities.reserve(nPoints);



    for (size_t i = 0; i < nPoints; ++i)
    {
        const double q =
            qMin +
            static_cast<double>(i)
            *
            qStep;



        qValues.push_back(q);



        intensities.push_back(
            calculate_point(
                atoms,
                formFactors,
                q));
    }
    std::cout
    << "Finished compute_curve(). "
    << "qValues = " << qValues.size()
    << ", intensities = " << intensities.size()
    << "\n";
}


double DebyeCalculator::calculate_shell_shell(
    const std::vector<ShellScatterer>& shell,
    double q
) const
{
    double intensity = 0.0;

    const size_t n = shell.size();


    for (size_t i = 0; i < n; i++)
    {
        const auto& si = shell[i];


        for (size_t j = i; j < n; j++)
        {
            const auto& sj = shell[j];


            const double dx =
                si.x - sj.x;

            const double dy =
                si.y - sj.y;

            const double dz =
                si.z - sj.z;


            const double r =
                std::sqrt(
                    dx*dx +
                    dy*dy +
                    dz*dz
                );


            const double qr =
                q * r;


            const double sinc =
                std::abs(qr) < 1e-12
                    ? 1.0
                    : std::sin(qr) / qr;


            const double contribution =
                si.amplitude *
                sj.amplitude *
                sinc;


            if (i == j)
            {
                intensity += contribution;
            }
            else
            {
                intensity +=
                    2.0 * contribution;
            }
        }
    }


    return intensity;
}


double DebyeCalculator::calculate_atom_shell(
    const std::vector<AtomData>& atoms,
    const std::vector<ShellScatterer>& shell,
    const std::vector<double>& amplitudes,
    double q
) const
{
    double intensity = 0.0;


    for (size_t i = 0; i < atoms.size(); i++)
    {
        const double fi = amplitudes[i];


        for (const auto& s : shell)
        {
            const double dx =
                atoms[i].x - s.x;

            const double dy =
                atoms[i].y - s.y;

            const double dz =
                atoms[i].z - s.z;


            const double r =
                std::sqrt(
                    dx*dx +
                    dy*dy +
                    dz*dz
                );


            const double qr =
                q*r;


            const double sinc =
                std::abs(qr)<1e-12
                ? 1.0
                : std::sin(qr)/qr;


            intensity +=
                2.0 *
                fi *
                s.amplitude *
                sinc;
        }
    }


    return intensity;
}


size_t DebyeCalculator::curve_size() const
{
    return intensities.size();
}



double DebyeCalculator::q_value(size_t i) const
{
    return qValues[i];
}



double DebyeCalculator::intensity_value(size_t i) const
{
    return intensities[i];
}