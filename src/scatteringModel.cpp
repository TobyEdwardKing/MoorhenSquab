
#include "scatteringModel.hpp"

#include <cmath>
#include <complex>
#include <iostream>
#include <algorithm>
#include <stdexcept>

namespace
{
    constexpr double PI = 3.14159265358979323846;

    // ------------------------------------------------------------
    // Associated Legendre polynomial P_l^m(x)
    // ------------------------------------------------------------
    double associated_legendre(
        int l,
        int m,
        double x)
    {
        if (m < 0 || m > l)
            return 0.0;

        x = std::max(-1.0, std::min(1.0, x));

        // P_m^m
        double pmm = 1.0;

        if (m > 0)
        {
            double somx2 =
                std::sqrt(
                    std::max(0.0, (1.0 - x) * (1.0 + x))
                );

            double fact = 1.0;

            for (int i = 1; i <= m; ++i)
            {
                pmm *=
                    -(fact) *
                    somx2;

                fact += 2.0;
            }
        }

        if (l == m)
            return pmm;

        // P_{m+1}^m
        double pmmp1 =
            x * (2.0 * m + 1.0) * pmm;

        if (l == m + 1)
            return pmmp1;

        // Recurrence for higher l
        double pll = 0.0;

        for (int ll = m + 2; ll <= l; ++ll)
        {
            pll =
                (
                    (2.0 * ll - 1.0) * x * pmmp1
                    - (ll + m - 1.0) * pmm
                )
                / (ll - m);

            pmm = pmmp1;
            pmmp1 = pll;
        }

        return pmmp1;
    }


    // ------------------------------------------------------------
    // Spherical harmonic Y_lm(theta, phi)
    //
    // Uses the Condon-Shortley convention.
    // ------------------------------------------------------------
    std::complex<double> spherical_harmonic(
        int l,
        int m,
        double theta,
        double phi)
    {
        if (m < 0)
        {
            const int absM = -m;

            const auto positive =
                spherical_harmonic(
                    l,
                    absM,
                    theta,
                    phi
                );

            return
                std::pow(-1.0, absM) *
                std::conj(positive);
        }

        const double x =
            std::cos(theta);

        const double P =
            associated_legendre(
                l,
                m,
                x
            );

        const double normalisation =
            std::sqrt(
                (
                    (2.0 * l + 1.0)
                    / (4.0 * PI)
                )
                *
                std::tgamma(l - m + 1.0)
                /
                std::tgamma(l + m + 1.0)
            );


        const double angle =
            static_cast<double>(m) * phi;

        return
            normalisation *
            P *
            std::complex<double>(
                std::cos(angle),
                std::sin(angle)
            );

    }


    // ------------------------------------------------------------
    // Spherical Bessel function j_l(x)
    // ------------------------------------------------------------
    double spherical_bessel_j(
    int l,
    double x
)
{
    if (l < 0)
        return 0.0;

    const double ax = std::abs(x);

    if (ax < 1.0e-12)
    {
        return l == 0 ? 1.0 : 0.0;
    }

    /*
     * Small-x power series.
     *
     * j_l(x) =
     *
     * x^l / (2l+1)!!
     *
     * times a convergent correction series.
     */
    if (ax < 1.0)
    {
        double prefactor = 1.0;

        for (int k = 1; k <= l; ++k)
        {
            prefactor *=
                x /
                static_cast<double>(2 * k + 1);
        }

        double sum = 1.0;
        double term = 1.0;
        const double x2 = x * x;

        for (int k = 1; k < 100; ++k)
        {
            term *=
                -x2 /
                (
                    2.0 *
                    static_cast<double>(k) *
                    static_cast<double>(2 * l + 2 * k + 1)
                );

            sum += term;

            if (std::abs(term) <
                std::abs(sum) * 1.0e-15)
            {
                break;
            }
        }

        return prefactor * sum;
    }

    /*
     * Stable downward recurrence.
     */
    const int nmax = l + 40;

    double jp1 = 0.0;
    double jn = 1.0;

    double jl = 0.0;
    double j0down = 0.0;

    for (int n = nmax; n >= 1; --n)
    {
        const double jm1 =
            (
                static_cast<double>(2 * n + 1) / x
            ) * jn - jp1;

        jp1 = jn;
        jn = jm1;

        if (n - 1 == l)
            jl = jm1;

        if (n - 1 == 0)
            j0down = jm1;
    }

    const double j0 =
        std::sin(x) / x;

    return jl * (j0 / j0down);
}


    // ------------------------------------------------------------
    // Multipole coefficient indexing
    //
    // l = 0:
    //   m = 0       -> 0
    //
    // l = 1:
    //   m = -1..1   -> 1..3
    //
    // l = 2:
    //   m = -2..2   -> 4..8
    // ------------------------------------------------------------
    int lm_index(
        int l,
        int m)
    {
        return l * l + l + m;
    }
}


// ============================================================================
// Existing atomic form-factor amplitudes
// ============================================================================

void ScatteringModel::compute_amplitudes(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table,
    double q)
{
    amplitudes_.clear();
    amplitudes_.reserve(atoms.size());

    for (const auto& atom : atoms)
    {
        const FormFactor* ff = nullptr;

        if (!atom.atomType.empty() &&
            table.contains_atom_type(atom.atomType))
        {
            ff =
                &table.get_atom_type(
                    atom.atomType
                );
        }
        else if (table.contains_element(atom.Z))
        {
            ff =
                &table.get_element(atom.Z);
        }
        else
        {
            amplitudes_.push_back(0.0);
            continue;
        }

        const double atomic =
            ff->evaluate(q);

        amplitudes_.push_back(atomic);
    }
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


// ============================================================================
// Atomic multipoles A_lm(q)
// ============================================================================

void ScatteringModel::compute_multipoles(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table,
    const std::vector<double>& qGrid,
    int L)
{
    L_ = L;

    const int nCoefficients =
        (L + 1) * (L + 1);

    multipoles_.clear();
    multipoles_.resize(qGrid.size());

    for (auto& coefficients : multipoles_)
    {
        coefficients.resize(
            nCoefficients,
            std::complex<double>(0.0, 0.0)
        );
    }


    // ------------------------------------------------------------
    // Centre the molecule on its geometric centre.
    //
    // The multipole expansion depends on the choice of origin,
    // so A and C must ultimately use the same origin.
    // ------------------------------------------------------------

    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;

    if (!atoms.empty())
    {
        for (const auto& atom : atoms)
        {
            cx += atom.x;
            cy += atom.y;
            cz += atom.z;
        }

        cx /= static_cast<double>(atoms.size());
        cy /= static_cast<double>(atoms.size());
        cz /= static_cast<double>(atoms.size());
    }


    // ------------------------------------------------------------
    // q loop
    // ------------------------------------------------------------

    for (size_t iq = 0;
         iq < qGrid.size();
         ++iq)
    {
        const double q =
            qGrid[iq];

        auto& coefficients =
            multipoles_[iq];


        // --------------------------------------------------------
        // Atom loop
        // --------------------------------------------------------

        for (const auto& atom : atoms)
        {
            // Coordinates relative to molecular centre
            const double x =
                static_cast<double>(atom.x) - cx;

            const double y =
                static_cast<double>(atom.y) - cy;

            const double z =
                static_cast<double>(atom.z) - cz;


            const double r =
                std::sqrt(
                    x * x +
                    y * y +
                    z * z
                );


            double theta = 0.0;
            double phi = 0.0;

            if (r > 1e-12)
            {
                theta =
                    std::acos(
                        std::max(
                            -1.0,
                            std::min(
                                1.0,
                                z / r
                            )
                        )
                    );

                phi =
                    std::atan2(y, x);
            }


            // ----------------------------------------------------
            // Select the same X-ray form factor used elsewhere.
            //
            // PEPSI atom type first.
            // Elemental fallback second.
            // ----------------------------------------------------

            const FormFactor* ff = nullptr;

            if (!atom.atomType.empty() &&
                table.contains_atom_type(atom.atomType))
            {
                ff =
                    &table.get_atom_type(
                        atom.atomType
                    );
            }
            else if (table.contains_element(atom.Z))
            {
                ff =
                    &table.get_element(atom.Z);
            }
            else
            {
                continue;
            }


            const double f =
                ff->evaluate(q);


            // ----------------------------------------------------
            // l loop
            // ----------------------------------------------------

            for (int l = 0;
                 l <= L;
                 ++l)
            {
                const double jl =
                    spherical_bessel_j(
                        l,
                        q * r
                    );


                // i^l
                const std::complex<double> phase =
                    std::pow(
                        std::complex<double>(
                            0.0,
                            1.0
                        ),
                        l
                    );


                // ------------------------------------------------
                // m loop
                // ------------------------------------------------

                for (int m = -l;
                     m <= l;
                     ++m)
                {
                    const std::complex<double> Y =
                        spherical_harmonic(
                            l,
                            m,
                            theta,
                            phi
                        );


                    const int index =
                        lm_index(l, m);


                    coefficients[index] +=
                        4.0 * PI *
                        phase *
                        f *
                        jl *
                        std::conj(Y);
                }
            }
        }
    }


    std::cout
        << "Atomic multipoles calculated:\n"
        << "  Atoms: "
        << atoms.size()
        << "\n"
        << "  q points: "
        << qGrid.size()
        << "\n"
        << "  L: "
        << L
        << "\n"
        << "  coefficients per q: "
        << nCoefficients
        << "\n"
        << "  Centre: ("
        << cx << ", "
        << cy << ", "
        << cz << ")\n";
}


// ============================================================================
// Intensity reconstructed from atomic multipoles
// ============================================================================

double ScatteringModel::intensity_from_multipoles(
    size_t qIndex,
    int maxL
) const
{
    if (qIndex >= multipoles_.size())
    {
        throw std::out_of_range(
            "Multipole q index out of range"
        );
    }

    if (maxL < 0)
    {
        throw std::invalid_argument(
            "Multipole L must be non-negative"
        );
    }

    if (maxL > L_)
    {
        throw std::invalid_argument(
            "Requested multipole L exceeds calculated L"
        );
    }

    double intensity = 0.0;

    for (int l = 0; l <= maxL; ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            const int index =
                lm_index(l, m);

            intensity +=
                std::norm(
                    multipoles_[qIndex][index]
                );
        }
    }

    return intensity / (4.0 * M_PI);
}

double ScatteringModel::combined_intensity_from_multipoles(
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
) const
{
    if (qIndex >= multipoles_.size() ||
        qIndex >= excludedVolumeCoefficients.size() ||
        qIndex >= hydrationCoefficients.size())
    {
        throw std::out_of_range(
            "Combined intensity q index out of range"
        );
    }

    if (maxL < 0 || maxL > L_)
    {
        throw std::invalid_argument(
            "Invalid L for combined intensity"
        );
    }

    const auto& atomic = multipoles_[qIndex];
    const auto& excluded = excludedVolumeCoefficients[qIndex];
    const auto& hydration = hydrationCoefficients[qIndex];

    const std::size_t requiredCoefficients =
        static_cast<std::size_t>((maxL + 1) * (maxL + 1));

    if (atomic.size() < requiredCoefficients ||
        excluded.size() < requiredCoefficients ||
        hydration.size() < requiredCoefficients)
    {
        throw std::runtime_error(
            "Insufficient multipole coefficients for combined intensity"
        );
    }

    if (qIndex == 0)
    {
        double atomicNorm = 0.0;
        double excludedNorm = 0.0;
        double hydrationNorm = 0.0;
        double combinedNorm = 0.0;

        for (int l = 0; l <= maxL; ++l)
        {
            for (int m = -l; m <= l; ++m)
            {
                const std::size_t index =
                    static_cast<std::size_t>(lm_index(l, m));

                const auto a = atomic[index];
                const auto v = solventDensity * excluded[index];
                const auto b = shellContrast * hydration[index];
                const auto total = a - v + b;
                if (l == 0 && m == 0)
                {
                    std::cout << "A00: " << a
                            << "\nV00 (scaled): " << v
                            << "\nB00 (scaled): " << b
                            << "\nTotal C00: " << total
                            << '\n';
                }
                atomicNorm += std::norm(a);
                excludedNorm += std::norm(v);
                hydrationNorm += std::norm(b);
                combinedNorm += std::norm(total);
            }
        }

        std::cout << "\n=== COEFFICIENT-LEVEL DIAGNOSTIC ===\n";
        std::cout << "qIndex: " << qIndex << '\n';
        std::cout << "maxL: " << maxL << '\n';
        std::cout << "Atomic intensity: "
                << atomicNorm / (4.0 * PI) << '\n';
        std::cout << "Scaled excluded-volume intensity: "
                << excludedNorm / (4.0 * PI) << '\n';
        std::cout << "Scaled hydration intensity: "
                << hydrationNorm / (4.0 * PI) << '\n';
        std::cout << "Combined intensity: "
                << combinedNorm / (4.0 * PI) << '\n';
        std::cout << "===================================\n";
    }

    double intensity = 0.0;

    for (int l = 0; l <= maxL; ++l)
    {
        for (int m = -l; m <= l; ++m)
        {
            const std::size_t index =
                static_cast<std::size_t>(lm_index(l, m));

            const std::complex<double> total =
                atomic[index]
                - solventDensity * excluded[index]
                + shellContrast * hydration[index];

            intensity += std::norm(total);
        }
    }

    return intensity / (4.0 * PI);
}


// ============================================================================
// Multipole accessor
// ============================================================================

const std::vector<
    std::vector<std::complex<double>>
>&
ScatteringModel::multipoles() const
{
    return multipoles_;
}
