
#include "excludedVolumeModel.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>

namespace
{
    constexpr double PI = 3.14159265358979323846;

    /*
     * Associated Legendre polynomial P_l^m(x)
     *
     * Uses the standard Condon-Shortley convention,
     * including the (-1)^m factor.
     */
    double associated_legendre(
        int l,
        int m,
        double x)
    {
        if (m < 0 || m > l)
            return 0.0;

        // P_m^m(x)
        double pmm = 1.0;

        if (m > 0)
        {
            const double somx2 =
                std::sqrt(
                    std::max(
                        0.0,
                        (1.0 - x) * (1.0 + x)
                    )
                );

            double fact = 1.0;

            for (int i = 1; i <= m; ++i)
            {
                pmm *= -fact * somx2;
                fact += 2.0;
            }
        }

        if (l == m)
            return pmm;

        // P_{m+1}^m(x)
        double pmmp1 =
            x * (2.0 * m + 1.0) * pmm;

        if (l == m + 1)
            return pmmp1;

        /*
         * Recurrence:
         *
         * P_l^m(x) =
         * [(2l-1)x P_{l-1}^m
         *  - (l+m-1)P_{l-2}^m] / (l-m)
         */
        double pll = 0.0;

        double pprev = pmm;
        double pcurr = pmmp1;

        for (int ll = m + 2; ll <= l; ++ll)
        {
            pll =
                (
                    (2.0 * ll - 1.0) * x * pcurr
                    - (ll + m - 1.0) * pprev
                )
                / (ll - m);

            pprev = pcurr;
            pcurr = pll;
        }

        return pcurr;
    }


    /*
     * Spherical harmonic Y_l^m(theta, phi)
     *
     * For m >= 0:
     *
     * Y_l^m =
     * sqrt[(2l+1)/(4pi) * (l-m)!/(l+m)!]
     * P_l^m(cos(theta)) exp(i m phi)
     *
     * For m < 0:
     *
     * Y_l^{-m} = (-1)^|m| conj(Y_l^|m|)
     */
    std::complex<double> spherical_harmonic(
        int l,
        int m,
        double theta,
        double phi)
    {
        if (std::abs(m) > l)
            return {0.0, 0.0};

        const double cosTheta =
            std::cos(theta);

        if (m >= 0)
        {
            const double P =
                associated_legendre(
                    l,
                    m,
                    cosTheta
                );

            /*
             * Calculate factorial ratio in logarithmic
             * form to avoid unnecessarily large factorials.
             */
            const double logNorm =
                0.5 *
                (
                    std::log(
                        (2.0 * l + 1.0) /
                        (4.0 * PI)
                    )
                    +
                    std::lgamma(l - m + 1.0)
                    -
                    std::lgamma(l + m + 1.0)
                );

            const double norm =
                std::exp(logNorm);

            const double phase =
                static_cast<double>(m) * phi;

            return
                norm *
                P *
                std::complex<double>(
                    std::cos(phase),
                    std::sin(phase)
                );
        }

        /*
         * Negative m:
         *
         * Y_l^{-m} = (-1)^|m| conj(Y_l^|m|)
         */
        const int mp = -m;

        const auto positive =
            spherical_harmonic(
                l,
                mp,
                theta,
                phi
            );

        const double sign =
            (mp % 2 == 0) ? 1.0 : -1.0;

        return sign * std::conj(positive);
    }


    
    /*
     * Spherical Bessel function j_l(x).
     *
     * Uses a small-x power series and downward recurrence,
     * matching the atomic multipole implementation.
     */
    double spherical_bessel_j(
        int l,
        double x)
    {
        if (l < 0)
            return 0.0;

        const double ax = std::abs(x);

        if (ax < 1.0e-12)
        {
            return l == 0 ? 1.0 : 0.0;
        }

        // Small-x power series.
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

        // Stable downward recurrence.
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

        const double j0 = std::sin(x) / x;

        return jl * (j0 / j0down);
    }


} // anonymous namespace


/*
 * Convert (l,m) into a flat coefficient index.
 *
 * Number of coefficients up to L:
 *
 *     (L+1)^2
 */
int ExcludedVolumeModel::lm_index(
    int l,
    int m)
{
    return l * l + l + m;
}


/*
 * Reference Gaussian dummy-atom form factor:
 *
 *     g_i(q) =
 *         V_i exp[-V_i^(2/3) q^2 / (4 pi)]
 *
 * V_i is in Angstrom^3.
 * q is in Angstrom^-1.
 */

double ExcludedVolumeModel::dummy_form_factor(
    double q,
    double volume)
{
    if (volume <= 0.0)
        return 0.0;

    return volume * std::exp(
        -std::pow(volume, 2.0 / 3.0) * q * q
        / (4.0 * PI)
    );
}


/*
 * Calculate excluded-volume multipole coefficients C_lm(q).
 */
void ExcludedVolumeModel::compute_coefficients(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table,
    const std::vector<double>& qGrid,
    int L,
    double radiusAdjustmentFraction)
{
    std::cout << ">>> ENTERED compute_coefficients <<<\n";

    if (L < 0)
        throw std::runtime_error(
            "ExcludedVolumeModel: L must be >= 0"
        );

    L_ = L;

    const int nCoefficients =
        (L_ + 1) * (L_ + 1);

    /*
     * Centre the coordinates using the arithmetic mean
     * of all atom coordinates.
     *
     * This MUST match the centring used by
     * ScatteringModel::compute_multipoles().
     *
     * Do not modify the original atom coordinates:
     * only the temporary coordinates used below are centred.
     */
    double centreX = 0.0;
    double centreY = 0.0;
    double centreZ = 0.0;

    for (const auto& atom : atoms)
    {
        centreX += static_cast<double>(atom.x);
        centreY += static_cast<double>(atom.y);
        centreZ += static_cast<double>(atom.z);
    }

    if (!atoms.empty())
    {
        const double invN =
            1.0 / static_cast<double>(atoms.size());

        centreX *= invN;
        centreY *= invN;
        centreZ *= invN;
    }

    std::cout
        << "ExcludedVolumeModel coordinate centre: ("
        << centreX << ", "
        << centreY << ", "
        << centreZ << ") A\n";


    /*
     * Allocate:
     *
     * coefficients_[qIndex][lmIndex]
     */
    coefficients_.clear();

    coefficients_.resize(
        qGrid.size(),
        std::vector<std::complex<double>>(
            nCoefficients,
            std::complex<double>(0.0, 0.0)
        )
    );


    /*
     * Main calculation:
     *
     * C_lm(q) =
     *
     * 4 pi i^l
     * sum_i [
     *     g_i(q)
     *     j_l(q r_i)
     *     Y_lm*(theta_i,phi_i)
     * ]
     */
    for (size_t iq = 0; iq < qGrid.size(); ++iq)
    {
        const double q =
            qGrid[iq];

        for (const auto& atom : atoms)
        {
            /*
             * Select excluded-volume parameters.
             *
             * First try the PEPSI unified atom type.
             * If there is no PEPSI type, fall back to
             * the elemental form factor table.
             */
            const FormFactor* ff = nullptr;

            if (!atom.atomType.empty() &&
                table.contains_atom_type(atom.atomType))
            {
                ff =
                    &table.get_atom_type(
                        atom.atomType
                    );

                static int debugAtomTypeCount = 0;

                if (debugAtomTypeCount < 10)
                {
                    std::cout
                        << "Using PEPSI atom type ["
                        << atom.atomType
                        << "] excl_volume = "
                        << ff->excl_volume
                        << "\n";

                    ++debugAtomTypeCount;
                }
            }
            else if (table.contains_element(atom.Z))
            {
                ff =
                    &table.get_element(atom.Z);

                static int debugElementCount = 0;

                if (debugElementCount < 10)
                {
                    std::cout
                        << "Using elemental fallback for "
                        << atom.atomName
                        << " (Z index "
                        << atom.Z
                        << ") excl_volume = "
                        << ff->excl_volume
                        << "\n";

                    ++debugElementCount;
                }
            }
            else
            {
                std::cerr
                    << "ExcludedVolumeModel: "
                    << "no parameters for atom "
                    << atom.atomName
                    << " in residue "
                    << atom.residueName
                    << "\n";

                continue;
            }

            const double G =
                1.0 + 3.0 * radiusAdjustmentFraction;

            const double volume =
                ff->excl_volume * G;

            if (volume <= 0.0)
                continue;

            const double g =
                dummy_form_factor(
                    q,
                    volume
                );

            /*
             * Convert CENTRED Cartesian coordinates
             * into spherical coordinates.
             *
             * This is the key change: subtract the same
             * arithmetic-mean centre used by the atomic
             * multipole calculation.
             *
             * Coordinates are in Angstroms.
             */
            const double x =
                static_cast<double>(atom.x) - centreX;

            const double y =
                static_cast<double>(atom.y) - centreY;

            const double z =
                static_cast<double>(atom.z) - centreZ;

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
                const double cosTheta =
                    std::clamp(z / r, -1.0, 1.0);

                theta = std::acos(cosTheta);
                phi = std::atan2(y, x);
            }

            /*
             * q*r is dimensionless:
             *
             * q -> Angstrom^-1
             * r -> Angstrom
             */
            const double qr =
                q * r;

            /*
             * i^l:
             *
             * l % 4 = 0 ->  1
             * l % 4 = 1 ->  i
             * l % 4 = 2 -> -1
             * l % 4 = 3 -> -i
             */
            for (int l = 0; l <= L_; ++l)
            {
                const double jl =
                    spherical_bessel_j(
                        l,
                        qr
                    );

                std::complex<double> iPower;

                switch (l % 4)
                {
                    case 0:
                        iPower = {1.0, 0.0};
                        break;

                    case 1:
                        iPower = {0.0, 1.0};
                        break;

                    case 2:
                        iPower = {-1.0, 0.0};
                        break;

                    default:
                        iPower = {0.0, -1.0};
                        break;
                }

                const std::complex<double> prefactor =
                    4.0 *
                    PI *
                    iPower *
                    g *
                    jl;

                for (int m = -l; m <= l; ++m)
                {
                    const auto Y =
                        spherical_harmonic(
                            l,
                            m,
                            theta,
                            phi
                        );

                    /*
                     * Equation uses Y_lm*.
                     */
                    const auto contribution =
                        prefactor *
                        std::conj(Y);

                    coefficients_[iq]
                        [lm_index(l, m)]
                        += contribution;
                }
            }
        }
    }

    std::cout
        << "ExcludedVolumeModel:\n"
        << "  Atoms: "
        << atoms.size()
        << "\n"
        << "  q points: "
        << qGrid.size()
        << "\n"
        << "  L: "
        << L_
        << "\n"
        << "  coefficients per q: "
        << nCoefficients
        << "\n";


        // ------------------------------------------------------------
        // Diagnostic: q = 0 excluded-volume normalisation
        //
        // Expected:
        //     C_00(0) = sqrt(4*pi) * sum_i V_i
        //
        // Do not compare a nonzero-q coefficient against this value.
        // ------------------------------------------------------------

        constexpr double Q_ZERO_TOLERANCE = 1.0e-12;

        size_t qZeroIndex = qGrid.size();

        for (size_t iq = 0; iq < qGrid.size(); ++iq)
        {
            if (std::abs(qGrid[iq]) <= Q_ZERO_TOLERANCE)
            {
                qZeroIndex = iq;
                break;
            }
        }

        if (qZeroIndex < qGrid.size())
        {
            double totalVolume = 0.0;

            for (const auto& atom : atoms)
            {
                const FormFactor* ff = nullptr;

                if (!atom.atomType.empty() &&
                    table.contains_atom_type(atom.atomType))
                {
                    ff = &table.get_atom_type(atom.atomType);
                }
                else if (table.contains_element(atom.Z))
                {
                    ff = &table.get_element(atom.Z);
                }

                if (ff != nullptr && ff->excl_volume > 0.0)
                    totalVolume += ff->excl_volume;
            }

            const double expectedC00 =
                std::sqrt(4.0 * PI) * totalVolume;

            const auto actualC00 =
                coefficients_[qZeroIndex][lm_index(0, 0)];

            const double error =
                std::abs(actualC00.real() - expectedC00);

            const double relativeError =
                expectedC00 > 1.0e-12
                    ? error / expectedC00
                    : 0.0;

            std::cout
                << "\nExcluded-volume q=0 normalisation check:\n"
                << "  q used: " << qGrid[qZeroIndex] << " A^-1\n"
                << "  Total excluded volume: " << totalVolume << " A^3\n"
                << "  Expected C_00: " << expectedC00 << "\n"
                << "  Actual C_00:   " << actualC00 << "\n"
                << "  Absolute error: " << error << "\n"
                << "  Relative error: " << relativeError << "\n";
        }
        else
        {
            std::cout
                << "\nExcluded-volume q=0 normalisation check skipped: "
                << "the supplied q grid does not contain q=0.\n";
        }
}


/*
 * Access calculated coefficients.
 */
const std::vector<
    std::vector<std::complex<double>>
>&
ExcludedVolumeModel::coefficients() const
{
    return coefficients_;
}
