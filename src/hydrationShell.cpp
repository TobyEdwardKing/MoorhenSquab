
#include "hydrationShell.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <iostream>
#include <utility>

// std::vector<HydrationPoint>


namespace
{
    constexpr double PI = 3.14159265358979323846;

    double associated_legendre(int l, int m, double x)
    {
        if (m < 0 || m > l)
            return 0.0;

        x = std::max(-1.0, std::min(1.0, x));

        double pmm = 1.0;

        if (m > 0)
        {
            const double somx2 =
                std::sqrt(std::max(0.0, (1.0 - x) * (1.0 + x)));

            double fact = 1.0;

            for (int i = 1; i <= m; ++i)
            {
                pmm *= -fact * somx2;
                fact += 2.0;
            }
        }

        if (l == m)
            return pmm;

        double pmmp1 = x * (2.0 * m + 1.0) * pmm;

        if (l == m + 1)
            return pmmp1;

        for (int ll = m + 2; ll <= l; ++ll)
        {
            const double pll =
                (
                    (2.0 * ll - 1.0) * x * pmmp1
                    - (ll + m - 1.0) * pmm
                ) / (ll - m);

            pmm = pmmp1;
            pmmp1 = pll;
        }

        return pmmp1;
    }

    std::complex<double> spherical_harmonic(
        int l, int m, double theta, double phi)
    {
        if (m < 0)
        {
            const int absM = -m;

            return std::pow(-1.0, absM) *
                std::conj(
                    spherical_harmonic(l, absM, theta, phi)
                );
        }

        const double P =
            associated_legendre(l, m, std::cos(theta));

        const double normalisation = std::sqrt(
            ((2.0 * l + 1.0) / (4.0 * PI)) *
            std::tgamma(l - m + 1.0) /
            std::tgamma(l + m + 1.0)
        );

        const double angle = static_cast<double>(m) * phi;

        return normalisation * P *
            std::complex<double>(std::cos(angle), std::sin(angle));
    }

    double spherical_bessel_j(int l, double x)
    {
        if (l < 0)
            return 0.0;

        if (std::abs(x) < 1.0e-12)
            return l == 0 ? 1.0 : 0.0;

        if (std::abs(x) < 1.0)
        {
            double prefactor = 1.0;

            for (int k = 1; k <= l; ++k)
                prefactor *= x / static_cast<double>(2 * k + 1);

            double sum = 1.0;
            double term = 1.0;
            const double x2 = x * x;

            for (int k = 1; k < 100; ++k)
            {
                term *= -x2 /
                    (2.0 * k * static_cast<double>(2 * l + 2 * k + 1));

                sum += term;

                if (std::abs(term) < std::abs(sum) * 1.0e-15)
                    break;
            }

            return prefactor * sum;
        }

        const int nmax = l + 40;

        double jp1 = 0.0;
        double jn = 1.0;
        double jl = 0.0;
        double j0down = 0.0;

        for (int n = nmax; n >= 1; --n)
        {
            const double jm1 =
                (static_cast<double>(2 * n + 1) / x) * jn - jp1;

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

    int lm_index(int l, int m)
    {
        return l * l + l + m;
    }
}

std::vector<HydrationPoint> HydrationShell::generate(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table,
    double shellThickness,
    double gridSpacing
) const
{
    if (shellThickness <= 0.0)
    {
        throw std::invalid_argument(
            "Hydration shell thickness must be positive."
        );
    }

    if (gridSpacing <= 0.0)
    {
        throw std::invalid_argument(
            "Hydration grid spacing must be positive."
        );
    }

    std::vector<HydrationPoint> shell;

    if (atoms.empty())
        return shell;


    // ---------------------------------------------------------
    // 1. Build the atomic-sphere representation of the protein.
    //
    // The atomic radii define the molecular surface used for
    // this initial geometric model.
    // ---------------------------------------------------------

    std::vector<AtomSphere> spheres;
    spheres.reserve(atoms.size());

    for (const auto& atom : atoms)
    {
        if (!table.contains_element(atom.Z))
            continue;

        const auto& element =
            table.get_element(atom.Z);

        if (element.atomic_radius <= 0.0)
            continue;

        spheres.push_back({
            static_cast<double>(atom.x),
            static_cast<double>(atom.y),
            static_cast<double>(atom.z),
            element.atomic_radius
        });
    }

    if (spheres.empty())
        return shell;


    // ---------------------------------------------------------
    // 2. Determine the bounding box of the protein plus shell.
    // ---------------------------------------------------------

    double minX = spheres.front().x - spheres.front().radius;
    double minY = spheres.front().y - spheres.front().radius;
    double minZ = spheres.front().z - spheres.front().radius;

    double maxX = spheres.front().x + spheres.front().radius;
    double maxY = spheres.front().y + spheres.front().radius;
    double maxZ = spheres.front().z + spheres.front().radius;

    for (const auto& sphere : spheres)
    {
        minX = std::min(minX, sphere.x - sphere.radius);
        minY = std::min(minY, sphere.y - sphere.radius);
        minZ = std::min(minZ, sphere.z - sphere.radius);

        maxX = std::max(maxX, sphere.x + sphere.radius);
        maxY = std::max(maxY, sphere.y + sphere.radius);
        maxZ = std::max(maxZ, sphere.z + sphere.radius);
    }

    // Extend the box by the shell thickness in every direction.
    minX -= shellThickness;
    minY -= shellThickness;
    minZ -= shellThickness;

    maxX += shellThickness;
    maxY += shellThickness;
    maxZ += shellThickness;


    // ---------------------------------------------------------
    // 3. Calculate grid dimensions.
    //
    // Each accepted grid cell represents gridSpacing^3 A^3.
    // ---------------------------------------------------------

    const auto gridCount = [gridSpacing](
        double minimum,
        double maximum
    ) -> std::size_t
    {
        return static_cast<std::size_t>(
            std::ceil((maximum - minimum) / gridSpacing)
        );
    };

    const std::size_t nx = gridCount(minX, maxX);
    const std::size_t ny = gridCount(minY, maxY);
    const std::size_t nz = gridCount(minZ, maxZ);

    const double cellVolume =
        gridSpacing * gridSpacing * gridSpacing;


    // ---------------------------------------------------------
    // 4. Classify each grid-cell centre.
    //
    // For each point, calculate its signed distance from the
    // union of atomic spheres:
    //
    // signedDistance = min(distance_to_atom - atom_radius)
    //
    // signedDistance < 0: inside the protein
    // signedDistance = 0: on the approximate surface
    // signedDistance > 0: outside the protein
    //
    // Keep points outside the protein and no farther than
    // shellThickness from its approximate surface.
    // ---------------------------------------------------------

    for (std::size_t ix = 0; ix < nx; ++ix)
    {
        const double x =
            minX + (static_cast<double>(ix) + 0.5)
                   * gridSpacing;

        for (std::size_t iy = 0; iy < ny; ++iy)
        {
            const double y =
                minY + (static_cast<double>(iy) + 0.5)
                       * gridSpacing;

            for (std::size_t iz = 0; iz < nz; ++iz)
            {
                const double z =
                    minZ + (static_cast<double>(iz) + 0.5)
                           * gridSpacing;

                double minimumSurfaceDistance =
                    INFINITY;

                for (const auto& sphere : spheres)
                {
                    const double dx = x - sphere.x;
                    const double dy = y - sphere.y;
                    const double dz = z - sphere.z;

                    const double distance =
                        std::sqrt(
                            dx * dx +
                            dy * dy +
                            dz * dz
                        );

                    const double surfaceDistance =
                        distance - sphere.radius;

                    minimumSurfaceDistance =
                        std::min(
                            minimumSurfaceDistance,
                            surfaceDistance
                        );

                    // This point is inside at least one sphere.
                    // No other sphere can make it external.
                    if (minimumSurfaceDistance < 0.0)
                        break;
                }

                if (
                    minimumSurfaceDistance >= 0.0 &&
                    minimumSurfaceDistance <= shellThickness
                )
                {
                    shell.push_back({
                        x,
                        y,
                        z,
                        cellVolume
                    });
                }
            }
        }
    }

    return shell;
}


// -------------------------------------------------------------
// Sum the volume represented by all accepted grid cells.
// -------------------------------------------------------------

double HydrationShell::total_volume(
    const std::vector<HydrationPoint>& shell
)
{
    double volume = 0.0;

    for (const auto& point : shell)
        volume += point.volume;

    return volume;
}


void HydrationShell::compute_multipoles(
    const std::vector<HydrationPoint>& shellPoints,
    const std::vector<double>& qGrid,
    const std::vector<double>& centre,
    int L)
{
    if (L < 0)
    {
        throw std::invalid_argument(
            "Hydration multipole L must be non-negative."
        );
    }

    if (centre.size() != 3)
    {
        throw std::invalid_argument(
            "Hydration multipole centre must contain x, y, z."
        );
    }

    if (std::any_of(qGrid.begin(), qGrid.end(),
                    [](double q) { return !std::isfinite(q) || q < 0.0; }))
    {
        throw std::invalid_argument(
            "Hydration q values must be finite and non-negative."
        );
    }

    L_ = L;

    const int nCoefficients = (L + 1) * (L + 1);

    multipoles_.assign(
        qGrid.size(),
        std::vector<std::complex<double>>(
            nCoefficients,
            std::complex<double>(0.0, 0.0)
        )
    );

    const double cx = centre[0];
    const double cy = centre[1];
    const double cz = centre[2];

    // Cache each point's radius and angular harmonics. These do
    // not depend on q, so they need to be calculated only once.
    struct PointData
    {
        double r;
        double volume;
        std::vector<std::complex<double>> conjugateY;
    };

    std::vector<PointData> points;
    points.reserve(shellPoints.size());

    for (const auto& point : shellPoints)
    {
        const double x = point.x - cx;
        const double y = point.y - cy;
        const double z = point.z - cz;

        const double r = std::sqrt(x * x + y * y + z * z);

        double theta = 0.0;
        double phi = 0.0;

        if (r > 1.0e-12)
        {
            theta = std::acos(
                std::max(-1.0, std::min(1.0, z / r))
            );
            phi = std::atan2(y, x);
        }

        PointData data;
        data.r = r;
        data.volume = point.volume;
        data.conjugateY.resize(nCoefficients);

        for (int l = 0; l <= L; ++l)
        {
            for (int m = -l; m <= l; ++m)
            {
                data.conjugateY[lm_index(l, m)] =
                    std::conj(spherical_harmonic(l, m, theta, phi));
            }
        }

        points.push_back(std::move(data));
    }

    // B_lm(q) = 4 pi i^l sum_p V_p j_l(q r_p) conj(Y_lm).
    for (std::size_t iq = 0; iq < qGrid.size(); ++iq)
    {
        const double q = qGrid[iq];
        auto& coefficients = multipoles_[iq];

        for (const auto& point : points)
        {
            for (int l = 0; l <= L; ++l)
            {
                const double jl = spherical_bessel_j(l, q * point.r);

                const std::complex<double> phase =
                    std::pow(std::complex<double>(0.0, 1.0), l);

                const std::complex<double> radial =
                    4.0 * PI * phase * point.volume * jl;

                for (int m = -l; m <= l; ++m)
                {
                    const int index = lm_index(l, m);

                    coefficients[index] +=
                        radial * point.conjugateY[index];
                }
            }
        }
    }

    std::cout
        << "Hydration-shell multipoles calculated:\n"
        << "  Shell points: " << shellPoints.size() << "\n"
        << "  q points: " << qGrid.size() << "\n"
        << "  L: " << L << "\n"
        << "  Coefficients per q: " << nCoefficients << "\n"
        << "  Centre: (" << cx << ", " << cy << ", " << cz << ")\n";
}