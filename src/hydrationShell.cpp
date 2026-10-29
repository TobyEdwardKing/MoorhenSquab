
#include "hydrationShell.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include <iostream>
#include <utility>
#include <limits>
#include <map>
#include <tuple>

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

    void spherical_bessel_all(
        int L,
        double x,
        std::vector<double>& values
    )
    {
        values.resize(L + 1);

        // Preserve the existing small-x treatment.
        if (std::abs(x) < 1.0)
        {
            for (int l = 0; l <= L; ++l)
                values[l] = spherical_bessel_j(l, x);

            return;
        }

        // One downward recurrence calculates every order up to L.
        const int nmax = L + 40;

        double jp1 = 0.0;
        double jn = 1.0;

        for (int n = nmax; n >= 1; --n)
        {
            const double jm1 =
                (static_cast<double>(2 * n + 1) / x) * jn - jp1;

            jp1 = jn;
            jn = jm1;

            if (n - 1 <= L)
                values[n - 1] = jm1;
        }

        // Normalise against whichever of exact j0 or j1 has
        // the larger magnitude, avoiding division by a value
        // close to a zero of either function.
        const double exactJ0 = std::sin(x) / x;
        const double exactJ1 =
            std::sin(x) / (x * x) - std::cos(x) / x;

        double normalisation;

        if (L == 0 || std::abs(exactJ0) >= std::abs(exactJ1))
        {
            normalisation = exactJ0 / values[0];
        }
        else
        {
            const double j1Downward =
                values.size() > 1 ? values[1] : 0.0;

            normalisation = exactJ1 / j1Downward;
        }

        for (double& value : values)
            value *= normalisation;
    }

    int lm_index(int l, int m)
    {
        return l * l + l + m;
    }

    // Natural cubic-spline second derivatives for equally or unequally
    // spaced x values. x must be strictly increasing.
    std::vector<double> spline_second_derivatives(
        const std::vector<double>& x,
        const std::vector<double>& y)
    {
        const std::size_t n = x.size();
        std::vector<double> second(n, 0.0);
        if (n < 3)
            return second;

        std::vector<double> work(n - 1, 0.0);
        for (std::size_t i = 1; i + 1 < n; ++i)
        {
            const double span = x[i + 1] - x[i - 1];
            const double left = x[i] - x[i - 1];
            const double right = x[i + 1] - x[i];
            if (left <= 0.0 || right <= 0.0)
                throw std::invalid_argument("Spline q samples must be strictly increasing.");

            const double sig = left / span;
            const double pivot = sig * second[i - 1] + 2.0;
            second[i] = (sig - 1.0) / pivot;
            const double slopeDifference =
                6.0 * ((y[i + 1] - y[i]) / right -
                       (y[i] - y[i - 1]) / left) / span;
            work[i] = (slopeDifference - sig * work[i - 1]) / pivot;
        }

        for (std::size_t k = n - 1; k-- > 0;)
            second[k] = second[k] * second[k + 1] + work[k];

        return second;
    }

    double spline_interpolate(
        const std::vector<double>& x,
        const std::vector<double>& y,
        const std::vector<double>& second,
        double target)
    {
        const std::size_t n = x.size();
        if (n == 1)
            return y[0];
        if (target <= x.front())
            return y.front();
        if (target >= x.back())
            return y.back();

        const auto upper = std::upper_bound(x.begin(), x.end(), target);
        const std::size_t hi = static_cast<std::size_t>(upper - x.begin());
        const std::size_t lo = hi - 1;
        const double h = x[hi] - x[lo];
        const double a = (x[hi] - target) / h;
        const double b = (target - x[lo]) / h;
        return a * y[lo] + b * y[hi] +
            ((a * a * a - a) * second[lo] +
             (b * b * b - b) * second[hi]) * h * h / 6.0;
    }
}

    void run_spherical_bessel_test()
    {
        constexpr int L = 12;
        constexpr double tolerance = 1.0e-10;

        const std::vector<double> testX = {
            1.0e-6,
            0.1,
            0.5,
            0.999,
            1.001,
            2.0,
            2.404825557695773, // First zero of j0
            2.404825567695773, // Just above that zero
            3.0,
            3.831705970207512, // First zero of j1
            5.0,
            10.0,
            19.0
        };

        std::vector<double> values;

        double maxDifference = 0.0;
        double maxKnownError = 0.0;
        double worstX = 0.0;
        int worstL = 0;
        bool passed = true;

        for (double x : testX)
        {
            spherical_bessel_all(L, x, values);

            // Compare the optimised implementation with the original.
            for (int l = 0; l <= L; ++l)
            {
                const double reference = spherical_bessel_j(l, x);
                const double difference = std::abs(values[l] - reference);

                if (difference > maxDifference)
                {
                    maxDifference = difference;
                    worstX = x;
                    worstL = l;
                }

                if (difference > tolerance)
                    passed = false;
            }

            // Independently check j0 and j1 against their exact formulas.
            const double exactJ0 = std::sin(x) / x;
            const double exactJ1 =
                std::sin(x) / (x * x) - std::cos(x) / x;

            const double errorJ0 = std::abs(values[0] - exactJ0);
            const double errorJ1 = std::abs(values[1] - exactJ1);

            maxKnownError = std::max(
                maxKnownError, std::max(errorJ0, errorJ1));

            if (errorJ0 > tolerance || errorJ1 > tolerance)
                passed = false;
        }

        std::cout << "\n=== SPHERICAL BESSEL TEST ===\n";
        std::cout << "L: " << L << '\n';
        std::cout << "Test points: " << testX.size() << '\n';
        std::cout << "Max difference vs original: "
                << maxDifference << '\n';
        std::cout << "Worst comparison: l=" << worstL
                << ", x=" << worstX << '\n';
        std::cout << "Max error vs exact j0/j1: "
                << maxKnownError << '\n';
        std::cout << "Result: "
                << (passed ? "PASS" : "FAIL") << '\n';
    }

std::vector<HydrationPoint> HydrationShell::generate(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& /*table*/,
    double gridSpacing
) const
{
    if (!std::isfinite(gridSpacing) || gridSpacing <= 0.0)
        throw std::invalid_argument("Hydration grid spacing must be positive and finite.");

    std::vector<HydrationPoint> shell;
    if (atoms.empty())
        return shell;

    // Pepsi-SAXS-style shell width based on the coordinate-derived
    // radius of gyration of the molecule.
    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (const auto& atom : atoms) {
        cx += atom.x;
        cy += atom.y;
        cz += atom.z;
    }
    cx /= static_cast<double>(atoms.size());
    cy /= static_cast<double>(atoms.size());
    cz /= static_cast<double>(atoms.size());

    double rg2 = 0.0;
    for (const auto& atom : atoms) {
        const double dx = atom.x - cx;
        const double dy = atom.y - cy;
        const double dz = atom.z - cz;
        rg2 += dx * dx + dy * dy + dz * dz;
    }
    const double molecularRg = std::sqrt(rg2 / static_cast<double>(atoms.size()));

    double shellThickness = 3.0;
    if (molecularRg >= 20.0)
        shellThickness = 5.0;
    else if (molecularRg > 15.0)
        shellThickness = 3.0 + (molecularRg - 15.0) * (2.0 / 5.0);

    lastShellThickness_ = shellThickness;
    lastGridSpacing_ = gridSpacing;

    // The paper pads the grid by at least 12 A in each direction.
    constexpr double PADDING = 12.0;
    constexpr double INNER_DISTANCE = 3.0;

    double minX = atoms.front().x, minY = atoms.front().y, minZ = atoms.front().z;
    double maxX = minX, maxY = minY, maxZ = minZ;
    for (const auto& atom : atoms) {
        minX = std::min(minX, static_cast<double>(atom.x));
        minY = std::min(minY, static_cast<double>(atom.y));
        minZ = std::min(minZ, static_cast<double>(atom.z));
        maxX = std::max(maxX, static_cast<double>(atom.x));
        maxY = std::max(maxY, static_cast<double>(atom.y));
        maxZ = std::max(maxZ, static_cast<double>(atom.z));
    }
    minX -= PADDING; minY -= PADDING; minZ -= PADDING;
    maxX += PADDING; maxY += PADDING; maxZ += PADDING;

    const auto gridCount = [gridSpacing](double lo, double hi) {
        return static_cast<std::size_t>(std::ceil((hi - lo) / gridSpacing));
    };
    const std::size_t nx = gridCount(minX, maxX);
    const std::size_t ny = gridCount(minY, maxY);
    const std::size_t nz = gridCount(minZ, maxZ);
    const double cellVolume = gridSpacing * gridSpacing * gridSpacing;

    // Assign each atom to a grid cell. Candidate cell centres only need
    // to inspect atoms assigned to that cell and its 26 neighbours.
    using Cell = std::tuple<int, int, int>;
    std::map<Cell, std::vector<std::size_t>> atomsByCell;
    for (std::size_t i = 0; i < atoms.size(); ++i) {
        const int ix = static_cast<int>(std::floor((atoms[i].x - minX) / gridSpacing));
        const int iy = static_cast<int>(std::floor((atoms[i].y - minY) / gridSpacing));
        const int iz = static_cast<int>(std::floor((atoms[i].z - minZ) / gridSpacing));
        atomsByCell[{ix, iy, iz}].push_back(i);
    }

    const double outerDistance = INNER_DISTANCE + shellThickness;
    for (std::size_t ix = 0; ix < nx; ++ix) {
        const double x = minX + (static_cast<double>(ix) + 0.5) * gridSpacing;
        for (std::size_t iy = 0; iy < ny; ++iy) {
            const double y = minY + (static_cast<double>(iy) + 0.5) * gridSpacing;
            for (std::size_t iz = 0; iz < nz; ++iz) {
                const double z = minZ + (static_cast<double>(iz) + 0.5) * gridSpacing;
                const int cellX = static_cast<int>(ix);
                const int cellY = static_cast<int>(iy);
                const int cellZ = static_cast<int>(iz);

                double minDistanceSquared = std::numeric_limits<double>::infinity();
                for (int dxCell = -1; dxCell <= 1; ++dxCell) {
                    for (int dyCell = -1; dyCell <= 1; ++dyCell) {
                        for (int dzCell = -1; dzCell <= 1; ++dzCell) {
                            const auto found = atomsByCell.find(
                                {cellX + dxCell, cellY + dyCell, cellZ + dzCell});
                            if (found == atomsByCell.end())
                                continue;

                            for (const std::size_t atomIndex : found->second) {
                                const double dx = x - atoms[atomIndex].x;
                                const double dy = y - atoms[atomIndex].y;
                                const double dz = z - atoms[atomIndex].z;
                                const double d2 = dx * dx + dy * dy + dz * dz;
                                minDistanceSquared = std::min(minDistanceSquared, d2);
                            }
                        }
                    }
                }

                if (!std::isfinite(minDistanceSquared))
                    continue;

                const double minDistance = std::sqrt(minDistanceSquared);
                if (minDistance >= INNER_DISTANCE && minDistance <= outerDistance) {
                    shell.push_back({x, y, z, cellVolume});
                }
            }
        }
    }

    std::cout << "Hydration shell generated:\n"
              << "  Coordinate-derived molecular Rg: " << molecularRg << " A\n"
              << "  Shell thickness: " << shellThickness << " A\n"
              << "  Grid-cell width: " << gridSpacing << " A\n"
              << "  Padding: " << PADDING << " A\n"
              << "  Shell points: " << shell.size() << "\n";

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

double HydrationShell::radius_of_gyration(
    const std::vector<HydrationPoint>& points
)
{
    if (points.empty())
        return 0.0;

    double cx = 0.0, cy = 0.0, cz = 0.0;
    for (const auto& point : points) {
        cx += point.x;
        cy += point.y;
        cz += point.z;
    }
    const double count = static_cast<double>(points.size());
    cx /= count; cy /= count; cz /= count;

    double sumSquaredDistance = 0.0;
    for (const auto& point : points) {
        const double dx = point.x - cx;
        const double dy = point.y - cy;
        const double dz = point.z - cz;
        sumSquaredDistance += dx * dx + dy * dy + dz * dz;
    }
    return std::sqrt(sumSquaredDistance / count);
}



void HydrationShell::compute_multipoles(
    const std::vector<HydrationPoint>& shellPoints,
    const FormFactorTable& table,
    const std::vector<double>& qGrid,
    const std::vector<double>& centre,
    int L,
    bool useInterpolation)
{
    if (L < 0)
        throw std::invalid_argument("Hydration multipole L must be non-negative.");
    if (centre.size() != 3)
        throw std::invalid_argument("Hydration multipole centre must contain x, y, z.");
    if (qGrid.empty())
        throw std::invalid_argument("Hydration q grid must not be empty.");
    if (std::any_of(qGrid.begin(), qGrid.end(), [](double q) {
            return !std::isfinite(q) || q < 0.0;
        }))
        throw std::invalid_argument("Hydration q values must be finite and non-negative.");
    if (!std::is_sorted(qGrid.begin(), qGrid.end()) ||
        std::adjacent_find(qGrid.begin(), qGrid.end(), std::greater_equal<double>()) != qGrid.end())
        throw std::invalid_argument("Hydration q values must be strictly increasing.");

    L_ = L;
    const std::size_t coefficientCount = static_cast<std::size_t>((L + 1) * (L + 1));
    const std::size_t nPoints = shellPoints.size();

    const double qMin = qGrid.front();
    const double qMax = qGrid.back();
    const bool hasQRange = qMax > qMin;

    const std::size_t sampleCount = useInterpolation
        ? (hasQRange
            ? std::max<std::size_t>(
                2, static_cast<std::size_t>(2 * L))
            : 1)
        : qGrid.size();

    std::vector<double> sampleQ(sampleCount, qMin);

    if (useInterpolation && hasQRange) {
        for (std::size_t i = 0; i < sampleCount; ++i) {
            const double t =
                static_cast<double>(i) /
                static_cast<double>(sampleCount - 1);

            sampleQ[i] = qMin + t * (qMax - qMin);
        }
    } else if (!useInterpolation) {
        sampleQ = qGrid;
    }

    const double cx = centre[0], cy = centre[1], cz = centre[2];
    std::vector<double> radii(nPoints), volumes(nPoints);
    std::vector<std::complex<double>> conjugateY(nPoints * coefficientCount);

    for (std::size_t p = 0; p < nPoints; ++p) {
        const auto& point = shellPoints[p];
        const double x = point.x - cx, y = point.y - cy, z = point.z - cz;
        const double r = std::sqrt(x * x + y * y + z * z);
        radii[p] = r;
        volumes[p] = point.volume;

        double theta = 0.0, phi = 0.0;
        if (r > 1.0e-12) {
            theta = std::acos(std::max(-1.0, std::min(1.0, z / r)));
            phi = std::atan2(y, x);
        }
        const std::size_t base = p * coefficientCount;
        for (int l = 0; l <= L; ++l) {
            for (int m = -l; m <= l; ++m) {
                const std::size_t index = static_cast<std::size_t>(lm_index(l, m));
                conjugateY[base + index] =
                    std::conj(spherical_harmonic(l, m, theta, phi));
            }
        }
    }

    std::vector<std::complex<double>> phases(static_cast<std::size_t>(L + 1));
    const std::complex<double> imaginaryUnit(0.0, 1.0);
    for (int l = 0; l <= L; ++l)
        phases[static_cast<std::size_t>(l)] = std::pow(imaginaryUnit, l);

    std::vector<std::vector<std::complex<double>>> sampledCoefficients(
        sampleCount,
        std::vector<std::complex<double>>(coefficientCount, {0.0, 0.0}));

    std::vector<double> besselValues(static_cast<std::size_t>(L + 1));
    for (std::size_t iq = 0; iq < sampleCount; ++iq) {
        const double q = sampleQ[iq];
        auto& coefficients = sampledCoefficients[iq];

        for (std::size_t p = 0; p < nPoints; ++p) {
            spherical_bessel_all(L, q * radii[p], besselValues);
            const std::size_t base = p * coefficientCount;

            for (int l = 0; l <= L; ++l) {
                const std::complex<double> radial =
                    4.0 * PI * phases[static_cast<std::size_t>(l)] *
                    volumes[p] * besselValues[static_cast<std::size_t>(l)];

                for (int m = -l; m <= l; ++m) {
                    const std::size_t index = static_cast<std::size_t>(lm_index(l, m));
                    coefficients[index] += radial * conjugateY[base + index];
                }
            }
        }
    }

    multipoles_.assign(
        qGrid.size(),
        std::vector<std::complex<double>>(
            coefficientCount, {0.0, 0.0}));

    if (!useInterpolation) {
        // Every requested q-value was calculated directly.
        multipoles_ = std::move(sampledCoefficients);
    } else if (!hasQRange) {
        // Single-q normalisation diagnostic.
        multipoles_[0] = sampledCoefficients[0];
    } else {
    // Natural cubic spline interpolation is applied independently to the
    // real and imaginary parts of each complex coefficient.
    const std::size_t n = sampleCount;
    std::vector<double> secondDerivative(n, 0.0);
    std::vector<double> u(n, 0.0);

    for (std::size_t coefficient = 0; coefficient < coefficientCount; ++coefficient) {
        // Build the natural-spline second derivatives for the real component.
        secondDerivative.assign(n, 0.0);
        u.assign(n, 0.0);
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const double sig = (sampleQ[i] - sampleQ[i - 1]) /
                               (sampleQ[i + 1] - sampleQ[i - 1]);
            const double p = sig * secondDerivative[i - 1] + 2.0;
            const double yPrev = sampledCoefficients[i - 1][coefficient].real();
            const double yHere = sampledCoefficients[i][coefficient].real();
            const double yNext = sampledCoefficients[i + 1][coefficient].real();
            const double slopeDelta =
                (yNext - yHere) / (sampleQ[i + 1] - sampleQ[i]) -
                (yHere - yPrev) / (sampleQ[i] - sampleQ[i - 1]);
            secondDerivative[i] = (sig - 1.0) / p;
            u[i] = (6.0 * slopeDelta / (sampleQ[i + 1] - sampleQ[i - 1]) -
                    sig * u[i - 1]) / p;
        }
        for (std::size_t k = n - 1; k-- > 0;)
            secondDerivative[k] = secondDerivative[k] * secondDerivative[k + 1] + u[k];

        std::vector<double> realSecond = secondDerivative;

        // Build the natural-spline second derivatives for the imaginary component.
        secondDerivative.assign(n, 0.0);
        u.assign(n, 0.0);
        for (std::size_t i = 1; i + 1 < n; ++i) {
            const double sig = (sampleQ[i] - sampleQ[i - 1]) /
                               (sampleQ[i + 1] - sampleQ[i - 1]);
            const double p = sig * secondDerivative[i - 1] + 2.0;
            const double yPrev = sampledCoefficients[i - 1][coefficient].imag();
            const double yHere = sampledCoefficients[i][coefficient].imag();
            const double yNext = sampledCoefficients[i + 1][coefficient].imag();
            const double slopeDelta =
                (yNext - yHere) / (sampleQ[i + 1] - sampleQ[i]) -
                (yHere - yPrev) / (sampleQ[i] - sampleQ[i - 1]);
            secondDerivative[i] = (sig - 1.0) / p;
            u[i] = (6.0 * slopeDelta / (sampleQ[i + 1] - sampleQ[i - 1]) -
                    sig * u[i - 1]) / p;
        }
        for (std::size_t k = n - 1; k-- > 0;)
            secondDerivative[k] = secondDerivative[k] * secondDerivative[k + 1] + u[k];
        std::vector<double> imagSecond = secondDerivative;

        for (std::size_t iq = 0; iq < qGrid.size(); ++iq) {
            const double q = qGrid[iq];
            auto upper = std::upper_bound(sampleQ.begin(), sampleQ.end(), q);
            std::size_t hi = upper == sampleQ.end() ? n - 1 :
                             static_cast<std::size_t>(upper - sampleQ.begin());
            if (hi == 0) hi = 1;
            const std::size_t lo = hi - 1;
            const double h = sampleQ[hi] - sampleQ[lo];
            const double a = (sampleQ[hi] - q) / h;
            const double b = (q - sampleQ[lo]) / h;

            const double realValue =
                a * sampledCoefficients[lo][coefficient].real() +
                b * sampledCoefficients[hi][coefficient].real() +
                ((a * a * a - a) * realSecond[lo] +
                 (b * b * b - b) * realSecond[hi]) * h * h / 6.0;
            const double imagValue =
                a * sampledCoefficients[lo][coefficient].imag() +
                b * sampledCoefficients[hi][coefficient].imag() +
                ((a * a * a - a) * imagSecond[lo] +
                 (b * b * b - b) * imagSecond[hi]) * h * h / 6.0;

            multipoles_[iq][coefficient] = {realValue, imagValue};
        }
    }
    }
    
    // CrabSAXS / Pepsi-SAXS united-group approximation for water.
    const double fO = table.get_element(7).evaluate(0.0);
    const double fH = table.get_element(0).evaluate(0.0);

    // Store the water form factor at every requested q value.
    const auto& oxygen = table.get_element(7);
    const auto& hydrogen = table.get_element(0);

    waterFormFactors_.resize(qGrid.size());

    constexpr double OH_DISTANCE = 0.9572;

    for (std::size_t iq = 0; iq < qGrid.size(); ++iq)
    {
        const double q = qGrid[iq];
        const double x = q * OH_DISTANCE;

        const double sinc = std::abs(x) < 1.0e-12
            ? 1.0
            : std::sin(x) / x;

        waterFormFactors_[iq] =
            oxygen.evaluate(q) + 2.0 * hydrogen.evaluate(q) * sinc;

        if (iq == 0 || iq == qGrid.size() / 2 ||
            iq + 1 == qGrid.size())
        {
            std::cout
                << "Water form factor: q = " << q
                << ", h(q) = " << waterFormFactors_[iq]
                << '\n';
        }
    }
    std::cout << "Hydration-shell multipoles calculated with q interpolation:\n"
              << "  Shell points: " << nPoints << "\n"
              << "  Requested q points: " << qGrid.size() << "\n"
              << "  Sampled q points: " << sampleCount << "\n"
              << "  L: " << L << "\n"
              << "  Coefficients per q: " << coefficientCount << "\n"
              << "  Centre: (" << cx << ", " << cy << ", " << cz << ")\n";
}

void HydrationShell::validate_interpolation(
    const std::vector<HydrationPoint>& shellPoints,
    const std::vector<double>& qGrid,
    const std::vector<double>& centre,
    int maxTests
) const
{
    if (L_ < 0)
        throw std::runtime_error(
            "Run compute_multipoles before validating interpolation.");

    if (centre.size() != 3)
        throw std::invalid_argument(
            "Interpolation validation centre must contain x, y, z.");

    if (qGrid.size() != multipoles_.size())
        throw std::invalid_argument(
            "Interpolation validation q grid does not match calculated coefficients.");

    if (qGrid.size() < 2 || maxTests <= 0 || shellPoints.empty())
        throw std::invalid_argument(
            "Interpolation validation requires shell points, q values and a positive test count.");

    const std::size_t coefficientCount =
        static_cast<std::size_t>((L_ + 1) * (L_ + 1));

    // This must match the sparse sample count in compute_multipoles().
    const std::size_t sampleCount =
        std::max<std::size_t>(2, static_cast<std::size_t>(2 * L_));

    const double qMin = qGrid.front();
    const double qMax = qGrid.back();

    if (!(qMax > qMin))
        throw std::invalid_argument(
            "Interpolation validation requires a non-zero q range.");

    std::vector<double> sampleQ(sampleCount);

    for (std::size_t i = 0; i < sampleCount; ++i)
    {
        const double t =
            static_cast<double>(i) /
            static_cast<double>(sampleCount - 1);

        sampleQ[i] = qMin + t * (qMax - qMin);
    }

    // Pick experimental-grid indices nearest the midpoints between
    // adjacent sparse samples.
    std::vector<std::size_t> testIndices;
    testIndices.reserve(
        std::min<std::size_t>(
            static_cast<std::size_t>(maxTests),
            sampleCount - 1));

    for (std::size_t i = 0;
         i + 1 < sampleCount &&
         testIndices.size() < static_cast<std::size_t>(maxTests);
         ++i)
    {
        const double midpoint = 0.5 * (sampleQ[i] + sampleQ[i + 1]);

        const auto upper =
            std::lower_bound(qGrid.begin(), qGrid.end(), midpoint);

        std::size_t index = 0;

        if (upper == qGrid.end())
        {
            index = qGrid.size() - 1;
        }
        else if (upper == qGrid.begin())
        {
            index = 0;
        }
        else
        {
            const std::size_t hi =
                static_cast<std::size_t>(upper - qGrid.begin());
            const std::size_t lo = hi - 1;

            index = std::abs(qGrid[lo] - midpoint) <=
                    std::abs(qGrid[hi] - midpoint)
                ? lo
                : hi;
        }

        if (testIndices.empty() || testIndices.back() != index)
            testIndices.push_back(index);
    }

    const double cx = centre[0];
    const double cy = centre[1];
    const double cz = centre[2];

    // Precompute the geometry-dependent terms exactly as in
    // compute_multipoles().
    const std::size_t nPoints = shellPoints.size();

    std::vector<double> radii(nPoints);
    std::vector<double> volumes(nPoints);

    std::vector<std::complex<double>> conjugateY(
        nPoints * coefficientCount);

    for (std::size_t p = 0; p < nPoints; ++p)
    {
        const auto& point = shellPoints[p];

        const double x = point.x - cx;
        const double y = point.y - cy;
        const double z = point.z - cz;

        const double r = std::sqrt(x * x + y * y + z * z);

        radii[p] = r;
        volumes[p] = point.volume;

        double theta = 0.0;
        double phi = 0.0;

        if (r > 1.0e-12)
        {
            theta = std::acos(
                std::max(-1.0, std::min(1.0, z / r)));
            phi = std::atan2(y, x);
        }

        const std::size_t base = p * coefficientCount;

        for (int l = 0; l <= L_; ++l)
        {
            for (int m = -l; m <= l; ++m)
            {
                const std::size_t index =
                    static_cast<std::size_t>(lm_index(l, m));

                conjugateY[base + index] =
                    std::conj(spherical_harmonic(l, m, theta, phi));
            }
        }
    }

    std::vector<std::complex<double>> phases(
        static_cast<std::size_t>(L_ + 1));

    const std::complex<double> imaginaryUnit(0.0, 1.0);

    for (int l = 0; l <= L_; ++l)
        phases[static_cast<std::size_t>(l)] =
            std::pow(imaginaryUnit, l);

    std::vector<double> besselValues(
        static_cast<std::size_t>(L_ + 1));

    double maximumNormalisedError = 0.0;
    double sumSquaredNormalisedErrors = 0.0;
    std::size_t coefficientComparisons = 0;

    std::cout << "\n=== HYDRATION INTERPOLATION VALIDATION ===\n"
              << "L: " << L_ << '\n'
              << "Sparse q samples: " << sampleCount << '\n'
              << "Shell points: " << nPoints << '\n'
              << "Test q values: " << testIndices.size() << '\n'
              << "q | max coefficient error / max direct coefficient"
              << " | RMS coefficient error / max direct coefficient\n";

    for (const std::size_t iq : testIndices)
    {
        const double q = qGrid[iq];

        std::vector<std::complex<double>> direct(
            coefficientCount, {0.0, 0.0});

        // Direct calculation: sum over every shell point at this q.
        // This does not use the sparse samples or spline interpolation.
        for (std::size_t p = 0; p < nPoints; ++p)
        {
            spherical_bessel_all(
                L_, q * radii[p], besselValues);

            const std::size_t base = p * coefficientCount;

            for (int l = 0; l <= L_; ++l)
            {
                const std::complex<double> radial =
                    4.0 * PI *
                    phases[static_cast<std::size_t>(l)] *
                    volumes[p] *
                    besselValues[static_cast<std::size_t>(l)];

                for (int m = -l; m <= l; ++m)
                {
                    const std::size_t index =
                        static_cast<std::size_t>(lm_index(l, m));

                    direct[index] += radial * conjugateY[base + index];
                }
            }
        }

        double referenceMagnitude = 0.0;

        for (const auto& value : direct)
            referenceMagnitude =
                std::max(referenceMagnitude, std::abs(value));

        if (referenceMagnitude <= 1.0e-30)
        {
            std::cout << q << " | skipped (direct coefficients are near zero)\n";
            continue;
        }

        double maxErrorAtQ = 0.0;
        double sumSquaredErrorAtQ = 0.0;

        for (std::size_t c = 0; c < coefficientCount; ++c)
        {
            const double error =
                std::abs(multipoles_[iq][c] - direct[c]);

            const double normalisedError =
                error / referenceMagnitude;

            maxErrorAtQ =
                std::max(maxErrorAtQ, normalisedError);

            sumSquaredErrorAtQ +=
                normalisedError * normalisedError;
        }

        const double rmsErrorAtQ = std::sqrt(
            sumSquaredErrorAtQ /
            static_cast<double>(coefficientCount));

        maximumNormalisedError =
            std::max(maximumNormalisedError, maxErrorAtQ);

        sumSquaredNormalisedErrors += sumSquaredErrorAtQ;
        coefficientComparisons += coefficientCount;

        std::cout << q << " | "
                  << maxErrorAtQ << " | "
                  << rmsErrorAtQ << '\n';
    }

    std::cout << "Maximum normalised coefficient error: "
              << maximumNormalisedError << '\n';

    if (coefficientComparisons > 0)
    {
        std::cout << "Overall RMS normalised coefficient error: "
                  << std::sqrt(
                         sumSquaredNormalisedErrors /
                         static_cast<double>(coefficientComparisons))
                  << '\n';
    }

    std::cout << "==========================================\n";
}


void HydrationShell::validate_combined_intensity_interpolation(
    const std::vector<HydrationPoint>& shellPoints,
    const std::vector<double>& qGrid,
    const std::vector<double>& centre,
    const std::vector<std::vector<std::complex<double>>>& atomic,
    const std::vector<std::vector<std::complex<double>>>& excluded,
    double solventDensity,
    double shellContrast,
    double radiusAdjustmentFraction,
    int maxTests
) const
{
    if (L_ < 0)
        throw std::runtime_error(
            "Run compute_multipoles before validating combined intensity.");

    if (!std::isfinite(radiusAdjustmentFraction) ||
        radiusAdjustmentFraction < -0.05 ||
        radiusAdjustmentFraction > 0.05)
        throw std::invalid_argument(
            "Radius adjustment fraction must be finite and within [-0.05, 0.05].");

    if (centre.size() != 3)
        throw std::invalid_argument(
            "Combined-intensity validation centre must contain x, y, z.");

    if (qGrid.size() != multipoles_.size() ||
        qGrid.size() != atomic.size() ||
        qGrid.size() != excluded.size())
    {
        throw std::invalid_argument(
            "Combined-intensity validation arrays have different q-grid sizes.");
    }

    if (qGrid.size() < 2 || maxTests <= 0 || shellPoints.empty())
        throw std::invalid_argument(
            "Combined-intensity validation requires shell points, q values "
            "and a positive test count.");

    const std::size_t coefficientCount =
        static_cast<std::size_t>((L_ + 1) * (L_ + 1));

    for (std::size_t iq = 0; iq < qGrid.size(); ++iq)
    {
        if (atomic[iq].size() != coefficientCount ||
            excluded[iq].size() != coefficientCount ||
            multipoles_[iq].size() != coefficientCount)
        {
            throw std::invalid_argument(
                "Combined-intensity validation coefficient counts do not match L.");
        }
    }

    const std::size_t sampleCount =
        std::max<std::size_t>(2, static_cast<std::size_t>(2 * L_));

    const double qMin = qGrid.front();
    const double qMax = qGrid.back();

    if (!(qMax > qMin))
        throw std::invalid_argument(
            "Combined-intensity validation requires a non-zero q range.");

    std::vector<double> sampleQ(sampleCount);

    for (std::size_t i = 0; i < sampleCount; ++i)
    {
        const double t = static_cast<double>(i) /
                         static_cast<double>(sampleCount - 1);

        sampleQ[i] = qMin + t * (qMax - qMin);
    }

    // Select experimental-grid indices nearest the sparse-sample midpoints.
    std::vector<std::size_t> testIndices;
    testIndices.reserve(
        std::min<std::size_t>(
            static_cast<std::size_t>(maxTests),
            sampleCount - 1));

    for (std::size_t i = 0;
         i + 1 < sampleCount &&
         testIndices.size() < static_cast<std::size_t>(maxTests);
         ++i)
    {
        const double midpoint = 0.5 * (sampleQ[i] + sampleQ[i + 1]);

        const auto upper =
            std::lower_bound(qGrid.begin(), qGrid.end(), midpoint);

        std::size_t index = 0;

        if (upper == qGrid.end())
        {
            index = qGrid.size() - 1;
        }
        else if (upper != qGrid.begin())
        {
            const std::size_t hi =
                static_cast<std::size_t>(upper - qGrid.begin());
            const std::size_t lo = hi - 1;

            index = std::abs(qGrid[lo] - midpoint) <=
                    std::abs(qGrid[hi] - midpoint)
                ? lo
                : hi;
        }

        if (testIndices.empty() || testIndices.back() != index)
            testIndices.push_back(index);
    }

    const double cx = centre[0];
    const double cy = centre[1];
    const double cz = centre[2];

    const std::size_t nPoints = shellPoints.size();

    std::vector<double> radii(nPoints);
    std::vector<double> volumes(nPoints);

    std::vector<std::complex<double>> conjugateY(
        nPoints * coefficientCount);

    // Precompute exactly the same geometry-dependent terms as the
    // production hydration calculation.
    for (std::size_t p = 0; p < nPoints; ++p)
    {
        const auto& point = shellPoints[p];

        const double x = point.x - cx;
        const double y = point.y - cy;
        const double z = point.z - cz;

        const double r = std::sqrt(x * x + y * y + z * z);

        radii[p] = r;
        volumes[p] = point.volume;

        double theta = 0.0;
        double phi = 0.0;

        if (r > 1.0e-12)
        {
            theta = std::acos(
                std::max(-1.0, std::min(1.0, z / r)));
            phi = std::atan2(y, x);
        }

        const std::size_t base = p * coefficientCount;

        for (int l = 0; l <= L_; ++l)
        {
            for (int m = -l; m <= l; ++m)
            {
                const std::size_t index =
                    static_cast<std::size_t>(lm_index(l, m));

                conjugateY[base + index] =
                    std::conj(spherical_harmonic(l, m, theta, phi));
            }
        }
    }

    std::vector<std::complex<double>> phases(
        static_cast<std::size_t>(L_ + 1));

    const std::complex<double> imaginaryUnit(0.0, 1.0);

    for (int l = 0; l <= L_; ++l)
        phases[static_cast<std::size_t>(l)] =
            std::pow(imaginaryUnit, l);

    std::vector<double> besselValues(
        static_cast<std::size_t>(L_ + 1));

    double maxAbsoluteIntensityError = 0.0;
    double maxRelativeIntensityError = 0.0;
    double sumSquaredRelativeErrors = 0.0;
    std::size_t compared = 0;

    const double G = 1.0 + 3.0 * radiusAdjustmentFraction;

    const auto combinedIntensity =
        [&](std::size_t iq,
            const std::vector<std::complex<double>>& hydration)
    {
        double intensity = 0.0;

        for (std::size_t c = 0; c < coefficientCount; ++c)
        {
            const std::complex<double> total =
                atomic[iq][c]
                - solventDensity * G * excluded[iq][c]
                + shellContrast * hydration[c];

            intensity += std::norm(total);
        }

        return intensity / (4.0 * PI);
    };

    std::cout << "\n=== COMBINED INTENSITY INTERPOLATION VALIDATION ===\n"
              << "L: " << L_ << '\n'
              << "Sparse q samples: " << sampleCount << '\n'
              << "Test q values: " << testIndices.size() << '\n'
              << "q | interpolated I | direct I | abs error | relative error\n";

    for (const std::size_t iq : testIndices)
    {
        const double q = qGrid[iq];

        // Calculate the hydration coefficients directly at this q.
        // Crucially, this local vector does not modify multipoles_.
        std::vector<std::complex<double>> direct(
            coefficientCount, {0.0, 0.0});

        for (std::size_t p = 0; p < nPoints; ++p)
        {
            spherical_bessel_all(
                L_, q * radii[p], besselValues);

            const std::size_t base = p * coefficientCount;

            for (int l = 0; l <= L_; ++l)
            {
                const std::complex<double> radial =
                    4.0 * PI
                    * phases[static_cast<std::size_t>(l)]
                    * volumes[p]
                    * besselValues[static_cast<std::size_t>(l)];

                for (int m = -l; m <= l; ++m)
                {
                    const std::size_t index =
                        static_cast<std::size_t>(lm_index(l, m));

                    direct[index] += radial * conjugateY[base + index];
                }
            }
        }
// Compare hydration coefficients before applying shell contrast.
double maxCoefficientError = 0.0;
double maxDirectCoefficient = 0.0;

double interpolatedHydrationNorm2 = 0.0;
double directHydrationNorm2 = 0.0;

for (std::size_t c = 0; c < coefficientCount; ++c)
{
    maxCoefficientError = std::max(
        maxCoefficientError,
        std::abs(multipoles_[iq][c] - direct[c]));

    maxDirectCoefficient = std::max(
        maxDirectCoefficient,
        std::abs(direct[c]));

    interpolatedHydrationNorm2 +=
        std::norm(multipoles_[iq][c]);

    directHydrationNorm2 +=
        std::norm(direct[c]);
}

    const double coefficientErrorRelative =
        maxCoefficientError /
        std::max(maxDirectCoefficient, 1.0e-30);

    const double interpolatedHydrationNorm =
        std::sqrt(interpolatedHydrationNorm2);

    const double directHydrationNorm =
        std::sqrt(directHydrationNorm2);

    // std::cout
    //     << "Hydration coefficients at q=" << q << '\n'
    //     << "  Max coefficient error (relative): "
    //     << coefficientErrorRelative << '\n'
    //     << "  Interpolated coefficient norm: "
    //     << interpolatedHydrationNorm << '\n'
    //     << "  Direct coefficient norm: "
    //     << directHydrationNorm << '\n'
    //     << "  Contrast-scaled interpolated norm: "
    //     << std::abs(shellContrast) * interpolatedHydrationNorm << '\n'
    //     << "  Contrast-scaled direct norm: "
    //     << std::abs(shellContrast) * directHydrationNorm << '\n';

        
        const double interpolatedIntensity =
            combinedIntensity(iq, multipoles_[iq]);

        const double directIntensity =
            combinedIntensity(iq, direct);

        const double absoluteError =
            std::abs(interpolatedIntensity - directIntensity);

        const double relativeError =
            absoluteError /
            std::max(std::abs(directIntensity), 1.0e-30);

        maxAbsoluteIntensityError =
            std::max(maxAbsoluteIntensityError, absoluteError);

        maxRelativeIntensityError =
            std::max(maxRelativeIntensityError, relativeError);

        sumSquaredRelativeErrors += relativeError * relativeError;
        ++compared;

        std::cout << q << " | "
                  << interpolatedIntensity << " | "
                  << directIntensity << " | "
                  << absoluteError << " | "
                  << relativeError << '\n';
    }

    std::cout << "Maximum absolute intensity error: "
              << maxAbsoluteIntensityError << '\n'
              << "Maximum relative intensity error: "
              << maxRelativeIntensityError << '\n';

    if (compared > 0)
    {
        std::cout << "RMS relative intensity error: "
                  << std::sqrt(
                         sumSquaredRelativeErrors /
                         static_cast<double>(compared))
                  << '\n';
    }

    std::cout << "==================================================\n";
}