#include <emscripten/bind.h>

#include "atomExtractor.hpp"
#include "formFactor.hpp"
#include "ApplicationState.hpp"
#include "DebyeCalculator.hpp"
#include <algorithm>
#include <cmath>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include "proteinAtomType.hpp"
#include "excludedVolumeModel.hpp"
#include "scatteringModel.hpp"
#include "hydrationShell.hpp"
#include <chrono>
ApplicationState app;
DebyeCalculator calculator;
ProteinAtomTypeTable proteinAtomTypes;
ScatteringModel scatteringModel;
ExcludedVolumeModel excludedVolumeModel;
HydrationShell hydrationShell;
std::vector<HydrationPoint> hydrationPoints;

std::vector<float> get_combined_multipole_intensities(
    double solventDensity
)
{
    const auto& atomic = scatteringModel.multipoles();
    const auto& excluded = excludedVolumeModel.coefficients();
    if (atomic.empty()) {
        throw std::runtime_error(
            "Atomic multipoles have not been calculated"
        );
    }

    if (excluded.empty()) {
        throw std::runtime_error(
            "Excluded-volume multipoles have not been calculated"
        );
    }
    if (atomic.size() != excluded.size()) {
        throw std::runtime_error(
            "Atomic and excluded-volume q-grid sizes do not match"
        );
    }

    std::vector<float> intensities;
    intensities.reserve(atomic.size());

    constexpr double FOUR_PI = 4.0 * 3.14159265358979323846;

    for (std::size_t qIndex = 0; qIndex < atomic.size(); ++qIndex) {
        const auto& atomicRow = atomic[qIndex];
        const auto& excludedRow = excluded[qIndex];

        if (atomicRow.size() != excludedRow.size()) {
            throw std::runtime_error(
                "Atomic and excluded-volume multipole counts do not match"
            );
        }

        double intensity = 0.0;

        for (std::size_t lm = 0; lm < atomicRow.size(); ++lm) {
            const std::complex<double> total =
                atomicRow[lm] - solventDensity * excludedRow[lm];

            intensity += std::norm(total);
        }

        intensity /= FOUR_PI;
        intensities.push_back(static_cast<float>(intensity));
    }

    return intensities;
}
void prepare_debye_gpu(
    const std::string& elementJson,
    const std::string& saxsAtomTypeJson,
    const std::string& residueAtomTypeJson,
    const std::string& contents,
    const std::string& filename,
    const std::vector<float>& qValues)
{
    std::cout << "ENTERED prepare_debye_gpu\n";

    FormFactorTable table;

    if (!table.load_element_table(elementJson))
    {
        throw std::runtime_error(
            "Failed loading element form factors");
    }

    if (!table.load_atomtype_table(saxsAtomTypeJson))
    {
        throw std::runtime_error(
            "Failed loading SAXS atom type form factors");
    }

    if (!proteinAtomTypes.load_from_string(residueAtomTypeJson))
    {
        throw std::runtime_error(
            "Failed loading residue atom types");
    }

    app.atoms =
        AtomExtractor::extract_atoms_from_string(
            contents,
            filename,
            proteinAtomTypes);

    /*
     * ------------------------------------------------------------
     * Convert q-grid to double precision
     * ------------------------------------------------------------
     */

    std::vector<double> qGrid(
        qValues.begin(),
        qValues.end()
    );

    /*
     * ------------------------------------------------------------
     * Atomic multipole calculation
     * ------------------------------------------------------------
     *
     * Calculate the full multipole expansion once at L = 12.
     * Lower orders can be reconstructed from these coefficients.
     *
     * The existing GPU Debye calculation remains the source
     * of the displayed SAXS curve.
     */
    constexpr int MULTIPOLE_L = 12;

    scatteringModel.compute_multipoles(
        app.atoms,
        table,
        qGrid,
        MULTIPOLE_L
    );

    excludedVolumeModel.compute_coefficients(
        app.atoms,
        table,
        qGrid,
        MULTIPOLE_L
    );

    std::cout
        << "Atomic multipoles calculated:\n"
        << "  Atoms: " << app.atoms.size() << "\n"
        << "  q points: " << qGrid.size() << "\n"
        << "  L: " << MULTIPOLE_L << "\n"
        << "  coefficients per q: "
        << (MULTIPOLE_L + 1) * (MULTIPOLE_L + 1)
        << "\n";

    /*
     * ------------------------------------------------------------
     * Excluded-volume multipoles across the full q-grid
     * ------------------------------------------------------------
     *
     * Use the same q-grid and maximum L as the atomic multipoles.
     * The coefficients are retained in the global
     * excludedVolumeModel and used by the combined-intensity method.
     */

    constexpr int EXCLUDED_VOLUME_L = MULTIPOLE_L;

    const auto& excludedVolumeCoefficients =
        excludedVolumeModel.coefficients();
    if (excludedVolumeCoefficients.size() != qGrid.size()) {
        throw std::runtime_error(
            "Excluded-volume q-grid size does not match requested q-grid"
        );
    }

    const std::size_t expectedCoefficientCount =
        (MULTIPOLE_L + 1) * (MULTIPOLE_L + 1);

    for (std::size_t i = 0; i < excludedVolumeCoefficients.size(); ++i) {
        if (excludedVolumeCoefficients[i].size() != expectedCoefficientCount) {
            throw std::runtime_error(
                "Incorrect excluded-volume coefficient count at q index " +
                std::to_string(i)
            );
        }
    }
    std::cout
        << "Excluded-volume multipoles calculated:\n"
        << "  Atoms: " << app.atoms.size() << "\n"
        << "  q points requested: " << qGrid.size() << "\n"
        << "  q rows stored: "
        << excludedVolumeCoefficients.size() << "\n"
        << "  L: " << EXCLUDED_VOLUME_L << "\n"
        << "  coefficients per q: "
        << (EXCLUDED_VOLUME_L + 1) * (EXCLUDED_VOLUME_L + 1)
        << "\n";

    /*
     * ------------------------------------------------------------
     * Preserve the q = 0 excluded-volume normalisation diagnostic
     * ------------------------------------------------------------
     *
     * The experimental q-grid does not start at zero, so evaluate
     * q = 0 separately. This does not replace the full-grid result.
     */

    const std::vector<double> testQ = {0.0};

    ExcludedVolumeModel excludedVolumeZeroModel;

    excludedVolumeZeroModel.compute_coefficients(
        app.atoms,
        table,
        testQ,
        EXCLUDED_VOLUME_L
    );

    const auto& zeroCoefficients =
        excludedVolumeZeroModel.coefficients();

    if (!zeroCoefficients.empty() &&
        !zeroCoefficients.front().empty())
    {
        const auto& c00 = zeroCoefficients.front()[0];

        std::cout
            << "Excluded-volume q=0 normalisation check:\n"
            << "  C00: ("
            << c00.real() << ", "
            << c00.imag() << ")\n";
    }
    /*
    * ------------------------------------------------------------
    * Generate hydration-shell geometry
    * ------------------------------------------------------------
    *
    * This stage generates the shell points and reports their
    * represented volume. It does not yet calculate hydration
    * multipoles or alter the displayed SAXS curve.
    */

    constexpr double HYDRATION_SHELL_THICKNESS = 5.0; // Angstroms
    constexpr double HYDRATION_GRID_SPACING = 1.0;     // Angstroms

    hydrationPoints = hydrationShell.generate(
        app.atoms,
        table,
        HYDRATION_SHELL_THICKNESS,
        HYDRATION_GRID_SPACING
    );

    const double hydrationVolume =
        HydrationShell::total_volume(hydrationPoints);

    std::cout
        << "\n===== HYDRATION SHELL GEOMETRY =====\n"
        << "  Protein atoms: " << app.atoms.size() << "\n"
        << "  Shell thickness: "
        << HYDRATION_SHELL_THICKNESS << " A\n"
        << "  Grid spacing: "
        << HYDRATION_GRID_SPACING << " A\n"
        << "  Shell points: " << hydrationPoints.size() << "\n"
        << "  Represented shell volume: "
        << hydrationVolume << " A^3\n"
        << "====================================\n";

    double cx = 0.0;
    double cy = 0.0;
    double cz = 0.0;

    for (const auto& atom : app.atoms)
    {
        cx += atom.x;
        cy += atom.y;
        cz += atom.z;
    }

    if (app.atoms.empty())
    {
        throw std::runtime_error(
            "Cannot calculate hydration centre: no atoms were extracted"
        );
    }

    cx /= static_cast<double>(app.atoms.size());
    cy /= static_cast<double>(app.atoms.size());
    cz /= static_cast<double>(app.atoms.size());

    const std::vector<double> hydrationCentre = {cx, cy, cz};

    // Benchmark hydration-shell multipoles at representative q-values.
    const std::size_t qCount = qGrid.size();

    if (qCount < 4)
        throw std::runtime_error(
            "Need at least four q-values for hydration benchmark."
        );

    const std::vector<double> benchmarkQ = {
        qGrid[0],
        qGrid[qCount / 3],
        qGrid[(2 * qCount) / 3],
        qGrid[qCount - 1]
    };

    for (int benchmarkL : {4, 8, 12})
    {
        const auto start = std::chrono::steady_clock::now();

        hydrationShell.compute_multipoles(
            hydrationPoints,
            benchmarkQ,
            hydrationCentre,
            benchmarkL
        );

        const auto stop = std::chrono::steady_clock::now();

        const double elapsedMs =
            std::chrono::duration<double, std::milli>(
                stop - start
            ).count();

        std::cout
            << "Hydration benchmark: L=" << benchmarkL
            << ", shell points=" << hydrationPoints.size()
            << ", q points=" << benchmarkQ.size()
            << ", coefficients/q="
            << (benchmarkL + 1) * (benchmarkL + 1)
            << ", elapsed=" << elapsedMs << " ms"
            << std::endl;
    }

    const auto& hydrationCoefficients = hydrationShell.coefficients();

    // Independent q=0 normalisation check.
    const std::vector<double> hydrationZeroQ = {0.0};

    hydrationShell.compute_multipoles(
        hydrationPoints,
        hydrationZeroQ,
        hydrationCentre,
        0
    );

    const auto& hydrationZeroCoefficients =
        hydrationShell.coefficients();

    const double expectedB00 =
        std::sqrt(4.0 * 3.14159265358979323846) *
        hydrationVolume;

    const double actualB00 =
        hydrationZeroCoefficients.at(0).at(0).real();

    const double relativeError =
        std::abs(actualB00 - expectedB00) /
        std::max(std::abs(expectedB00), 1.0e-30);

    std::cout
        << "Hydration q=0 B00 check:\n"
        << "  Expected: " << expectedB00 << "\n"
        << "  Actual: " << actualB00 << "\n"
        << "  Relative error: " << relativeError << '\n';
    /*
     * ------------------------------------------------------------
     * Existing GPU Debye calculation
     * ------------------------------------------------------------
     *
     * Unchanged: this continues to supply the displayed SAXS curve.
     */

    calculator.prepare_buffers(
        app.atoms,
        table,
        qValues);

    std::cout << "positions after prepare = "
              << calculator.get_positions().size()
              << "\n";

    std::cout << "prepare_debye_gpu GPU amplitudes = "
              << calculator.get_amplitudes().size()
              << "\n";

    std::cout << "prepare_debye_gpu atomic inputs = "
              << app.atoms.size()
              << "\n";

    std::cout << "prepare_debye_gpu q values = "
              << qValues.size()
              << "\n";
}


/*
 * ------------------------------------------------------------
 * Return a multipole intensity curve for a requested L.
 *
 * The underlying atomic multipole coefficients were calculated
 * once at L = 12.
 * ------------------------------------------------------------
 */

std::vector<float> get_multipole_intensities(
    int requestedL)
{
    const auto& multipoles =
        scatteringModel.multipoles();

    std::vector<float> intensities;

    intensities.reserve(
        multipoles.size()
    );

    for (size_t iq = 0;
         iq < multipoles.size();
         ++iq)
    {
        const double intensity =
            scatteringModel.intensity_from_multipoles(
                iq,
                requestedL
            );

        intensities.push_back(
            static_cast<float>(intensity)
        );
    }

    return intensities;
}


std::vector<float> get_positions()
{
    return calculator.get_positions();
}


std::vector<float> get_amplitudes()
{
    return calculator.get_amplitudes();
}


EMSCRIPTEN_BINDINGS(saxs)
{
    emscripten::register_vector<float>(
        "FloatVector");

    emscripten::function(
        "prepare_debye_gpu",
        &prepare_debye_gpu);

    emscripten::function(
        "get_positions",
        &get_positions);

    emscripten::function(
        "get_amplitudes",
        &get_amplitudes);

    emscripten::function(
        "get_multipole_intensities",
        &get_multipole_intensities);

    emscripten::function(
        "get_combined_multipole_intensities",
        &get_combined_multipole_intensities
    );
}
