
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
#include <iomanip>

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


std::vector<float> get_hydration_combined_multipole_intensities(
    double solventDensity,
    double shellContrast,
    double radiusAdjustmentFraction
)
{
    // PEPSI's standard search range for delta-r / r_m.
    if (!std::isfinite(radiusAdjustmentFraction) ||
        radiusAdjustmentFraction < -0.05 ||
        radiusAdjustmentFraction > 0.05)
    {
        throw std::runtime_error(
            "Radius adjustment fraction must be between -0.05 and 0.05"
        );
    }

    if (!std::isfinite(solventDensity) ||
        !std::isfinite(shellContrast))
    {
        throw std::runtime_error(
            "Solvent density and shell contrast must be finite"
        );
    }

    const auto& atomic = scatteringModel.multipoles();
    const auto& excluded = excludedVolumeModel.coefficients();
    const auto& hydration = hydrationShell.coefficients();

    if (atomic.empty()) {
        throw std::runtime_error(
            "Atomic multipoles have not been calculated"
        );
    }

    if (excluded.size() != atomic.size() ||
        hydration.size() != atomic.size())
    {
        throw std::runtime_error(
            "Atomic, excluded-volume and hydration q-grid sizes do not match"
        );
    }

    /*
     * PEPSI's simplified, q-independent first-order radius correction.
     *
     * C_lm(q) is calculated once using the reference dummy-atom
     * form factors. Apply the first-order correction to the
     * excluded-volume amplitude here, once.
     */
    const double G = 1.0 + 3.0 * radiusAdjustmentFraction;

    constexpr double FOUR_PI = 4.0 * 3.14159265358979323846;

    /*
     * Water molecular form factor at q=0, in electron units:
     * oxygen contributes 8 electrons and the two hydrogens
     * contribute 1 electron each.
     *
     * HydrationShell::compute_multipoles() already includes the
     * represented grid-cell volume in B_lm(q). Therefore, do not
     * multiply B_lm(q) by grid_spacing()^3 again here.
     *
     * Normalising f_water(q) by f_water(0) gives a dimensionless
     * molecular form-factor shape. shellContrast is then applied
     * exactly once to the volume-weighted shell coefficients.
     */
    constexpr double WATER_FORM_FACTOR_AT_ZERO = 10.0;

    std::vector<float> intensities;
    intensities.reserve(atomic.size());

    for (std::size_t qIndex = 0; qIndex < atomic.size(); ++qIndex)
    {
        const auto& atomicRow = atomic[qIndex];
        const auto& excludedRow = excluded[qIndex];
        const auto& hydrationRow = hydration[qIndex];

        if (atomicRow.size() != excludedRow.size() ||
            atomicRow.size() != hydrationRow.size())
        {
            throw std::runtime_error(
                "Atomic, excluded-volume and hydration multipole counts do not match"
            );
        }

        /*
         * Dimensionless water form-factor shape at this q.
         * At q=0, this ratio must equal 1.
         */
        const double waterFormFactorShape =
            hydrationShell.water_form_factor(qIndex) /
            WATER_FORM_FACTOR_AT_ZERO;

        double intensity = 0.0;

        for (std::size_t lm = 0; lm < atomicRow.size(); ++lm)
        {
            /*
             * Total amplitude:
             *
             * A_lm(q) - rho * G * C_lm(q)
             *          + delta_rho * [f_water(q)/f_water(0)] * B_lm(q)
             *
             * B_lm already contains the shell grid-cell volume.
             * No additional grid-volume or bulk-density multiplier
             * belongs in this expression under this convention.
             */
            const std::complex<double> base =
                atomicRow[lm]
                - solventDensity * G * excludedRow[lm];

            const std::complex<double> hydrationContribution =
                shellContrast *
                waterFormFactorShape *
                hydrationRow[lm];

            const std::complex<double> total =
                base + hydrationContribution;

            intensity += std::norm(total);
        }

        intensity /= FOUR_PI;
        intensities.push_back(static_cast<float>(intensity));
    }

    return intensities;
}

static int selectedMultipoleL = 0;
int get_multipole_L()
{
    return selectedMultipoleL;
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
     * Convert q-grid to double precision.
     */
    std::vector<double> qGrid(
        qValues.begin(),
        qValues.end()
    );

    /*
     * Generate hydration-shell geometry before choosing L.
     */
    constexpr double HYDRATION_GRID_SPACING = 3.5;

    hydrationPoints = hydrationShell.generate(
        app.atoms,
        table,
        HYDRATION_GRID_SPACING
    );

    if (hydrationPoints.empty())
    {
        throw std::runtime_error(
            "Hydration-shell generation returned no points."
        );
    }

    const double hydrationVolume =
        HydrationShell::total_volume(hydrationPoints);

    const double hydrationRg =
        HydrationShell::radius_of_gyration(hydrationPoints);

    // Use the same arithmetic-mean atomic centre as the other
    // multipole calculations.
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
            "Cannot calculate hydration centre: no atoms were extracted."
        );
    }

    cx /= static_cast<double>(app.atoms.size());
    cy /= static_cast<double>(app.atoms.size());
    cz /= static_cast<double>(app.atoms.size());

    const std::vector<double> hydrationCentre = {cx, cy, cz};

    if (qGrid.empty() || qGrid.back() <= 0.0)
    {
        throw std::runtime_error(
            "Cannot select multipole order: invalid q-grid."
        );
    }

    // Estimate L from q_max * Rg of the generated hydration shell.
    const int MULTIPOLE_L = std::max(
        1,
        static_cast<int>(std::ceil(qGrid.back() * hydrationRg))
    );

    selectedMultipoleL = MULTIPOLE_L;

    std::cout
        << "\n===== HYDRATION SHELL GEOMETRY =====\n"
        << "  Grid-cell width: "
        << hydrationShell.grid_spacing() << " A\n"
        << "  Shell thickness: "
        << hydrationShell.shell_thickness() << " A\n"
        << "  Shell points: " << hydrationPoints.size() << "\n"
        << "  Represented shell volume: "
        << hydrationVolume << " A^3\n"
        << "  Shell Rg: " << hydrationRg << " A\n"
        << "  q_max: " << qGrid.back() << " A^-1\n"
        << "  Selected L: " << MULTIPOLE_L << "\n"
        << "====================================\n";

    /*
     * Atomic multipole calculation.
     */
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
     * Validate excluded-volume multipoles.
     */
    const int EXCLUDED_VOLUME_L = MULTIPOLE_L;

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
     * q=0 excluded-volume normalisation diagnostic.
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

    run_spherical_bessel_test();

    // Validate the 2L-sample interpolation against direct coefficient
    // calculations at 16 q-values spanning the experimental range.
    constexpr std::size_t VALIDATION_Q_COUNT = 16;

    std::vector<double> validationQ(VALIDATION_Q_COUNT);

    for (std::size_t i = 0; i < VALIDATION_Q_COUNT; ++i) {
        const double t =
            static_cast<double>(i) /
            static_cast<double>(VALIDATION_Q_COUNT - 1);

        validationQ[i] =
            qGrid.front() + t * (qGrid.back() - qGrid.front());
    }

    using Clock = std::chrono::steady_clock;

    // Existing 2L-sample interpolation.
    const auto interpStart = Clock::now();

    hydrationShell.compute_multipoles(
        hydrationPoints,
        table,
        validationQ,
        hydrationCentre,
        MULTIPOLE_L,
        true);

    const auto interpolated = hydrationShell.coefficients();

    const std::vector<std::vector<std::complex<double>>>
        interpolatedCopy(interpolated);

    const double interpMs =
        std::chrono::duration<double, std::milli>(
            Clock::now() - interpStart).count();

    // Direct calculation at every validation q-value.
    const auto directStart = Clock::now();

    hydrationShell.compute_multipoles(
        hydrationPoints,
        table,
        validationQ,
        hydrationCentre,
        MULTIPOLE_L,
        false);

    const auto direct = hydrationShell.coefficients();

    const std::vector<std::vector<std::complex<double>>>
        directCopy(direct);

    const double directMs =
        std::chrono::duration<double, std::milli>(
            Clock::now() - directStart).count();

    // Compare all complex coefficients across all validation q-values.
    double maxAbsError = 0.0;
    double maxDirectMagnitude = 0.0;
    long double sumErrorSquared = 0.0L;
    long double sumDirectSquared = 0.0L;
    std::size_t compared = 0;

    for (std::size_t iq = 0; iq < validationQ.size(); ++iq) {
        for (std::size_t lm = 0;
            lm < directCopy[iq].size(); ++lm) {

            const double error = std::abs(
                interpolatedCopy[iq][lm] - directCopy[iq][lm]);

            const double reference =
                std::abs(directCopy[iq][lm]);

            maxAbsError = std::max(maxAbsError, error);
            maxDirectMagnitude =
                std::max(maxDirectMagnitude, reference);

            sumErrorSquared +=
                static_cast<long double>(error) * error;

            sumDirectSquared +=
                static_cast<long double>(reference) * reference;

            ++compared;
        }
    }

    const double maxRelativeError =
        maxDirectMagnitude > 0.0
            ? maxAbsError / maxDirectMagnitude
            : maxAbsError;

    const double rmsRelativeError =
        sumDirectSquared > 0.0L
            ? std::sqrt(
                static_cast<double>(
                    sumErrorSquared / sumDirectSquared))
            : std::sqrt(
                static_cast<double>(
                    sumErrorSquared /
                    std::max<std::size_t>(compared, 1)));

    std::cout << std::fixed << std::setprecision(3)
            << "[Hydration validation] L=" << MULTIPOLE_L
            << ", q samples=" << VALIDATION_Q_COUNT
            << ", coefficients=" << compared << '\n'
            << "  Interpolated: " << interpMs << " ms\n"
            << "  Direct:       " << directMs << " ms\n"
            << std::scientific << std::setprecision(3)
            << "  Max error / max direct magnitude: "
            << maxRelativeError << '\n'
            << "  Relative RMS coefficient error:   "
            << rmsRelativeError << '\n';

    /*
     * Independent q=0 hydration normalisation check.
     *
     * This checks the volume-weighted geometric coefficients.
     * It does not validate the water form factor or the physical
     * hydration contrast.
     */
    const std::vector<double> hydrationZeroQ = {0.0};

    hydrationShell.compute_multipoles(
        hydrationPoints,
        table,
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
     * Calculate hydration-shell multipoles on the full q-grid.
     * This must follow the q=0 diagnostic because that diagnostic
     * overwrites the shell's stored coefficients.
     */
    hydrationShell.compute_multipoles(
        hydrationPoints,
        table,
        qGrid,
        hydrationCentre,
        MULTIPOLE_L
    );

    hydrationShell.validate_interpolation(
        hydrationPoints,
        qGrid,
        hydrationCentre,
        25
    );

    const double diagnosticSolventDensity = 0.334;

    const std::vector<double> diagnosticShellContrasts = {
        0.0, 0.005, 0.01, 0.015, 0.02, 0.025, 0.03
    };

    for (double contrast : diagnosticShellContrasts)
    {
        std::cout
            << "\nValidating combined intensity interpolation"
            << " for shell contrast = " << contrast << " e-/A^3\n";

        hydrationShell.validate_combined_intensity_interpolation(
            hydrationPoints,
            qGrid,
            hydrationCentre,
            scatteringModel.multipoles(),
            excludedVolumeModel.coefficients(),
            diagnosticSolventDensity,
            contrast,
            0.0, // radius adjustment fraction
            25
        );
    }

    const auto& fullHydrationCoefficients =
        hydrationShell.coefficients();

    if (fullHydrationCoefficients.size() != qGrid.size())
    {
        throw std::runtime_error(
            "Hydration-shell q-grid size does not match requested q-grid"
        );
    }

    const std::size_t expectedHydrationCoefficientCount =
        (MULTIPOLE_L + 1) * (MULTIPOLE_L + 1);

    for (std::size_t i = 0;
         i < fullHydrationCoefficients.size();
         ++i)
    {
        if (fullHydrationCoefficients[i].size() !=
            expectedHydrationCoefficientCount)
        {
            throw std::runtime_error(
                "Incorrect hydration-shell coefficient count at q index " +
                std::to_string(i)
            );
        }
    }

    std::cout
        << "Hydration-shell multipoles calculated:\n"
        << "  Shell points: " << hydrationPoints.size() << "\n"
        << "  q points: " << fullHydrationCoefficients.size() << "\n"
        << "  L: " << MULTIPOLE_L << "\n"
        << "  coefficients per q: "
        << expectedHydrationCoefficientCount << "\n";

    /*
     * Existing GPU Debye calculation remains unchanged.
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
 * Return a multipole intensity curve for a requested L.
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

    emscripten::function(
        "get_hydration_combined_multipole_intensities",
        &get_hydration_combined_multipole_intensities
    );

    emscripten::function(
        "get_multipole_L",
        &get_multipole_L
    );
}
