// Extra Embind entry points compiled into the npm build only.
//
// They wrap Moorhen Squab's functions from src/main.cpp in try/catch so a C++
// exception comes back to JavaScript as plain text. (Emscripten 6 no longer
// exports getExceptionMessage via -sEXPORT_EXCEPTION_HANDLING_HELPERS, so
// without this a thrown std::runtime_error reaches JS as an opaque object.)
//
// Contract: each squab_* call sets the last error ("" on success); the JS
// wrapper reads it with squab_last_error() and throws if it is non-empty.

#include <emscripten/bind.h>

#include <exception>
#include <type_traits>
#include <string>
#include <vector>

// Defined in ../src/main.cpp (external linkage).
void prepare_debye_gpu(const std::string& elementJson,
                       const std::string& saxsAtomTypeJson,
                       const std::string& residueAtomTypeJson,
                       const std::string& contents,
                       const std::string& filename,
                       const std::vector<float>& qValues);
std::vector<float> get_multipole_intensities(int requestedL);
std::vector<float> get_combined_multipole_intensities(double solventDensity);
std::vector<float> get_hydration_combined_multipole_intensities(
    double solventDensity, double shellContrast, double radiusAdjustmentFraction);

namespace {

std::string lastError;

template <class F>
auto guarded(F&& f) -> decltype(f())
{
    lastError.clear();
    try {
        return f();
    } catch (const std::exception& e) {
        lastError = e.what();
    } catch (...) {
        lastError = "unknown C++ exception";
    }
    if constexpr (!std::is_void_v<decltype(f())>) {
        return {};
    }
}

void squab_prepare(const std::string& elementJson,
                   const std::string& saxsAtomTypeJson,
                   const std::string& residueAtomTypeJson,
                   const std::string& contents,
                   const std::string& filename,
                   const std::vector<float>& qValues)
{
    guarded([&] {
        prepare_debye_gpu(elementJson, saxsAtomTypeJson, residueAtomTypeJson,
                          contents, filename, qValues);
    });
}

std::vector<float> squab_multipole_intensities(int requestedL)
{
    return guarded([&] { return get_multipole_intensities(requestedL); });
}

std::vector<float> squab_combined_multipole_intensities(double solventDensity)
{
    return guarded([&] { return get_combined_multipole_intensities(solventDensity); });
}

std::vector<float> squab_hydration_combined_multipole_intensities(
    double solventDensity, double shellContrast, double radiusAdjustmentFraction)
{
    return guarded([&] {
        return get_hydration_combined_multipole_intensities(
            solventDensity, shellContrast, radiusAdjustmentFraction);
    });
}

std::string squab_last_error()
{
    return lastError;
}

} // namespace

EMSCRIPTEN_BINDINGS(squab_npm)
{
    emscripten::function("squab_prepare", &squab_prepare);
    emscripten::function("squab_multipole_intensities", &squab_multipole_intensities);
    emscripten::function("squab_combined_multipole_intensities", &squab_combined_multipole_intensities);
    emscripten::function("squab_hydration_combined_multipole_intensities",
                         &squab_hydration_combined_multipole_intensities);
    emscripten::function("squab_last_error", &squab_last_error);
}
