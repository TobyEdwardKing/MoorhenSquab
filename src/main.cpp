
#include <emscripten/bind.h>

#include "atomExtractor.hpp"
#include "formFactor.hpp"
#include "ApplicationState.hpp"
#include "DebyeCalculator.hpp"
#include <algorithm>
#include <sstream>
#include <iostream>
#include "proteinAtomType.hpp"


ApplicationState app;
DebyeCalculator calculator;
ProteinAtomTypeTable proteinAtomTypes;
// int atom_count(
//     const std::string& text,
//     const std::string& filename)
// {
//     app.atoms =
//         AtomExtractor::extract_atoms_from_string(
//             text,
//             filename);

//     return app.atoms.size();
// }


// std::string atom_summary(
//     const std::string& text,
//     const std::string& filename)
// {
//     return atom_summary_from_string(
//         text,
//         filename);
// }


// std::string form_factor_summary(
//     const std::string& text,
//     const std::string& filename,
//     const std::string& jsonText)
// {
//     app.atoms =
//         AtomExtractor::extract_atoms_from_string(
//             text,
//             filename);


//     FormFactorTable table;

//     if (!table.load_from_string(jsonText))
//     {
//         return "Failed to load form factors";
//     }


//     std::stringstream output;

//     output << "Form factor test\n";
//     output << "Atoms: "
//            << app.atoms.size()
//            << "\n\n";


//     for (size_t i = 0;
//          i < std::min<size_t>(10, app.atoms.size());
//          i++)
//     {
//         const auto& atom = app.atoms[i];


//         if (!table.contains(atom.Z))
//         {
//             output
//                 << i
//                 << " missing form factor\n";

//             continue;
//         }


//         const auto& ff =
//             table[atom.Z];


//         output
//             << i
//             << " lookup="
//             << atom.Z
//             << " atomic_number="
//             << ff.atomic_number
//             << " f(0)="
//             << ff.evaluate(0.0)
//             << "\n";
//     }


//     return output.str();
// }


// std::string form_factor_test(
//     const std::string& json,
//     const std::string& contents,
//     const std::string& filename)
// {
//     FormFactorTable table;

//     if (!table.load_from_string(json))
//         return "Failed loading form factors";


//     app.atoms =
//         AtomExtractor::extract_atoms_from_string(
//             contents,
//             filename);


//     std::vector<uint32_t> Zs;

//     for (const auto& atom : app.atoms)
//     {
//         Zs.push_back(atom.Z);
//     }


//     return table.test_lookup(Zs);
// }


// void compute_debye_curve(
//     const std::string& json,
//     const std::string& contents,
//     const std::string& filename,
//     double qMin,
//     double qMax,
//     double qStep)
// {
//     FormFactorTable table;

//     if (!table.load_from_string(json))
//         return;

//     app.atoms =
//         AtomExtractor::extract_atoms_from_string(
//             contents,
//             filename);

//     // calculator.compute_curve(
//     //     app.atoms,
//     //     table,
//     //     qMin,
//     //     qMax,
//     //     qStep);
//     // std::cout
//     //     << "Calculator reports "
//     //     << calculator.curve_size()
//     //     << " points\n";
// }

// int curve_size()
// {
//     return calculator.curve_size();
// }

// double curve_q(int i)
// {
//     return calculator.q_value(i);
// }

// double curve_intensity(int i)
// {
//     return calculator.intensity_value(i);
// }

void prepare_debye_gpu(
    const std::string& elementJson,
    const std::string& saxsAtomTypeJson,
    const std::string& residueAtomTypeJson,
    const std::string& contents,
    const std::string& filename,
    double q)
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

    calculator.prepare_buffers(
        app.atoms,
        table,
        q);
    std::cout << "positions after prepare = "
            << calculator.get_positions().size()
            << "\n";

    std::cout
        << "prepare_debye_gpu GPU scatterers = "
        << calculator.get_amplitudes().size()
        << "\n";
    std::cout
        << "prepare_debye_gpu atomic inputs = "
        << app.atoms.size()
        << "\n";
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

    // emscripten::function(
    //     "atom_count",
    //     &atom_count);


    // emscripten::function(
    //     "atom_summary",
    //     &atom_summary);


    // emscripten::function(
    //     "form_factor_summary",
    //     &form_factor_summary);


    // emscripten::function(
    //     "form_factor_test",
    //     &form_factor_test);


    // emscripten::function(
    //     "compute_debye_curve",
    //     &compute_debye_curve);

    // emscripten::function(
    //     "curve_size",
    //     &curve_size);

    // emscripten::function(
    //     "curve_q",
    //     &curve_q);

    // emscripten::function(
    //     "curve_intensity",
    //     &curve_intensity);
    emscripten::function(
        "prepare_debye_gpu",
        &prepare_debye_gpu);

    emscripten::function(
        "get_positions",
        &get_positions);

    emscripten::function(
        "get_amplitudes",
        &get_amplitudes);
}

    