#ifndef ATOMEXTRACTOR_HPP
#define ATOMEXTRACTOR_HPP

#include <cstdint>
#include <string>
#include <vector>
#include "proteinAtomType.hpp"
struct AtomData
{
    float x;
    float y;
    float z;

    uint32_t Z;
    std::string residueName;    // ARG
    std::string atomName;       // NH1
    std::string atomType;       // N-H+
};


class AtomExtractor
{
public:

    static std::vector<AtomData>
    extract_atoms_from_string(
        const std::string& contents,
        const std::string& filename,
        const ProteinAtomTypeTable& proteinAtomTypes);
};

// Free function declaration
std::string atom_summary_from_string(
    const std::string& contents,
    const std::string& filename,
    const ProteinAtomTypeTable& proteinAtomTypes);

#endif