#ifndef ATOMEXTRACTOR_HPP
#define ATOMEXTRACTOR_HPP

#include <cstdint>
#include <string>
#include <vector>

struct AtomData
{
    float x;
    float y;
    float z;
    uint32_t Z;
};

class AtomExtractor
{
public:

    static std::vector<AtomData>
    extract_atoms_from_string(
        const std::string& contents,
        const std::string& filename);
};

// Free function declaration
std::string atom_summary_from_string(
    const std::string& contents,
    const std::string& filename);

#endif