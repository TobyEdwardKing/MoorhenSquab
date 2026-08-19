#ifndef PROTEIN_ATOM_TYPE_HPP
#define PROTEIN_ATOM_TYPE_HPP

#include <string>
#include <unordered_map>

class ProteinAtomTypeTable
{
public:

    bool load_from_string(
        const std::string& jsonString);

    bool contains(
        const std::string& residueName,
        const std::string& atomName) const;

    std::string lookup(
        const std::string& residueName,
        const std::string& atomName) const;

private:

    std::unordered_map<
        std::string,
        std::unordered_map<std::string, std::string>
    > residueTable_;
};

#endif