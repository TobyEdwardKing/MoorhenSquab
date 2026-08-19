#include "proteinAtomType.hpp"

#include <nlohmann/json.hpp>

#include <stdexcept>

using json = nlohmann::json;

bool ProteinAtomTypeTable::load_from_string(
    const std::string& jsonContents)
{
    json j = json::parse(jsonContents);

    residueTable_.clear();

    for (const auto& [residueName, residueData] : j.items())
    {
        std::unordered_map<std::string, std::string> atomMap;

        for (const auto& [atomName, atomType] : residueData.items())
        {
            atomMap.emplace(
                atomName,
                atomType.get<std::string>());
        }

        residueTable_.emplace(
            residueName,
            std::move(atomMap));
    }

    return true;
}

bool ProteinAtomTypeTable::contains(
    const std::string& residueName,
    const std::string& atomName) const
{
    auto residueIt = residueTable_.find(residueName);

    if (residueIt == residueTable_.end())
        return false;

    return residueIt->second.find(atomName)
        != residueIt->second.end();
}

std::string ProteinAtomTypeTable::lookup(
    const std::string& residueName,
    const std::string& atomName) const
{
    auto residueIt = residueTable_.find(residueName);

    if (residueIt == residueTable_.end())
    {
        throw std::runtime_error(
            "Unknown residue: " + residueName);
    }

    auto atomIt = residueIt->second.find(atomName);

    if (atomIt == residueIt->second.end())
    {
        throw std::runtime_error(
            "Unknown atom '" + atomName +
            "' in residue '" + residueName + "'");
    }

    return atomIt->second;
}