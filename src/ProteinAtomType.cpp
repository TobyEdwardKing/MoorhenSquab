#include "proteinAtomType.hpp"

#include <nlohmann/json.hpp>

using json = nlohmann::json;


bool ProteinAtomTypeTable::load_from_string(
    const std::string& jsonString)
{
    residueTable_.clear();

    json j = json::parse(jsonString);

    for (const auto& [residueName, residueData] : j.items())
    {
        std::unordered_map<std::string, std::string> atomTable;

        for (const auto& [atomName, atomType] : residueData.items())
        {
            atomTable.emplace(
                atomName,
                atomType.get<std::string>()
            );
        }

        residueTable_.emplace(
            residueName,
            std::move(atomTable)
        );
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
        return "";
    
    auto atomIt = residueIt->second.find(atomName);

    if (atomIt == residueIt->second.end())
        return "";

    return atomIt->second;
}