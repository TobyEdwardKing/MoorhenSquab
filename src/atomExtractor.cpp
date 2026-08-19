#include "atomExtractor.hpp"

#include <gemmi/model.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/cif.hpp>

#include <utility>
#include <algorithm>
#include <cstdio>
#include <exception>
#include <sstream>

static bool ends_with(
    const std::string& value,
    const std::string& ending)
{
    if (ending.size() > value.size())
        return false;

    return std::equal(
        ending.rbegin(),
        ending.rend(),
        value.rbegin());
}


// static void extract_structure(
//     const gemmi::Model& model,
//     std::vector<AtomData>& atoms,
//     const ProteinAtomTypeTable& proteinAtomTypes)
// {
//     for (const auto& chain : model.chains)
//     {
//         for (const auto& residue : chain.residues)
//         {
//             const std::string residueName = residue.name;

//             for (const auto& atom : residue.atoms)
//             {
//                 const std::string atomName = atom.name;

//                 std::string atomType =
//                     proteinAtomTypes.lookup(
//                         residueName,
//                         atomName
//                     );

//                 atoms.push_back({
//                     static_cast<float>(atom.pos.x),
//                     static_cast<float>(atom.pos.y),
//                     static_cast<float>(atom.pos.z),

//                     static_cast<uint32_t>(
//                         atom.element.atomic_number() - 1
//                     ),

//                     residueName,
//                     atomName,
//                     atomType
//                 });
//             }
//         }
//     }
// }
static void extract_structure(
    const gemmi::Model& model,
    std::vector<AtomData>& atoms,
    const ProteinAtomTypeTable& proteinAtomTypes)
{
    for (const auto& chain : model.chains)
    {
        for (const auto& residue : chain.residues)
        {
            const std::string residueName = residue.name;

            // Gemmi tells us what kind of structural entity this residue
            // belongs to: Polymer, NonPolymer, Branched, Water, etc.
            const bool isPolymer =
                residue.entity_type == gemmi::EntityType::Polymer;

            // Ignore water for now.
            //
            // If you later want explicit solvent/hydration-shell treatment,
            // handle it separately.
            if (residue.entity_type == gemmi::EntityType::Water)
                continue;
            printf(
                    "Residue %-4s entity_type=%d het_flag=%c\n",
                    residue.name.c_str(),
                    static_cast<int>(residue.entity_type),
                    residue.het_flag ? residue.het_flag : ' '
                );
            for (const auto& atom : residue.atoms)
            {
                const std::string atomName = atom.name;

                std::string atomType;

            if (isPolymer)
                {
                    if (proteinAtomTypes.contains(residueName, atomName))
                    {
                        atomType =
                            proteinAtomTypes.lookup(
                                residueName,
                                atomName
                            );
                    }
                    else if (atomName == "OXT")
                    {
                        // Legitimate C-terminal atom.
                        // Leave atomType empty so scattering uses elemental FF.
                    }
                    else
                    {
                        throw std::runtime_error(
                            "Unknown polymer atom: " +
                            residueName + " " +
                            atomName
                        );
                    }
                }

                // For non-polymer atoms atomType deliberately remains empty.
                //
                // Their elemental identity comes directly from Gemmi:
                //
                //     atom.element.atomic_number()
                //
                // This handles arbitrary ligand naming conventions such as:
                //
                //     MYR: C1, C2, O1, ...
                //     9AZ: C1, C2, N7, O18, ...
                //
                // without needing ligand-specific atom-name tables.

                atoms.push_back({
                    static_cast<float>(atom.pos.x),
                    static_cast<float>(atom.pos.y),
                    static_cast<float>(atom.pos.z),

                    static_cast<uint32_t>(
                        atom.element.atomic_number() - 1
                    ),

                    residueName,
                    atomName,
                    atomType
                });
            }
        }
    }
}

std::vector<AtomData>
AtomExtractor::extract_atoms_from_string(
    const std::string& contents,
    const std::string& filename,
    const ProteinAtomTypeTable& proteinAtomTypes)
{
    gemmi::Structure structure;


    try
    {
        if (ends_with(filename, ".pdb") ||
            ends_with(filename, ".ent"))
        {
            structure =
                gemmi::read_pdb_string(
                    contents,
                    filename);
        }
        else
        {
            auto doc = gemmi::cif::read_string(contents);

            printf("Number of CIF blocks: %zu\n", doc.blocks.size());

            structure = gemmi::make_structure(std::move(doc));

            printf("Structure created successfully\n");
        }
    }
    catch (const std::exception& e)
    {
        printf("Gemmi exception: %s\n", e.what());
        return {};
    }


    std::vector<AtomData> atoms;


    for (const auto& model : structure.models)
    {
        extract_structure(
            model,
            atoms,
            proteinAtomTypes
        );
    }


    return atoms;
}

std::string atom_summary_from_string(
    const std::string& contents,
    const std::string& filename,
    const ProteinAtomTypeTable& proteinAtomTypes)
{
    auto atoms =
        AtomExtractor::extract_atoms_from_string(
            contents,
            filename,
            proteinAtomTypes);

    std::ostringstream out;

    out << "Number of atoms: "
        << atoms.size()
        << "\n\n";

    out << "First "
        << std::min<size_t>(10, atoms.size())
        << " atoms:\n\n";

    for (size_t i = 0;
         i < std::min<size_t>(10, atoms.size());
         ++i)
    {
        const auto& a = atoms[i];

        out
            << i
            << "  Z=" << a.Z + 1
            << "  ("
            << a.x << ", "
            << a.y << ", "
            << a.z << ")"
            << "\n";
    }

    return out.str();
}