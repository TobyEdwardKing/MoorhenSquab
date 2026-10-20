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
#include <stdexcept>
#include <map>

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


/*
 * Extract atoms from a Gemmi model.
 *
 * Important:
 *   - Atom identity is determined from atom.element, NOT from
 *     residue/atom names.
 *   - This means arbitrary valid PDB atom names are accepted:
 *
 *       GLU HN
 *       GLU HA
 *       GLU HB1
 *       GLU OE1
 *       ALA CA
 *       ...
 *
 *   - Water is ignored.
 *   - Atoms with an unknown elemental identity are skipped.
 *
 * The ProteinAtomTypeTable argument is retained in the function
 * signature for compatibility with the existing API, but is not
 * used here. Elemental atomic form factors do not require
 * residue-specific atom names.
 */

static void extract_structure(
    const gemmi::Model& model,
    std::vector<AtomData>& atoms,
    const ProteinAtomTypeTable& proteinAtomTypes)
{
    // Exclude explicit hydrogens because PEPSI unified atom types
    // already account for attached hydrogens.
    constexpr bool includeHydrogens = false;

    for (const auto& chain : model.chains)
    {
        // Find the last polymer residue in this chain.
        // Do not use residue numbering: numbering can contain gaps.
        const gemmi::Residue* lastPolymerResidue = nullptr;

        for (const auto& residue : chain.residues)
        {
            if (residue.entity_type != gemmi::EntityType::Water &&
                residue.entity_type == gemmi::EntityType::Polymer)
            {
                lastPolymerResidue = &residue;
            }
        }

        for (const auto& residue : chain.residues)
        {
            const std::string residueName = residue.name;

            if (residue.entity_type == gemmi::EntityType::Water)
                continue;

            // Apply this rule only to the final polymer residue and
            // only when both terminal oxygen names are present.
            bool hasOT1 = false;
            bool hasOT2 = false;

            if (&residue == lastPolymerResidue)
            {
                for (const auto& atom : residue.atoms)
                {
                    if (atom.name == "OT1") hasOT1 = true;
                    if (atom.name == "OT2") hasOT2 = true;
                }
            }

            const bool hasTerminalOxygenPair = hasOT1 && hasOT2;

            for (const auto& atom : residue.atoms)
            {
                const std::string atomName = atom.name;

                const unsigned int atomicNumber =
                    atom.element.atomic_number();

                if (atomicNumber == 0)
                {
                    printf(
                        "Skipping atom with unknown element: %s %s\n",
                        residueName.c_str(),
                        atomName.c_str());
                    continue;
                }

                if (!includeHydrogens && atomicNumber == 1)
                    continue;

                const uint32_t zIndex =
                    static_cast<uint32_t>(atomicNumber - 1);

                std::string atomType;

                // Special case: terminal carboxyl oxygen pair.
                // Both oxygens use the existing O(resonance) type.
                if (hasTerminalOxygenPair &&
                    (atomName == "OT1" || atomName == "OT2"))
                {
                    atomType = "O(resonance)";
                }
                // Otherwise retain the existing residue/atom lookup.
                else if (atomicNumber != 1 &&
                         proteinAtomTypes.contains(
                             residueName, atomName))
                {
                    atomType = proteinAtomTypes.lookup(
                        residueName, atomName);
                }

                atoms.push_back({
                    static_cast<float>(atom.pos.x),
                    static_cast<float>(atom.pos.y),
                    static_cast<float>(atom.pos.z),
                    zIndex,
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
        /*
         * Determine the input format from the filename.
         */
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
            /*
             * Assume CIF/mmCIF for anything that isn't PDB/ENT.
             */
            auto doc =
                gemmi::cif::read_string(contents);

            printf(
                "Number of CIF blocks: %zu\n",
                doc.blocks.size());

            structure =
                gemmi::make_structure(std::move(doc));

            printf(
                "Structure created successfully\n");
        }
    }
    catch (const std::exception& e)
    {
        printf(
            "Gemmi exception: %s\n",
            e.what());

        return {};
    }


    std::vector<AtomData> atoms;


    /*
     * Extract atoms from every model.
     *
     * NOTE:
     * For an ensemble/NMR structure this will currently concatenate
     * every model. If you later want one model only, this should be
     * changed.
     */
    for (const auto& model : structure.models)
    {
        extract_structure(
            model,
            atoms,
            proteinAtomTypes
        );
    }

    // Diagnostic: identify heavy atoms without a PEPSI unified type.
    std::map<std::string, size_t> unmappedHeavyAtoms;

    for (const auto& atom : atoms)
    {
        if (atom.atomType.empty())
        {
            const std::string key =
                atom.residueName + ":" + atom.atomName;

            ++unmappedHeavyAtoms[key];
        }
    }

    std::printf(
            "Heavy atoms without a PEPSI unified type: %zu unique names\n",
            unmappedHeavyAtoms.size()
        );

        for (const auto& entry : unmappedHeavyAtoms)
        {
            std::printf(
                "  %s: %zu atoms\n",
                entry.first.c_str(),
                entry.second
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