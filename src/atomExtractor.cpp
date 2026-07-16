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


static void extract_structure(
    const gemmi::Model& model,
    std::vector<AtomData>& atoms)
{
    for (const auto& chain : model.chains)
    {
        for (const auto& residue : chain.residues)
        {
            for (const auto& atom : residue.atoms)
            {
                atoms.push_back({
                    static_cast<float>(atom.pos.x),
                    static_cast<float>(atom.pos.y),
                    static_cast<float>(atom.pos.z),
                    static_cast<uint32_t>(
                        atom.element.atomic_number() - 1)
                });
            }
        }
    }
}


std::vector<AtomData>
AtomExtractor::extract_atoms_from_string(
    const std::string& contents,
    const std::string& filename)
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
            atoms);
    }


    return atoms;
}

std::string atom_summary_from_string(
    const std::string& contents,
    const std::string& filename)
{
    auto atoms =
        AtomExtractor::extract_atoms_from_string(
            contents,
            filename);

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