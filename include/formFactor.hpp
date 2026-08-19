#ifndef FORMFACTOR_HPP
#define FORMFACTOR_HPP

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

struct FormFactor
{
    std::array<double, 5> a{};
    std::array<double, 5> b{};
    double c = 0.0;

    double coherent_scattering_length = 0.0;

    int atomic_number = 0;
    double atomic_radius = 0.0;

    // Excluded solvent volume (Å³)
    double excluded_volume = 0.0;

    double evaluate(double q) const;

    double excluded_amplitude(double q) const;
};

class FormFactorTable
{
public:

    //
    // Loading
    //

    // Loads the existing element_info.json
    bool load_element_table(
        const std::string& jsonContents);

    // Loads the new saxs_atom_types.json
    bool load_atomtype_table(
        const std::string& jsonContents);


    //
    // Element lookup
    //

    bool contains_element(
        uint32_t atomicNumber) const;

    const FormFactor& get_element(
        uint32_t atomicNumber) const;


    //
    // Protein atom-type lookup
    //

    bool contains_atom_type(
        const std::string& atomType) const;

    const FormFactor& get_atom_type(
        const std::string& atomType) const;


    //
    // Diagnostics
    //

    size_t number_of_elements() const;

    size_t number_of_atom_types() const;

    std::string test_lookup(
        const std::vector<uint32_t>& atomicNumbers) const;

private:

    // Indexed by atomic number - 1
    std::vector<FormFactor> element_table_;

    // Indexed by PEPSI atom type string
    std::unordered_map<std::string, FormFactor> atomtype_table_;
};

#endif