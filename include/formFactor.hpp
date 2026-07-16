#ifndef FORMFACTOR_HPP
#define FORMFACTOR_HPP

#include <vector>
#include <array>
#include <string>
#include <unordered_map>
#include <stdexcept>


struct FormFactor
{
    std::array<double,5> a{};
    std::array<double,5> b{};
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

    bool load_from_string(
        const std::string& jsonContents);


    bool contains(
        uint32_t index) const;


    const FormFactor& operator[](
        uint32_t index) const;


    size_t size() const;

    std::string test_lookup(
        const std::vector<uint32_t>& atomicNumbers
    ) const;
private:

    std::vector<FormFactor> table_;

};

#endif
