#include "formFactor.hpp"

#include <nlohmann/json.hpp>

#include <sstream>
#include <cmath>
#include <iomanip>
#include <iostream>

using json = nlohmann::json;



double FormFactor::evaluate(double q) const
{
    const double s =
        q / (4.0 * M_PI);

    double result = c;


    for (int i = 0; i < 5; i++)
    {
        result +=
            a[i] *
            std::exp(-b[i] * s * s);
    }

    return result;
}



bool FormFactorTable::load_element_table(
    const std::string& jsonContents)
{
    json j = json::parse(jsonContents);

    element_table_.clear();

    for (auto& [element, value] : j.items())
    {
        if (value.size() != 15)
        {
            std::cerr
                << "Skipping invalid entry: "
                << element
                << std::endl;
            continue;
        }

        FormFactor ff;

        // a1-a5
        for (int i = 0; i < 5; i++)
        {
            ff.a[i] =
                value[i].is_number()
                ? value[i].get<double>()
                : 0.0;
        }

        // c
        ff.c =
            value[5].is_number()
            ? value[5].get<double>()
            : 0.0;

        // b1-b5
        for (int i = 0; i < 5; i++)
        {
            ff.b[i] =
                value[6 + i].is_number()
                ? value[6 + i].get<double>()
                : 0.0;
        }

        ff.coherent_scattering_length =
            value[11].is_number()
            ? value[11].get<double>()
            : 0.0;

        ff.atomic_number =
            value[12].is_number()
            ? value[12].get<int>()
            : 0;

        ff.atomic_radius =
            value[13].is_number()
            ? value[13].get<double>()
            : 0.0;

        ff.excl_radius = ff.atomic_radius;

        ff.excl_volume =
            value[14].is_number()
            ? value[14].get<double>()
            : 0.0;

        // Store using the JSON element symbol, e.g. "O", "H", "C".
        element_table_[element] = ff;
    }

    std::cout
        << "Loaded "
        << element_table_.size()
        << " elemental form factor entries."
        << std::endl;

    return true;
}



bool FormFactorTable::contains_element(
    const std::string& element) const
{
    return element_table_.find(element)
        != element_table_.end();
}



const FormFactor& FormFactorTable::get_element(
    const std::string& element) const
{
    auto it = element_table_.find(element);

    if (it == element_table_.end())
    {
        throw std::runtime_error(
            "Missing elemental form factor: " + element);
    }

    return it->second;
}



size_t FormFactorTable::number_of_elements() const
{
    return element_table_.size();
}

size_t FormFactorTable::number_of_atom_types() const
{
    return atomtype_table_.size();
}

std::string FormFactorTable::test_lookup(
    const std::vector<std::string>& elements
) const
{
    std::ostringstream out;

    out << "Form factor lookup test\n";
    out << "======================\n\n";

    int count = 0;

    for (const std::string& element : elements)
    {
        if (count >= 10)
            break;

        out << "Atom " << count << "\n";

        out << "  Element: "
            << element
            << "\n";

        if (!contains_element(element))
        {
            out << "  ERROR: no form factor found\n\n";
            count++;
            continue;
        }

        const auto& ff = get_element(element);

        out << "  Atomic number: "
            << ff.atomic_number
            << "\n";

        out << "  f(0): "
            << std::fixed
            << std::setprecision(5)
            << ff.evaluate(0.0)
            << "\n";

        out << "  Radius: "
            << ff.atomic_radius
            << "\n\n";

        count++;
    }

    return out.str();
}



bool FormFactorTable::load_atomtype_table(
    const std::string& jsonContents)
{
    json j = json::parse(jsonContents);

    atomtype_table_.clear();

    for (auto& [atomType, value] : j.items())
    {
        FormFactor ff;

        // a coefficients
        if (!value.contains("a") || value["a"].size() != 5)
            throw std::runtime_error(
                "Invalid atom type entry: " + atomType);

        for (int i = 0; i < 5; i++)
            ff.a[i] = value["a"][i].get<double>();

        // b coefficients
        if (!value.contains("b") || value["b"].size() != 5)
            throw std::runtime_error(
                "Invalid atom type entry: " + atomType);

        for (int i = 0; i < 5; i++)
            ff.b[i] = value["b"][i].get<double>();

        ff.c = value["c"].get<double>();

        if (!value.contains("excl_radius") ||
            !value.contains("excl_volume"))
        {
            throw std::runtime_error(
                "Missing excluded-volume parameters for atom type: "
                + atomType);
        }

        ff.excl_radius =
            value["excl_radius"].get<double>();

        ff.excl_volume =
            value["excl_volume"].get<double>();

        atomtype_table_[atomType] = ff;
    }

    std::cout
        << "Loaded "
        << atomtype_table_.size()
        << " PEPSI atom types."
        << std::endl;

    return true;
}

bool FormFactorTable::contains_atom_type(
    const std::string& atomType) const
{
    return atomtype_table_.find(atomType)
        != atomtype_table_.end();
}

const FormFactor&
FormFactorTable::get_atom_type(
    const std::string& atomType) const
{
    auto it = atomtype_table_.find(atomType);

    if (it == atomtype_table_.end())
    {
        throw std::runtime_error(
            "Unknown SAXS atom type: " + atomType);
    }

    return it->second;
}