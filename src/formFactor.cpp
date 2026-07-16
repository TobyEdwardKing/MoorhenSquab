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

double FormFactor::excluded_amplitude(double q) const
{
    const double qr =
        q * std::cbrt(
            (3.0 * excluded_volume) /
            (4.0 * M_PI)
        );

    if (std::abs(qr) < 1e-8)
        return excluded_volume;

    return excluded_volume *
           3.0 *
           (std::sin(qr) - qr * std::cos(qr)) /
           (qr * qr * qr);
}

bool FormFactorTable::load_from_string(
    const std::string& jsonContents)
{
    json j = json::parse(jsonContents);

    table_.clear();

    // reserve enough space for Z=1..118
    table_.resize(118);


    for (auto& [element, value] : j.items())
    {
        for (size_t i = 0; i < value.size(); i++)
        {
            if (value[i].is_null())
            {
                std::cerr
                    << "NULL value found for element "
                    << element
                    << " index "
                    << i
                    << std::endl;

                continue;
            }
        }
        if (value.size() != 14)
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
            value[5].get<double>();


        // b1-b5
        for (int i = 0; i < 5; i++)
        {
            ff.b[i] =
                value[6+i].is_number()
                ? value[6+i].get<double>()
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
                
        ff.excluded_volume =
            (4.0 / 3.0) *
            M_PI *
            std::pow(ff.atomic_radius, 3.0);

        // IMPORTANT:
        // your vector index is Z-1
        int index =
            ff.atomic_number - 1;


        if (index >= 0)
        {
            table_[index] = ff;
        }
    }


    std::cout
        << "Loaded "
        << table_.size()
        << " form factor entries."
        << std::endl;


    return true;
}



bool FormFactorTable::contains(
    uint32_t index) const
{
    return index < table_.size()
        && table_[index].atomic_number != 0;
}



const FormFactor& FormFactorTable::operator[](
    uint32_t index) const
{
    if (!contains(index))
    {
        throw std::runtime_error(
            "Missing form factor index");
    }


    return table_[index];
}



size_t FormFactorTable::size() const
{
    return table_.size();
}


std::string FormFactorTable::test_lookup(
    const std::vector<uint32_t>& atomicNumbers
) const
{
    std::ostringstream out;

    out << "Form factor lookup test\n";
    out << "======================\n\n";

    int count = 0;

    for (uint32_t atomicNumber : atomicNumbers)
    {
        if (count >= 10)
            break;


        uint32_t index =
            atomicNumber - 1;


        out << "Atom "
            << count
            << "\n";


        out << "  Atomic number: "
            << atomicNumber
            << "\n";


        out << "  JSON index: "
            << index
            << "\n";


        if (!contains(index))
        {
            out << "  ERROR: no form factor found\n\n";
            count++;
            continue;
        }


        const auto& ff =
            (*this)[index];


        out << "  Stored Z: "
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