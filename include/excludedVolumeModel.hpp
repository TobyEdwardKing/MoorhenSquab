#ifndef EXCLUDEDVOLUMEMODEL_HPP
#define EXCLUDEDVOLUMEMODEL_HPP

#include <vector>
#include <complex>
#include <cstdint>

#include "formFactor.hpp"
#include "atomExtractor.hpp"


struct ComplexCoefficient
{
    std::complex<double> value{0.0, 0.0};
};


class ExcludedVolumeModel
{
public:

    ExcludedVolumeModel() = default;

    //
    // Calculate reference excluded-volume
    // multipole coefficients C_lm(q).
    //
    void compute_coefficients(
        const std::vector<AtomData>& atoms,
        const FormFactorTable& table,
        const std::vector<double>& qGrid,
        int L,
        double radiusAdjustmentFraction = 0.0
    );
    

    //
    // Access calculated coefficients.
    //
    const std::vector<std::vector<std::complex<double>>>&
    coefficients() const;


private:

    int L_ = 0;

    //
    // First index:
    //     q point
    //
    // Second index:
    //     flattened (l,m)
    //
    std::vector<
        std::vector<std::complex<double>>
    > coefficients_;

    //
    // Convert (l,m) into a flat array index.
    //
    static int lm_index(
        int l,
        int m
    );

    //
    // Calculate the Gaussian dummy-atom
    // form factor for one excluded volume.
    //
    static double dummy_form_factor(
        double q,
        double volume
    );
};

#endif