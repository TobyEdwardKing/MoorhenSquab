#pragma once

#include <vector>

#include "atomExtractor.hpp"
#include "formFactor.hpp"


class ApplicationState
{

public:

    std::vector<AtomData> atoms;

    FormFactorTable formFactors;


};  