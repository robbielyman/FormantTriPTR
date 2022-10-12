// PluginFormantTriPTR.hpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#pragma once

#include "SC_PlugIn.hpp"

namespace FormantTriPTR {

class FormantTriPTR : public SCUnit {
public:
    FormantTriPTR();

    // Destructor
    // ~FormantTriPTR();

private:
    // Calc function
    void next(int nSamples);

    // Member variables
};

} // namespace FormantTriPTR
