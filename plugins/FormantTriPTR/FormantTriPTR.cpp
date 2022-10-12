// PluginFormantTriPTR.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "FormantTriPTR.hpp"

static InterfaceTable* ft;

namespace FormantTriPTR {

FormantTriPTR::FormantTriPTR() {
    mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next>();
    next(1);
}

void FormantTriPTR::next(int nSamples) {
    const float* input = in(0);
    const float* gain = in(1);
    float* outbuf = out(0);

    // simple gain function
    for (int i = 0; i < nSamples; ++i) {
        outbuf[i] = input[i] * gain[i];
    }
}

} // namespace FormantTriPTR

PluginLoad(FormantTriPTRUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<FormantTriPTR::FormantTriPTR>(ft, "FormantTriPTR", false);
}
