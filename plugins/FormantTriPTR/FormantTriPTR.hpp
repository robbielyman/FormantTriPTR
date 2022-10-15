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
    void next_aaaka(int nSamples);
    void next_aaakk(int nSamples);
    void next_aakka(int nSamples);
    void next_aakkk(int nSamples);
    void next_akaka(int nSamples);
    void next_akakk(int nSamples);
    void next_akkka(int nSamples);
    void next_akkkk(int nSamples);
    void next_kaaka(int nSamples);
    void next_kaakk(int nSamples);
    void next_kakka(int nSamples);
    void next_kakkk(int nSamples);
    void next_kkaka(int nSamples);
    void next_kkakk(int nSamples);
    void next_kkkka(int nSamples);
    void next_kkkkk(int nSamples);
    double algorithm(double p, double t0, double t2, double t3,
            double w, double v, double a, double b, double c,
            double dc, double p1, double p2, double p3);
    double algorithm_retrig(double p, double t0, double t2, double t3,
            double w, double v, double a, double b, double c,
            double dc, double p1, double p2, double p3, double offset, double turnaround);

    // Member variables
    double mPhase   = 0;
    int mCounter = 0;
    double mSync    = 0;
    double mPhaseIn = 0;
    double mOffset = -1.f;
    double mLast = 0.f;
};

} // namespace FormantTriPTR
