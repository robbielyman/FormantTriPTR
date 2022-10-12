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
    void next_aaaaa(int nSamples);
    void next_aaaak(int nSamples);
    void next_aaaka(int nSamples);
    void next_aaakk(int nSamples);
    void next_aakaa(int nSamples);
    void next_aakak(int nSamples);
    void next_aakka(int nSamples);
    void next_aakkk(int nSamples);
    void next_akaaa(int nSamples);
    void next_akaak(int nSamples);
    void next_akaka(int nSamples);
    void next_akakk(int nSamples);
    void next_akkaa(int nSamples);
    void next_akkak(int nSamples);
    void next_akkka(int nSamples);
    void next_akkkk(int nSamples);
    void next_kaaaa(int nSamples);
    void next_kaaak(int nSamples);
    void next_kaaka(int nSamples);
    void next_kaakk(int nSamples);
    void next_kakaa(int nSamples);
    void next_kakak(int nSamples);
    void next_kakka(int nSamples);
    void next_kakkk(int nSamples);
    void next_kkaaa(int nSamples);
    void next_kkaak(int nSamples);
    void next_kkaka(int nSamples);
    void next_kkakk(int nSamples);
    void next_kkkaa(int nSamples);
    void next_kkkak(int nSamples);
    void next_kkkka(int nSamples);
    void next_kkkkk(int nSamples);
    double rise_prime(double p, double t0, double t2, double t3,
            double w, double v, double a, double b, double c,
            double dc, double p1, double p2, double p3);
    double rise(double p, double t0, double t2, double t3,
            double w, double a, double b, double c,
            double dc, double p1, double p2, double p3);
    double fall(double p, double t0, double t2, double t3,
            double w, double v, double a, double b, double c,
            double dc, double p1, double p2, double p3);

    // Member variables
    double mPhase   = 0;
    double mCounter = 0;
    double mSync    = 0;
    double mPhaseIn = 0;
    double mLast    = -1;
    double mOffset  = -1;
    // double mTurnaround = -1;
};

} // namespace FormantTriPTR
