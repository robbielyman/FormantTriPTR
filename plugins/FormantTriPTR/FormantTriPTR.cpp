// PluginFormantTriPTR.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "FormantTriPTR.hpp"
#include "math.h"

static InterfaceTable* ft;

namespace FormantTriPTR {

FormantTriPTR::FormantTriPTR() {
    if (isAudioRateIn(0)) {
        if (isAudioRateIn(1)) {
            if (isAudioRateIn(2)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aaaka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aaakk>();
                    }
            }
            else {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aakka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aakkk>();
                    }
            }
        }
        else {
            if (isAudioRateIn(2)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akaka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akakk>();
                    }
            }
            else {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akkka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akkkk>();
                    }
            }
        }
    }
    else {
        if (isAudioRateIn(1)) {
            if (isAudioRateIn(2)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kaaka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kaakk>();
                    }
            }
            else {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kakka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kakkk>();
                    }
            }
        }
        else {
            if (isAudioRateIn(2)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkaka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkakk>();
                    }
            }
            else {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkkka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkkkk>();
                    }
            }
        }
    }
    next_kkkkk(1);
}

inline double FormantTriPTR::algorithm_retrig(double p, double t0, double t2, double t3,
        double w, double v, double a, double b, double c, double dc,
        double p1, double p2, double p3, double offset, double turnaround) {
    if (p < turnaround) {
        if (p < t0) {
            // y = B*x - B*DC - 1 - 0.5*P3*x^4
            return b*p - b*dc - 0.5f * p3 * p * p * p * p + offset;
        }
        else if (p < t2) {
            // y = B*x - B*DC - 1 + P3*x^4 - 0.5*P2*x^3 + 0.75*P1*x^2 - 0.5*C*x + 0.125*C*T0
            return b*p - b*dc + p3 * p * p * p * p
                - 0.5f * p2 * p * p * p + 0.75f * p1 * p * p - 0.5 * c * p
                + 0.125f * c * t0 + offset;
        }
        else if (p < t3) {
            // y = B*x - B*DC - 1 - 0.5*P3*x^4 + 0.5*P2*x^3 - 2.25*P1*x^2 + 3.5*C*x - 1.875*C*T0
            return b*p - b*dc - 0.5f * p3 * p * p * p * p
                + 0.5f * p2 * p * p * p - 2.25f * p1 * p * p
                + 3.5f * c * p - 1.875f * c * t0 + offset;
        }
        else {
            // y = A*x - A*DC - 1
            return a * p - a * dc + offset;
        }
    }
    else {
        p = p - turnaround;
        if (p < t0) {
            // y = A*x - A*DC + 1 + 0.5*P3*x^4
            return a * p - a * dc + 1.f + 0.5f * p3 * p * p * p * p;
        }
        else if (p < t2) {
            // y = A*x - A*DC + 1 - P3*x^4 + 0.5*P2*x^3 - 0.75*P1*x^2 + 0.5*C*x - 0.125*C*T0
            return a * p - a * dc + 1.f - p3 * p * p * p * p
                + 0.5f * p2 * p * p * p - 0.75f * p1 * p * p + 0.5 * c * p
                - 0.125f * c * t0;
        }
        else if (p < t3) {
            // y = A*x - A*DC + 1 + 0.5*P3*x^4 - P2*x^3 + 4.5*P1*x^2 - 7*C*x + 3.75*C*T0
            return a * p - a * dc + 1.f + 0.5f * p3 * p * p * p * p
                - 0.5f * p2 * p * p * p + 2.25f * p1 * p * p
                - 3.5f * c * p + 1.875f * c * t0;
        }
        else if (p < v) {
            // y = B*x - B*DC + 1
            return b * p - b * dc + 1.f;
        }
        else {
            double p = p - v;
            if (p < t0) {
                // y = B*x - 1 - 1.5*B*T0 - 0.5*w*P3*x^4
                return b * p - 1.f - 1.5f * b * t0 - 0.5f * w * p3 * p * p * p * p;
            }
            else if (p < t2) {
                // y = 0.5*B*x - 1 - 1.375*B*T0 + 0.75*w*P1*x^2 - 0.5*w*P2*x^3 + w*P3*x^4
                return 0.5f * b * p - 1.f - 1.375f * b * t0 + 0.75 * w * p1 * p * p
                    - 0.5 * w * p2 * p * p * p + w * p3 * p * p * p * p;
            }
            else if (p < t3) {
                // y = 4.5*B*x - 1 - 3.375*B*T0 - 2.25*w*P1*x^2 + 0.5*w*P2*x^3 - 0.5*w*P3*x^4
                return 4.5f * b * p - 1.f - 3.375f * b * t0 - 2.25f * w * p1 * p * p
                    + 0.5f * w * p2 * p * p * p - 0.5f * w * p3 * p * p * p * p;
            }
            else {
                return -1.f;
            }
        }
    }
}

inline double FormantTriPTR::algorithm(double p, double t0, double t2, double t3,
        double w, double v, double a, double b, double c,
        double dc, double p1, double p2, double p3) {
    if (p < t0) {
        // y = 0.5*v*p3*x^4 - 1
        return 0.5 * v * p3 * p * p * p * p - 1.f;
    }
    else if (p < t2) {
        // y = v*p3*x^4 - 0.5*v*p2*x^3 + .75*v*p1*x^2 - .125*A*T0 + 0.5*A*x - 1
        return v * p3 * p * p * p * p - 0.5f * v * p2 * p * p * p
            + 0.75f * v * p1 * p * p - .125f * a * t0 + 0.5 * a * p - 1.f;
    }
    else if (p < t3) {
        // y = -0.5*v*p3*x^4 + 0.5*v*p2*x^3 - 2.25*v*p1*x^2 + 1.875*A*T0 - 3.5*A*x - 1
        return -0.5f * v * p3 * p * p * p * p + 0.5f * v * p2 * p * p * p
            - 2.25f * v * p1 * p * p + 1.875f * a * t0 - 3.5f * a * p - 1.f;
    }
    else if (p < w) {
        return a * p - a * dc - 1.f;
    }
    else {
        p = p - w;
        if (p < t0) {
            // y = A*x - A*DC + 1 + 0.5*P3*x^4
            return a * p - a * dc + 1.f + 0.5f * p3 * p * p * p * p;
        }
        else if (p < t2) {
            // y = A*x - A*DC + 1 - P3*x^4 + 0.5*P2*x^3 - 0.75*P1*x^2 + 0.5*C*x - 0.125*C*T0
            return a * p - a * dc + 1.f - p3 * p * p * p * p
                + 0.5f * p2 * p * p * p - 0.75f * p1 * p * p + 0.5 * c * p
                - 0.125f * c * t0;
        }
        else if (p < t3) {
            // y = A*x - A*DC + 1 + 0.5*P3*x^4 - P2*x^3 + 4.5*P1*x^2 - 7*C*x + 3.75*C*T0
            return a * p - a * dc + 1.f + 0.5f * p3 * p * p * p * p
                - 0.5f * p2 * p * p * p + 2.25f * p1 * p * p
                - 3.5f * c * p + 1.875f * c * t0;
        }
        else if (p < v) {
            // y = B*x - B*DC + 1
            return b * p - b * dc + 1.f;
        }
        else {
            p = p - v;
            if (p < t0) {
                // y = B*x - 1 - 1.5*B*T0 - 0.5*w*P3*x^4
                return b * p - 1.f - 1.5f * b * t0 - 0.5f * w * p3 * p * p * p * p;
            }
            else if (p < t2) {
                // y = 0.5*B*x - 1 - 1.375*B*T0 + 0.75*w*P1*x^2 - 0.5*w*P2*x^3 + w*P3*x^4
                return 0.5f * b * p - 1.f - 1.375f * b * t0 + 0.75 * w * p1 * p * p
                    - 0.5 * w * p2 * p * p * p + w * p3 * p * p * p * p;
            }
            else if (p < t3) {
                // y = 4.5*B*x - 1 - 3.375*B*T0 - 2.25*w*P1*x^2 + 0.5*w*P2*x^3 - 0.5*w*P3*x^4
                return 4.5f * b * p - 1.f - 3.375f * b * t0 - 2.25f * w * p1 * p * p
                    + 0.5f * w * p2 * p * p * p - 0.5f * w * p3 * p * p * p * p;
            }
            else {
                return -1.f;
            }
        }
    }
}

void FormantTriPTR::next_aaaka(int nSamples) {
    const float* freq    = in(0);
    const float* formant = in(1);
    const float* width   = in(2);
    const float phase   = in0(3);
    const float* sync    = in(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.0001f : formant[i];
        double ratio    = freq[i] / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync[i] > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync[i];
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_aaakk(int nSamples) {
    const float* freq    = in(0);
    const float* formant = in(1);
    const float* width   = in(2);
    const float phase   = in0(3);
    const float sync    = in0(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.0001f : formant[i];
        double ratio    = freq[i] / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync;
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_aakka(int nSamples) {
    const float* freq    = in(0);
    const float* formant = in(1);
    const float width   = in0(2);
    const float phase   = in0(3);
    const float* sync    = in(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.0001f : formant[i];
        double ratio    = freq[i] / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync[i] > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync[i];
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_aakkk(int nSamples) {
    const float* freq    = in(0);
    const float* formant = in(1);
    const float width   = in0(2);
    const float phase   = in0(3);
    const float sync    = in0(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.0001f : formant[i];
        double ratio    = freq[i] / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync;
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_akaka(int nSamples) {
    const float* freq    = in(0);
    const float formant = in0(1);
    const float* width   = in(2);
    const float phase   = in0(3);
    const float* sync    = in(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant == 0.f ? 0.0001f : formant;
        double ratio    = freq[i] / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync[i] > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync[i];
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_akakk(int nSamples) {
    const float* freq    = in(0);
    const float formant = in0(1);
    const float* width   = in(2);
    const float phase   = in0(3);
    const float sync    = in0(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant == 0.f ? 0.0001f : formant;
        double ratio    = freq[i] / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync;
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_akkka(int nSamples) {
    const float* freq    = in(0);
    const float formant = in0(1);
    const float width   = in0(2);
    const float phase   = in0(3);
    const float* sync    = in(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant == 0.f ? 0.0001f : formant;
        double ratio    = freq[i] / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync[i] > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync[i];
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_akkkk(int nSamples) {
    const float* freq    = in(0);
    const float formant = in0(1);
    const float width   = in0(2);
    const float phase   = in0(3);
    const float sync    = in0(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq[i] * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant == 0.f ? 0.0001f : formant;
        double ratio    = freq[i] / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync;
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_kaaka(int nSamples) {
    const float freq    = in0(0);
    const float* formant = in(1);
    const float* width   = in(2);
    const float phase   = in0(3);
    const float* sync    = in(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.0001f : formant[i];
        double ratio    = freq / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync[i] > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync[i];
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_kaakk(int nSamples) {
    const float freq    = in0(0);
    const float* formant = in(1);
    const float* width   = in(2);
    const float phase   = in0(3);
    const float sync    = in0(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.0001f : formant[i];
        double ratio    = freq / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync;
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_kakka(int nSamples) {
    const float freq    = in0(0);
    const float* formant = in(1);
    const float width   = in0(2);
    const float phase   = in0(3);
    const float* sync    = in(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.0001f : formant[i];
        double ratio    = freq / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync[i] > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync[i];
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_kakkk(int nSamples) {
    const float freq    = in0(0);
    const float* formant = in(1);
    const float width   = in0(2);
    const float phase   = in0(3);
    const float sync    = in0(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.0001f : formant[i];
        double ratio    = freq / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync;
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_kkaka(int nSamples) {
    const float freq    = in0(0);
    const float formant = in0(1);
    const float* width   = in(2);
    const float phase   = in0(3);
    const float* sync    = in(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant == 0.f ? 0.0001f : formant;
        double ratio    = freq / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync[i] > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync[i];
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_kkakk(int nSamples) {
    const float freq    = in0(0);
    const float formant = in0(1);
    const float* width   = in(2);
    const float phase   = in0(3);
    const float sync    = in0(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        double step     = freq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant == 0.f ? 0.0001f : formant;
        double ratio    = freq / form;
        ratio           = ratio < 0.f ? -ratio : ratio;
        double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
        w               = w * ratio;
        double v        = ratio - w;
        double a        = 2.f / w;
        double b        = -2.f / v;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        if (sync > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync;
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_kkkka(int nSamples) {
    const float freq    = in0(0);
    const float formant = in0(1);
    const float width   = in0(2);
    const float phase   = in0(3);
    const float* sync    = in(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double step     = freq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.0001f : formant;
    double ratio    = freq / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
    w               = w * ratio;
    double v        = ratio - w;
    double a        = 2.f / w;
    double b        = -2.f / v;
    double c        = 0.5f * a * b;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        if (sync[i] > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter > turnaround && pos <= step) {
            if (counter > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync[i];
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

void FormantTriPTR::next_kkkkk(int nSamples) {
    const float freq    = in0(0);
    const float formant = in0(1);
    const float width   = in0(2);
    const float phase   = in0(3);
    const float sync    = in0(4);
    float* outbuf       = out(0);

    double pos      = mPhase;
    int counter     = mCounter;
    double lastsync = mSync;
    double phasein  = mPhaseIn;
    double last     = mLast;
    double offset   = mOffset;
    double rateinv  = 1 / sampleRate();
    double step     = freq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.0001f : formant;
    double ratio    = freq / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
    w               = w * ratio;
    double v        = ratio - w;
    double a        = 2.f / w;
    double b        = -2.f / v;
    double c        = 0.5f * a * b;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    double out1     = 0.f;
    for (int i=0; i < nSamples; ++i) {
        if (sync > 0.f && lastsync <= 0.f) {
            pos = 0.f;
            counter = 0;
            offset = -1.f;
        }
        else {
            pos = pos + phase - phasein + step;
        }
        while (pos > 1.f) {
            pos -= 1.f;
            counter++;
        }
        while (pos < 0.f) {
            pos += 1.f;
            counter--;
        }
        while (counter < 0) {
            counter += ceil(ratio);
        }
        double turnaround = (1.f - offset) * 0.5 * w;
        if (counter + pos > turnaround && pos <= step) {
            if (counter + pos > ratio) {
                counter = 0;
                offset = -1.f;
            }
            else {
                offset = last;
                counter = 0;
                turnaround = (1.f - offset) * 0.5 * w;
            }

        }
        if (offset <= -1.f) {
            out1 = algorithm(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
        }
        else {
            out1 = algorithm_retrig(pos + counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3, offset, turnaround);
        }
        outbuf[i]   = out1;
        phasein     = phase;
        lastsync    = sync;
        last        = out1;
    }
    mPhase = pos;
    mCounter = counter;
    mPhaseIn = phasein;
    mSync = lastsync;
    mOffset = offset;
    mLast = last;
}

} // namespace FormantTriPTR

PluginLoad(FormantTriPTRUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<FormantTriPTR::FormantTriPTR>(ft, "FormantTriPTR", false);
}
