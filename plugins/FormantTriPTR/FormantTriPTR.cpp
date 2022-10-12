// PluginFormantTriPTR.cpp
// Rylee Alanza Lyman (ryleealanza@gmail.com)

#include "SC_PlugIn.hpp"
#include "FormantTriPTR.hpp"

static InterfaceTable* ft;

namespace FormantTriPTR {

FormantTriPTR::FormantTriPTR() {
    if (isAudioRateIn(0)) {
        if (isAudioRateIn(1)) {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aaaaa>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aaaak>();
                    }
                }
                else {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aaaka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aaakk>();
                    }
                }
            }
            else {
                if (isAudioRateIn(3)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aakaa>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aakak>();
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
        }
        else {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akaaa>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akaak>();
                    }
                }
                else {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akaka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akakk>();
                    }
                }
            }
            else {
                if (isAudioRateIn(3)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akkaa>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akkak>();
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
    }
    else {
        if (isAudioRateIn(1)) {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kaaaa>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kaaak>();
                    }
                }
                else {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kaaka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kaakk>();
                    }
                }
            }
            else {
                if (isAudioRateIn(3)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kakaa>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kakak>();
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
        }
        else {
            if (isAudioRateIn(2)) {
                if (isAudioRateIn(3)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkaaa>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkaak>();
                    }
                }
                else {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkaka>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkakk>();
                    }
                }
            }
            else {
                if (isAudioRateIn(3)) {
                    if (isAudioRateIn(4)) {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkkaa>();
                    }
                    else {
                        mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkkak>();
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
    }
    next_kkkkk(1);
}

inline double FormantTriPTR::rise_prime(double p, double t0, double t2, double t3,
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
    else {
        return a * p - a * dc - 1.f;
    }
    /* else {
        double p = p - w;
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
            double p = p - 1.f + w;
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
    */
}

inline double FormantTriPTR::rise(double p, double t0, double t2, double t3,
        double w, double a, double b, double c,
        double dc, double p1, double p2, double p3) {
    if (p < t0) {
        // y = B*x - B*DC - 1 - 0.5*P3*x^4
        return b * p - b * dc - 1.f - 0.5f * p3 * p * p * p * p;
    }
    else if (p < t2) {
        // y = B*x - B*DC - 1 + P3*x^4 - 0.5*P2*x^3 + 0.75*P1*x^2 - 0.5*C*x + 0.125*C*T0
        return b * p - b * dc - 1.f + p3 * p * p * p * p
            - 0.5f * p2 * p * p * p + 0.75f * p1 * p * p - 0.5f * c * p
            + 0.125f * c * t0;
    }
    else if (p < t3) {
        // y = B*x - B*DC - 1 - 0.5*P3*x^4 + 0.5*P2*x^3 - 2.25*P1*x^2 + 3.5*C*x - 1.875*C*T0
        return b * p - b * dc - 1.f - 0.5f * p3 * p * p * p * p
            + 0.5f * p2 * p * p * p - 2.25f * p1 * p * p
            + 3.5f * c * p - 1.875f * c * t0;
    }
    else {
        return a * p - a * dc - 1.f;
    }
}
    
inline double FormantTriPTR::fall(double p, double t0, double t2, double t3,
        double w, double v, double a, double b, double c,
        double dc, double p1, double p2, double p3) {
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

void FormantTriPTR::next_aaaaa(int nSamples) {
    const float* freq   = in(0);
    const float* formant = in(1);
    const float* width  = in(2);
    const float* phase  = in(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_aaaak(int nSamples) {
    const float* freq   = in(0);
    const float* formant = in(1);
    const float* width  = in(2);
    const float* phase  = in(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_aaaka(int nSamples) {
    const float* freq   = in(0);
    const float* formant = in(1);
    const float* width  = in(2);
    const float phase  = in0(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_aaakk(int nSamples) {
    const float* freq   = in(0);
    const float* formant = in(1);
    const float* width  = in(2);
    const float phase  = in0(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_aakaa(int nSamples) {
    const float* freq   = in(0);
    const float* formant = in(1);
    const float width  = in0(2);
    const float* phase  = in(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_aakak(int nSamples) {
    const float* freq   = in(0);
    const float* formant = in(1);
    const float width  = in0(2);
    const float* phase  = in(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_aakka(int nSamples) {
    const float* freq   = in(0);
    const float* formant = in(1);
    const float width  = in0(2);
    const float phase  = in0(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_aakkk(int nSamples) {
    const float* freq   = in(0);
    const float* formant = in(1);
    const float width  = in0(2);
    const float phase  = in0(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_akaaa(int nSamples) {
    const float* freq   = in(0);
    const float formant = in0(1);
    const float* width  = in(2);
    const float* phase  = in(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_akaak(int nSamples) {
    const float* freq   = in(0);
    const float formant = in0(1);
    const float* width  = in(2);
    const float* phase  = in(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double form     = formant == 0.f ? 0.01f : formant;
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_akaka(int nSamples) {
    const float* freq   = in(0);
    const float formant = in0(1);
    const float* width  = in(2);
    const float phase  = in0(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_akakk(int nSamples) {
    const float* freq   = in(0);
    const float formant = in0(1);
    const float* width  = in(2);
    const float phase  = in0(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_akkaa(int nSamples) {
    const float* freq   = in(0);
    const float formant = in0(1);
    const float width  = in0(2);
    const float* phase  = in(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_akkak(int nSamples) {
    const float* freq   = in(0);
    const float formant = in0(1);
    const float width  = in0(2);
    const float* phase  = in(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_akkka(int nSamples) {
    const float* freq   = in(0);
    const float formant = in0(1);
    const float width  = in0(2);
    const float phase  = in0(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_akkkk(int nSamples) {
    const float* freq   = in(0);
    const float formant = in0(1);
    const float width  = in0(2);
    const float phase  = in0(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    for(int i = 0; i < nSamples; ++i) {
        double posfreq  = freq[i] < 0.f ? -1.f * freq[i] : freq[i];
        double step     = posfreq * rateinv;
        double step2    = 2.f * step;
        double step3    = 3.f * step;
        double samples  = 1.f / step;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kaaaa(int nSamples) {
    const float freq   = in0(0);
    const float* formant = in(1);
    const float* width  = in(2);
    const float* phase  = in(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    for(int i = 0; i < nSamples; ++i) {
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kaaak(int nSamples) {
    const float freq   = in0(0);
    const float* formant = in(1);
    const float* width  = in(2);
    const float* phase  = in(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    for(int i = 0; i < nSamples; ++i) {
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kaaka(int nSamples) {
    const float freq   = in0(0);
    const float* formant = in(1);
    const float* width  = in(2);
    const float phase  = in0(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    for(int i = 0; i < nSamples; ++i) {
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kaakk(int nSamples) {
    const float freq   = in0(0);
    const float* formant = in(1);
    const float* width  = in(2);
    const float phase  = in0(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    for(int i = 0; i < nSamples; ++i) {
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kakaa(int nSamples) {
    const float freq   = in0(0);
    const float* formant = in(1);
    const float width  = in0(2);
    const float* phase  = in(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    for(int i = 0; i < nSamples; ++i) {
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kakak(int nSamples) {
    const float freq   = in0(0);
    const float* formant = in(1);
    const float width  = in0(2);
    const float* phase  = in(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    for(int i = 0; i < nSamples; ++i) {
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kakka(int nSamples) {
    const float freq   = in0(0);
    const float* formant = in(1);
    const float width  = in0(2);
    const float phase  = in0(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    for(int i = 0; i < nSamples; ++i) {
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kakkk(int nSamples) {
    const float freq   = in0(0);
    const float* formant = in(1);
    const float width  = in0(2);
    const float phase  = in0(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    for(int i = 0; i < nSamples; ++i) {
        double form     = formant[i] == 0.f ? 0.01f : formant[i];
        form = form < 0.f ? -1.f * form : form;
        double ratio    = posfreq / form;
        double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kkaaa(int nSamples) {
    const float freq   = in0(0);
    const float formant = in0(1);
    const float* width  = in(2);
    const float* phase  = in(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    double ratio    = posfreq / form;
    for(int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kkaak(int nSamples) {
    const float freq   = in0(0);
    const float formant = in0(1);
    const float* width  = in(2);
    const float* phase  = in(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    double ratio    = posfreq / form;
    for(int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kkaka(int nSamples) {
    const float freq   = in0(0);
    const float formant = in0(1);
    const float* width  = in(2);
    const float phase  = in0(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    double ratio    = posfreq / form;
    for(int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kkakk(int nSamples) {
    const float freq   = in0(0);
    const float formant = in0(1);
    const float* width  = in(2);
    const float phase  = in0(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    double ratio    = posfreq / form;
    for(int i = 0; i < nSamples; ++i) {
        double w        = width[i] < 0.01f ? 0.01f : width[i] > 0.99f ? 0.99f : width[i];
        w = w * ratio;
        double a        = 2.f / w;
        double b        = 2.f / (w - ratio);
        double v        = ratio - w;
        double c        = 0.5f * a * b;
        double dc       = 1.5f * step;
        double p1       = c * samples;
        double p2       = p1 * samples;
        double p3       = p2 * samples / 12.f;
        double turnaround = (1.f - offset - b * dc) * 0.5 * w;
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kkkaa(int nSamples) {
    const float freq   = in0(0);
    const float formant = in0(1);
    const float width  = in0(2);
    const float* phase  = in(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    double ratio    = posfreq / form;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    w = w * ratio;
    double a        = 2.f / w;
    double b        = 2.f / (w - ratio);
    double v        = ratio - w;
    double c        = 0.5f * a * b;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    double turnaround = (1.f - offset - b * dc) * 0.5 * w;
    for(int i = 0; i < nSamples; ++i) {
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kkkak(int nSamples) {
    const float freq   = in0(0);
    const float formant = in0(1);
    const float width  = in0(2);
    const float* phase  = in(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    double ratio    = posfreq / form;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    w = w * ratio;
    double a        = 2.f / w;
    double b        = 2.f / (w - ratio);
    double v        = ratio - w;
    double c        = 0.5f * a * b;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    double turnaround = (1.f - offset - b * dc) * 0.5 * w;
    for(int i = 0; i < nSamples; ++i) {
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase[i] - phasein) + step;
            counter = counter + (phase[i] - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase[i];
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kkkka(int nSamples) {
    const float freq   = in0(0);
    const float formant = in0(1);
    const float width  = in0(2);
    const float phase  = in0(3);
    const float* sync   = in(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    double ratio    = posfreq / form;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    w = w * ratio;
    double a        = 2.f / w;
    double b        = 2.f / (w - ratio);
    double v        = ratio - w;
    double c        = 0.5f * a * b;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    double turnaround = (1.f - offset - b * dc) * 0.5 * w;
    for(int i = 0; i < nSamples; ++i) {
        double out1     = 0.f;
        if (sync[i] > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync[i];
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

void FormantTriPTR::next_kkkkk(int nSamples) {
    const float freq   = in0(0);
    const float formant = in0(1);
    const float width  = in0(2);
    const float phase  = in0(3);
    const float sync   = in0(4);
    float* outbuf = out(0);

    double pos = mPhase;
    double counter  = mCounter;
    double phasein  = mPhaseIn;
    double lastsync = mSync;
    double rateinv  = 1 / sampleRate();
    double lastval  = mLast;
    double offset   = mOffset;
    double posfreq  = freq < 0.f ? -1.f * freq : freq;
    double step     = posfreq * rateinv;
    double step2    = 2.f * step;
    double step3    = 3.f * step;
    double samples  = 1.f / step;
    double form     = formant == 0.f ? 0.01f : formant;
    form = form < 0.f ? -1.f * form : form;
    double ratio    = posfreq / form;
    double w        = width < 0.01f ? 0.01f : width > 0.99f ? 0.99f : width;
    w = w * ratio;
    double a        = 2.f / w;
    double b        = 2.f / (w - ratio);
    double v        = ratio - w;
    double c        = 0.5f * a * b;
    double dc       = 1.5f * step;
    double p1       = c * samples;
    double p2       = p1 * samples;
    double p3       = p2 * samples / 12.f;
    double turnaround = (1.f - offset - b * dc) * 0.5 * w;
    for(int i = 0; i < nSamples; ++i) {
        double out1     = 0.f;
        if (sync > 0.f && lastsync == 0.f) {
            pos = 0.f;
            counter = pos;
        }
        else {
            pos = pos + (phase - phasein) + step;
            counter = counter + (phase - phasein) + step;
            while (pos > 1.f) {
                pos -= 1.f;
            }
            while (pos < 0.f) {
                pos += 1.f;
            }
            if (counter <= turnaround) {
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // retrigger while falling
            else if (pos <= step) {
                if (ratio > 1.f) {
                    offset = lastval + b * dc;
                }
                else {
                    offset = lastval;
                }
                turnaround = (1.f - offset - b * dc) * 0.5 * w;
                counter = step;
                if (ratio > 1.f) {
                    out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
                }
                else {
                    out1 = offset + 1 + rise_prime(counter, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
                }
            }
            // should we start rising again?
            else if (ratio > 1.0f && counter > turnaround + v) {
                while (counter > turnaround + v) {
                    counter -= turnaround + v;
                }
                offset = -1.f + b * dc;
                turnaround = w;
                out1 = offset + 1 + rise(counter, step, step2, step3, w, a, b, c, dc, p1, p2, p3);
            }
            // falling
            else {
                out1 = fall(counter - turnaround, step, step2, step3, w, v, a, b, c, dc, p1, p2, p3);
            }
        }
        outbuf[i] = out1;
        phasein   = phase;
        lastsync  = sync;
        lastval   = out1;
    }
    mPhase  = pos;
    mSync   = lastsync;
    mPhaseIn = phasein;
    mCounter = counter;
    mOffset = offset;
    mLast   = lastval;
    // mTurnaround = turnaround;
}

} // namespace FormantTriPTR

PluginLoad(FormantTriPTRUGens) {
    // Plugin magic
    ft = inTable;
    registerUnit<FormantTriPTR::FormantTriPTR>(ft, "FormantTriPTR", false);
}
