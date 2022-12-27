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
        if (isAudioRateIn(4)) mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aaaka>();
        else                  mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aaakk>();
      } else {
        if (isAudioRateIn(4)) mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aakka>();
        else                  mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_aakkk>();
      }
    } else {
      if (isAudioRateIn(2)) {
        if (isAudioRateIn(4)) mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akaka>();
        else                  mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akakk>();

      } else {
        if (isAudioRateIn(4)) mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akkka>();
        else                  mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_akkkk>();

      }
    }
  } else {
    if (isAudioRateIn(1)) {
      if (isAudioRateIn(2)) {
        if (isAudioRateIn(4)) mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kaaka>();
        else                  mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kaakk>();
      } else {
        if (isAudioRateIn(4)) mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kakka>();
        else                  mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kakkk>();
      }
    } else {
      if (isAudioRateIn(2)) {
        if (isAudioRateIn(4)) mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkaka>();
        else                  mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkakk>();

      } else {
        if (isAudioRateIn(4)) mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkkka>();
        else                  mCalcFunc = make_calc_function<FormantTriPTR, &FormantTriPTR::next_kkkkk>();
      }
    }
  }
  next_kkkkk(1);
}

inline double FormantTriPTR::algorithm_retrig(double p, double h, double w, double v, double offset, double turnaround) {
  double p0 = 1.f / h;
  double p2 = p0 * p0;
  double p3 = p0 * p0 * p0 / 12.f;
  double a = 2.f / w;
  double b = -2.f / v;
  if (p < turnaround) {
    if (p < h) {
      // y = (a-b)x^4/24h^3 - 1 + bx - 3bh/2
      return 0.5f * p3 * (a - b) * p * p * p * p + offset + b * p - 1.5f * b * h;
    } else if (p < 2*h) {
      // y = -(a-b)x^4/12h^3 - 1 + (a-b)x^3/2h^2 - 3(a-b)x^2/4h + (a+b)x/2 - ah/8 - 11bh/8
      return -p3 * (a - b) * p * p * p * p + offset + 0.5f * p2 * (a - b) * p * p * p - 0.75f * p0 * (a - b) * p * p + 0.5f * (a + b) * p - .125f * a * h - 1.375f * b * h; 
    } else if (p < 3*h) {
      // y = (a-b)x^4/24h^3 - 1 - (a-b)x^3/2h^2 + 9(a-b)x^2/4h - 7ax/2 + 9bx/2 + 15ah/8 - 27bh/8
      return 0.5f * p3 * (a - b) * p * p * p * p + offset - 0.5f * p2 * (a - b) * p * p * p + 2.25f * p0 * (a - b) * p * p - 3.5f * a * p + 4.5f * b * p + 1.875f * a * h - 3.375f * b * h;
    } else {
      // y = ax - 1 - 3ah/2 
      return a * p + offset - 1.5f * a * h; 
    }
  } else {
    p = p - turnaround;
    if (p < h) {
      // y = -(a-b)x^4/24h^3 + 1 + ax - 3ah/2
      return -0.5f * p3 * (a - b) * p * p * p * p + 1.f + a * p - 1.5f * a * h;
    } else if (p < 2*h) {
      // y = (a-b)x^4/12h^3 + 1 - (a-b)x^3/2h^2 + 3(a-b)x^2/4h + (a+b)x/2 - 11ah/8 - bh/8
      return p3 * (a - b) * p * p * p * p + 1.f - 0.5f * p2 * (a - b) * p * p * p + 0.75f * p0 * (a - b) * p * p + 0.5f * (a + b) * p - 1.375f * a * h - .125 * b * h;
    } else if (p < 3*h) {
      // y = -(a-b)x^4/24h^3 + 1 + (a-b)x^3/2h^2 - 9(a-b)x^2/4h + 9ax/2 - 7bx/2 - 27ah/8 + 15bh/8
      return -0.5f * p3 * (a - b) * p * p * p * p + 1.f + 0.5f * p2 * (a - b) * p * p * p - 2.25f * p0 * (a - b) * p * p + 4.5f * a * p - 3.5f * b * p - 3.375f * a * h + 1.875f * b * h;
    } else if (p < v) {
      // y = bx + 1 - 3bh/2
      return b * p + 1 - 1.5f * b * h;
    } else {
      p = p - v;
      if (p < h) {
        // y = -bx^4/24h^3 - 1 + bx - 3bh/2
        return - 0.5f * p3 * b * p * p * p * p - 1.f + b * p - 1.5f * b * h;
      } else if (p < 2*h) {
        // y = bx^4/12h^3 - 1 - bx^3/2h^2 + 3bx^2/4h + bx/2 - 11bh/8
        return p3 * b * p * p * p * p - 1.f - 0.5f * p2 * b * p * p * p + 0.75f * p0 * b * p * p + 0.5f * b * p - 1.375f * b * h;
      } else if (p < 3*h) {
        // y = -bx^4/24h^3 - 1 + bx^3/2h^2 - 9bx^2/4h + 9bx/2 - 27bh/8
        return -0.5f * p2 * b * p * p * p * p  - 1.f + 0.5f * p2 * b * p * p * p - 2.25f * p0 * b * p * p + 4.5f * b * p - 3.375f * b * h; 
      } else {
        return -1.f; 
      }
    }
  }
}

inline double FormantTriPTR::algorithm(double p, double h, double w, double v) {
  double p0 = 1.f / h;
  double p2 = p0 * p0;
  double p3 = p0 * p0 * p0 / 12.f;
  double a = 2.f / w;
  double b = -2.f / v;
  if (p < h) {
    // y = ax^4/24h^3 - 1
    return 0.5f * p3 * a * p * p * p * p - 1.f;
  } else if (p < 2*h) {
    // y = -ax^4/12h^3 - 1 + ax^3/2h^2 - 3ax^2/4h + ax/2 - ah/8
    return -p3 * a * p * p * p * p - 1.f + 0.5f * p2 * a * p * p * p - 0.75f * p0 * a * p * p + 0.5f * a * p - 0.125f * a * h;
  } else if (p < 3*h) {
    // y = ax^4/24h^3 - 1 - ax^3/2h^2 + 9ax^2/4h - 7ax/2 + 15ah/8
    return 0.5f * p3 * a * p * p * p * p - 1.f - 0.5f * p2 * a * p * p * p + 2.25f * p0 * a * p * p - 3.5f * a * p + 1.875f * a * h;
  } else if (p < w) {
    // y = ax - 1 - 3ah/2 
    return a * p - 1.f - 1.5f * a * h;
  } else {
    p = p - w;
    if (p < h) {
      // y = -(a-b)x^4/24h^3 + 1 + ax - 3ah/2
      return -0.5f * p3 * (a - b) * p * p * p * p + 1.f + a * p - 1.5f * a * h;
    } else if (p < 2*h) {
      // y = (a-b)x^4/12h^3 + 1 - (a-b)x^3/2h^2 + 3(a-b)x^2/4h + (a+b)x/2 - 11ah/8 - bh/8
      return p3 * (a - b) * p * p * p * p + 1.f - 0.5f * p2 * (a - b) * p * p * p + 0.75f * p0 * (a - b) * p * p + 0.5f * (a + b) * p - 1.375f * a * h - .125f * b * h;
    } else if (p < 3*h) {
      // y = -(a-b)x^4/24h^3 + 1 + (a-b)x^3/2h^2 - 9(a-b)x^2/4h + 9ax/2 - 7bx/2 - 27ah/8 + 15bh/8
      return -0.5f * p3 * (a - b) * p * p * p * p + 1.f + 0.5f * p2 * (a - b) * p * p * p - 2.25f * p0 * (a - b) * p * p + 4.5f * a * p - 3.5f * b * p - 3.375f * a * h + 1.875 * b * h;
    } else if (p < v) {
      // y = bx + 1 - 3bh/2
      return b * p + 1 - 1.5f * b * h;
    } else {
      p = p - v;
      if (p < h) {
        // y = -bx^4/24h^3 - 1 + bx - 3bh/2
        return -0.5f * p3 * b * p * p * p * p - 1.f + b * p - 1.5f * b * h;
      } else if (p < 2*h) {
        // y = bx^4/12h^3 - 1 - bx^3/2h^2 + 3bx^2/4h + bx/2 - 11bh/8
        return p3 * b * p * p * p * p - 1.f - 0.5f * p2 * b * p * p * p + 0.75f * p0 * b * p * p + 0.5f * b * p - 1.375f * b * h;
      } else if (p < 3*h) {
        // y = -bx^4/24h^3 - 1 + bx^3/2h^2 - 9bx^2/4h + 9bx/2 - 27bh/8
        return -0.5f * p3 * b * p * p * p * p - 1.f + 0.5f * p2 * b * p * p * p - 2.25f * p0 * b * p * p + 4.5 * b * p - 3.375 * b * h;
      } else {
        // y = -1
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
    double form     = formant[i] == 0.f ? 0.0001f : formant[i];
    double ratio    = freq[i] / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
    w               = w * ratio;
    double v        = ratio - w;
    if (sync[i] > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant[i] == 0.f ? 0.0001f : formant[i];
    double ratio    = freq[i] / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
    w               = w * ratio;
    double v        = ratio - w;
    if (sync > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant[i] == 0.f ? 0.0001f : formant[i];
    double ratio    = freq[i] / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
    w               = w * ratio;
    double v        = ratio - w;
    if (sync[i] > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant[i] == 0.f ? 0.0001f : formant[i];
    double ratio    = freq[i] / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
    w               = w * ratio;
    double v        = ratio - w;
    if (sync > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant == 0.f ? 0.0001f : formant;
    double ratio    = freq[i] / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
    w               = w * ratio;
    double v        = ratio - w;
    if (sync[i] > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant == 0.f ? 0.0001f : formant;
    double ratio    = freq[i] / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
    w               = w * ratio;
    double v        = ratio - w;
    if (sync > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant == 0.f ? 0.0001f : formant;
    double ratio    = freq[i] / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
    w               = w * ratio;
    double v        = ratio - w;
    if (sync[i] > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant == 0.f ? 0.0001f : formant;
    double ratio    = freq[i] / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
    w               = w * ratio;
    double v        = ratio - w;
    if (sync > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant[i] == 0.f ? 0.0001f : formant[i];
    double ratio    = freq / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
    w               = w * ratio;
    double v        = ratio - w;
    if (sync[i] > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant[i] == 0.f ? 0.0001f : formant[i];
    double ratio    = freq / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
    w               = w * ratio;
    double v        = ratio - w;
    if (sync > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant[i] == 0.f ? 0.0001f : formant[i];
    double ratio    = freq / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
    w               = w * ratio;
    double v        = ratio - w;
    if (sync[i] > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant[i] == 0.f ? 0.0001f : formant[i];
    double ratio    = freq / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
    w               = w * ratio;
    double v        = ratio - w;
    if (sync > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant == 0.f ? 0.0001f : formant;
    double ratio    = freq / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
    w               = w * ratio;
    double v        = ratio - w;
    if (sync[i] > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
    double form     = formant == 0.f ? 0.0001f : formant;
    double ratio    = freq / form;
    ratio           = ratio < 0.f ? -ratio : ratio;
    double w        = width[i] <= 0.01f ? 0.01f : width[i] >= 0.99f ? 0.99f : width[i];
    w               = w * ratio;
    double v        = ratio - w;
    if (sync > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
  double form     = formant == 0.f ? 0.0001f : formant;
  double ratio    = freq / form;
  ratio           = ratio < 0.f ? -ratio : ratio;
  double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
  w               = w * ratio;
  double v        = ratio - w;
  double out1     = 0.f;
  for (int i=0; i < nSamples; ++i) {
    if (sync[i] > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
    if (counter + pos> turnaround && pos <= step) {
      if (counter + pos > ratio) {
        counter = 0;
        offset = -1.f;
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
  double form     = formant == 0.f ? 0.0001f : formant;
  double ratio    = freq / form;
  ratio           = ratio < 0.f ? -ratio : ratio;
  double w        = width <= 0.01f ? 0.01f : width >= 0.99f ? 0.99f : width;
  w               = w * ratio;
  double v        = ratio - w;
  double out1     = 0.f;
  for (int i=0; i < nSamples; ++i) {
    if (sync > 0.f && lastsync <= 0.f) {
      pos = 0.f;
      counter = 0;
      offset = -1.f;
    } else {
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
      } else {
        offset = last;
        counter = 0;
        turnaround = (1.f - offset) * 0.5 * w;
      }

    }
    if (offset <= -1.f) {
      out1 = algorithm(pos + counter, step, w, v);
    } else {
      out1 = algorithm_retrig(pos + counter, step, w, v, offset, turnaround);
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
