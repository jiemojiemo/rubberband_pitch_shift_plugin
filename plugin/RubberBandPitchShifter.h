//
// Created by user on 4/7/22.
//

#ifndef PITCHSHIFTERPLUGIN_RUBBERBANDPITCHSHIFTER_H
#define PITCHSHIFTERPLUGIN_RUBBERBANDPITCHSHIFTER_H

#pragma once
#include "RubberBandStretcher.h"
#include "base/RingBuffer.h"
#include <iostream>
using RubberBand::RubberBandStretcher;
using RubberBand::RingBuffer;

// copy from rubberband/ladspa-lv2 and delete some unused functions
class RubberBandPitchShifter {
public:
  RubberBandPitchShifter(int sampleRate, size_t channels)
      : m_latency(nullptr), m_cents(nullptr), m_semitones(nullptr),
        m_octaves(nullptr), m_crispness(nullptr), m_formant(nullptr),
        m_wetDry(nullptr), m_ratio(1.0), m_prevRatio(1.0),
        m_currentCrispness(-1), m_currentFormant(false), m_blockSize(1024),
        m_reserve(8192), m_bufsize(0), m_minfill(0),
        m_stretcher(new RubberBandStretcher(
            sampleRate, channels,
            RubberBandStretcher::OptionProcessRealTime |
                RubberBandStretcher::OptionPitchHighConsistency)),
        m_sampleRate(sampleRate), m_channels(channels) {
    m_input = new float *[m_channels];
    m_output = new float *[m_channels];

    m_outputBuffer = new RingBuffer<float> *[m_channels];
    m_delayMixBuffer = new RingBuffer<float> *[m_channels];
    m_scratch = new float *[m_channels];
    m_inptrs = new float *[m_channels];

    m_bufsize = m_blockSize + m_reserve + 8192;

    for (size_t c = 0; c < m_channels; ++c) {

      m_input[c] = 0;
      m_output[c] = 0;

      m_outputBuffer[c] = new RingBuffer<float>(m_bufsize);
      m_delayMixBuffer[c] = new RingBuffer<float>(m_bufsize);

      m_scratch[c] = new float[m_bufsize];
      for (size_t i = 0; i < m_bufsize; ++i) {
        m_scratch[c][i] = 0.f;
      }

      m_inptrs[c] = 0;
    }

    activateImpl();
  }

  ~RubberBandPitchShifter()
  {
    delete m_stretcher;
    for (size_t c = 0; c < m_channels; ++c) {
      delete m_outputBuffer[c];
      delete m_delayMixBuffer[c];
      delete[] m_scratch[c];
    }
    delete[] m_outputBuffer;
    delete[] m_delayMixBuffer;
    delete[] m_inptrs;
    delete[] m_scratch;
    delete[] m_output;
    delete[] m_input;
  }

  void activateImpl()
  {
    updateRatio();
    m_prevRatio = m_ratio;
    m_stretcher->reset();
    m_stretcher->setPitchScale(m_ratio);

    for (size_t c = 0; c < m_channels; ++c) {
      m_outputBuffer[c]->reset();
    }

    for (size_t c = 0; c < m_channels; ++c) {
      m_delayMixBuffer[c]->reset();
      m_delayMixBuffer[c]->zero(getLatency());
    }

    for (size_t c = 0; c < m_channels; ++c) {
      for (size_t i = 0; i < m_bufsize; ++i) {
        m_scratch[c][i] = 0.f;
      }
    }

    m_minfill = 0;

    m_stretcher->process(m_scratch, m_reserve, false);
  }

  void runImpl(uint32_t insamples)
  {
    for (size_t c = 0; c < m_channels; ++c) {
      m_delayMixBuffer[c]->write(m_input[c], insamples);
    }

    size_t offset = 0;

    // We have to break up the input into chunks like this because
    // insamples could be arbitrarily large and our output buffer is
    // of limited size

    while (offset < insamples) {

      size_t block = m_blockSize;
      if (offset + block > insamples) {
        block = insamples - offset;
      }

      runImpl(block, offset);

      offset += block;
    }

    float mix = 0.0;
    if (m_wetDry) mix = *m_wetDry;

    for (size_t c = 0; c < m_channels; ++c) {
      if (mix > 0.0) {
        for (size_t i = 0; i < insamples; ++i) {
          float dry = m_delayMixBuffer[c]->readOne();
          m_output[c][i] *= (1.0 - mix);
          m_output[c][i] += dry * mix;
        }
      } else {
        m_delayMixBuffer[c]->skip(insamples);
      }
    }
  }

  void runImpl(uint32_t insamples, uint32_t offset)
  {
    updateRatio();
    if (m_ratio != m_prevRatio) {
      m_stretcher->setPitchScale(m_ratio);
      m_prevRatio = m_ratio;
    }

    if (m_latency) {
      *m_latency = getLatency();
    }

        updateCrispness();
        updateFormant();

    const int samples = insamples;
    int processed = 0;
    size_t outTotal = 0;

    while (processed < samples) {

      // never feed more than the minimum necessary number of
      // samples at a time; ensures nothing will overflow internally
      // and we don't need to call setMaxProcessSize

      int toCauseProcessing = m_stretcher->getSamplesRequired();
      int inchunk = std::min(samples - processed, toCauseProcessing);

      for (size_t c = 0; c < m_channels; ++c) {
        m_inptrs[c] = &(m_input[c][offset + processed]);
      }

      m_stretcher->process(m_inptrs, inchunk, false);

      processed += inchunk;

      int avail = m_stretcher->available();
      int writable = m_outputBuffer[0]->getWriteSpace();

      int outchunk = avail;
      if (outchunk > writable) {
        std::cerr << "RubberBandPitchShifter::runImpl: buffer is not large enough: size = " << m_outputBuffer[0]->getSize() << ", chunk = " << outchunk << ", space = " << writable << " (buffer contains " << m_outputBuffer[0]->getReadSpace() << " unread)" << std::endl;
        outchunk = writable;
      }

      size_t actual = m_stretcher->retrieve(m_scratch, outchunk);
      outTotal += actual;

      for (size_t c = 0; c < m_channels; ++c) {
        m_outputBuffer[c]->write(m_scratch[c], actual);
      }
    }

    for (size_t c = 0; c < m_channels; ++c) {
      int toRead = m_outputBuffer[c]->getReadSpace();
      if (toRead < samples && c == 0) {
        std::cerr << "RubberBandPitchShifter::runImpl: buffer underrun: required = " << samples << ", available = " << toRead << std::endl;
      }
      int chunk = std::min(toRead, samples);
      m_outputBuffer[c]->read(&(m_output[c][offset]), chunk);
    }

    size_t fill = m_outputBuffer[0]->getReadSpace();
    if (fill < m_minfill || m_minfill == 0) {
      m_minfill = fill;
      //        cerr << "minfill = " << m_minfill << endl;
    }
  }


  int getLatency() const
  {
    return m_reserve;
  }

  void updateRatio()
  {
    // The octaves, semitones, and cents parameters are supposed to be
    // integral: we want to enforce that, just to avoid
    // inconsistencies between hosts if some respect the hints more
    // than others

#ifdef RB_PLUGIN_LADSPA

    // But we don't want to change the long-standing behaviour of the
    // LADSPA plugin, so let's leave this as-is and only do "the right
    // thing" for LV2
    double oct = (m_octaves ? *m_octaves : 0.0);
    oct += (m_semitones ? *m_semitones : 0.0) / 12;
    oct += (m_cents ? *m_cents : 0.0) / 1200;
    m_ratio = pow(2.0, oct);

#else

    // LV2

    double octaves = round(m_octaves ? *m_octaves : 0.0);
    if (octaves < -2.0) octaves = -2.0;
    if (octaves >  2.0) octaves =  2.0;

    double semitones = round(m_semitones ? *m_semitones : 0.0);
    if (semitones < -12.0) semitones = -12.0;
    if (semitones >  12.0) semitones =  12.0;

    double cents = round(m_cents ? *m_cents : 0.0);
    if (cents < -100.0) cents = -100.0;
    if (cents >  100.0) cents =  100.0;

    m_ratio = pow(2.0,
                  octaves +
                      semitones / 12.0 +
                      cents / 1200.0);
#endif
  }

  void updateCrispness()
  {
    if (!m_crispness) return;

    int c = lrintf(*m_crispness);
    if (c == m_currentCrispness) return;
    if (c < 0 || c > 3) return;
    RubberBandStretcher *s = m_stretcher;

    switch (c) {
    case 0:
      s->setPhaseOption(RubberBandStretcher::OptionPhaseIndependent);
      s->setTransientsOption(RubberBandStretcher::OptionTransientsSmooth);
      break;
    case 1:
      s->setPhaseOption(RubberBandStretcher::OptionPhaseLaminar);
      s->setTransientsOption(RubberBandStretcher::OptionTransientsSmooth);
      break;
    case 2:
      s->setPhaseOption(RubberBandStretcher::OptionPhaseLaminar);
      s->setTransientsOption(RubberBandStretcher::OptionTransientsMixed);
      break;
    case 3:
      s->setPhaseOption(RubberBandStretcher::OptionPhaseLaminar);
      s->setTransientsOption(RubberBandStretcher::OptionTransientsCrisp);
      break;
    }

    m_currentCrispness = c;
  }

  void updateFormant()
  {
    if (!m_formant) return;

    bool f = (*m_formant > 0.5f);
    if (f == m_currentFormant) return;

    RubberBandStretcher *s = m_stretcher;

    s->setFormantOption(f ?
                          RubberBandStretcher::OptionFormantPreserved :
                          RubberBandStretcher::OptionFormantShifted);

    m_currentFormant = f;
  }

  size_t m_blockSize;
  size_t m_reserve;
  size_t m_bufsize;
  size_t m_minfill;

  RubberBand::RubberBandStretcher *m_stretcher;
  RubberBand::RingBuffer<float> **m_outputBuffer;
  RubberBand::RingBuffer<float> **m_delayMixBuffer;
  float **m_scratch;
  float **m_inptrs;

  int m_sampleRate;
  size_t m_channels;

  float **m_input;
  float **m_output;
  float *m_latency;
  float *m_cents;
  float *m_semitones;
  float *m_octaves;
  float *m_crispness;
  float *m_formant;
  float *m_wetDry;
  double m_ratio;
  double m_prevRatio;
  int m_currentCrispness;
  bool m_currentFormant;
};


#endif // PITCHSHIFTERPLUGIN_RUBBERBANDPITCHSHIFTER_H
