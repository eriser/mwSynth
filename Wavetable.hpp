#pragma once
#include "Math.hpp"

class mwInterpolator final
{
    friend class mwWaveTable;

    float** data;
    int m_FilterSize;
    int m_FilterPhases;
    float m_FilterPhasesF;

public:

    mwInterpolator();
    ~mwInterpolator();
    int Setup(int filterSize, int filterPhases);
};

#define IIR_FILTER_SIZE 11

class mwWaveSynthContext final
{
    friend class mwWaveTable;
    float m_Phase;
    float x[IIR_FILTER_SIZE];
    float y[IIR_FILTER_SIZE];

public:
    const mwInterpolator* pInterpolator;

    mwWaveSynthContext();

    /**
     * Downsamples two samples into one using IIR lowpass filter.
     * @param x0,x1 Input samples.
     */
    float Downsample(float* input);
};

class mwWaveTable
{
    int m_RootSize; // size of 0th mipmap in samples
    int m_MipsNum;  // total number of mipmaps = log2(m_RootSize)+1
    float** m_ppData;

public:
    mwWaveTable();
    ~mwWaveTable();

    void Release();

    /*
        order:
        0 -> 1 sample, 1 -> 2 samples, 2 -> 4 samples, ...
    */
    int LoadData(float* pData, int order, const mwInterpolator* pInterpolator);

    /**
     * Interpolate sample from specified mipmap level.
     * Simple FPU version.
     * @param ratio sampling ratio: 1.0 - normal, 2.0 - 2x less harmonics, etc...
     * @param phase sampling point [0.0 .. 1.0)
     */
    float Sample(float ratio, float phase, const mwInterpolator* pInterpolator) const;
    float Sample_SSE(float ratio, float phase, const mwInterpolator* pInterpolator) const;
    float Sample_AVX(float ratio, float phase, const mwInterpolator* pInterpolator) const;

    /**
     * Synthesize samples buffer.
     * @param samplesNum   Number of samples to generate.
     * @param pFreq        Signal frequency for each sample.
     * @param pCtx         Synthesis context.
     * @param[out] pOutput Buffer to write.
     */
    void Synth(size_t samplesNum, const float* pFreq, mwWaveSynthContext* pCtx, float* pOutput);

    /**
     * SSE version of @p Synth method.
     */
    void Synth_SSE(size_t samplesNum, const float* pFreq, mwWaveSynthContext* pCtx, float* pOutput);
};