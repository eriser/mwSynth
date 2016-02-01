#pragma once

#include "Interpolator.hpp"
#include <vector>

namespace mvSynth {

#define ALIGNED __declspec(align(16))

#define IIR_FILTER_SIZE 12
#define MAX_VOICES 16

#define MIPMAP_BLEND_TRESHOLD 0.98f

class ALIGNED WaveTableContext final
{
    friend class WaveTable;

    ALIGNED static const double a[];
    ALIGNED static const double b[];
    ALIGNED double mX[IIR_FILTER_SIZE];
    ALIGNED double mY[IIR_FILTER_SIZE];
    std::vector<float> mPhases;

public:
    WaveTableContext();

    /**
     * Clear filter history.
     */
    void Reset();

    /**
     * Set subvoices number and set initial phases.
     * @param numVoices Number of subvoices
     */
    void Init(size_t numVoices, float* newPhases);

    /**
     * Downsamples two samples into one using IIR lowpass filter.
     * @param input Array of 2 input samples.
     */
    float Downsample(float* input);

    float Downsample_SSE(float* input);
};

class WaveTable
{
    int mRootSize; // size of 0th mipmap in samples
    float mRootSizeF;

    int mMipsNum;  // total number of mipmaps = log2(m_RootSize)+1
    float** mData;

public:
    WaveTable();
    ~WaveTable();

    void Release();

    /*
        order:
        0 -> 1 sample, 1 -> 2 samples, 2 -> 4 samples, ...
    */
    int LoadData(float* pData, int order, const Interpolator& interpolator);

    /**
     * Interpolate sample from specified mipmap level.
     * Simple FPU version.
     * @param mipmap        mipmap index
     * @param phase         sampling point [0.0 .. 1.0)
     * @param pInterpolator interpolator configuration pointer
     */
    float Sample_FPU(int mipmap, float phase, const Interpolator& interpolator) const;

    /**
     * SSE version of @p Sample method.
     */
    __m128 Sample_SSE(int mipmap, __m128 phase, const Interpolator& interpolator) const;

    /**
     * AVX version of @p Sample method.
     */
    float Sample_AVX(int mipmap, float phase, const Interpolator& interpolator) const;


    /**
     * Synthesize samples buffer.
     * @param samplesNum   Number of samples to generate.
     * @param freqBuff     Signal frequency for each sample.
     * @param ctx          Synthesiser context.
     * @param interpolator Wavetable interpolator.
     * @param[out] output  Buffer to write.
     */
    void Synth_FPU(size_t samplesNum, const float* freqBuff, WaveTableContext& ctx,
                   const Interpolator& interpolator, float* output) const;

    /**
     * SSE version of @p Synth method.
     */
    void Synth_SSE(size_t samplesNum, const float* freqBuff, WaveTableContext& ctx,
                   const Interpolator& interpolator, float* output) const;

    /**
     * AVX version of @p Synth method.
     */
    void Synth_AVX(size_t samplesNum, const float* freqBuff, WaveTableContext& ctx,
                   const Interpolator& interpolator, float* output) const;
};

void* malloc_simd(const size_t size);
void free_simd(void* ptr);

} // namespace mvSynth