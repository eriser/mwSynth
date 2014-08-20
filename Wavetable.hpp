#pragma once

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
     * @param input Array of 2 input samples.
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
     * @param mipmap        mipmap index
     * @param phase         sampling point [0.0 .. 1.0)
     * @param pInterpolator interpolator configuration pointer
     */
    float Sample_FPU(int mipmap, float phase, const mwInterpolator* pInterpolator) const;

    /**
     * SSE version of @p Sample method.
     */
    __m128 Sample_SSE(int mipmap, __m128 phase, const mwInterpolator* pInterpolator) const;

    /**
     * AVX version of @p Sample method.
     */
    float Sample_AVX(int mipmap, float phase, const mwInterpolator* pInterpolator) const;


    /**
     * Synthesize samples buffer.
     * @param samplesNum   Number of samples to generate.
     * @param pFreq        Signal frequency for each sample.
     * @param pCtx         Synthesis context.
     * @param[out] pOutput Buffer to write.
     */
    void Synth(size_t samplesNum, const float* pFreq, mwWaveSynthContext* pCtx, float* pOutput);

    /**
     * FPU version of @p Synth method.
     */
    void Synth_FPU(size_t samplesNum, const float* pFreq, mwWaveSynthContext* pCtx, float* pOutput);

    /**
     * SSE version of @p Synth method.
     */
    void Synth_SSE(size_t samplesNum, const float* pFreq, mwWaveSynthContext* pCtx, float* pOutput);
};

void* malloc_simd(const size_t size);
void free_simd(void* ptr);
