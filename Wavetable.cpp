#include "stdafx.hpp"
#include "Wavetable.hpp"
#include "Math.hpp"

void free_simd(void* ptr)
{
#if defined WIN32
    _aligned_free(ptr);
#else
    free(ptr);
#endif
}

void* malloc_simd(const size_t size)
{
#if defined WIN32
    return _aligned_malloc(size, 16);
#elif defined __linux__
    return memalign(16, size);
#elif defined __MACH__
    return malloc(size);
#else
    return valloc(size);
#endif
}

mwWaveSynthContext::mwWaveSynthContext()
{
    m_Phase = 0.0f;
    for (int i = 0; i < IIR_FILTER_SIZE; ++i)
    {
        x[i] = 0.0f;
        y[i] = 0.0f;
    }
}

float mwWaveSynthContext::Downsample(float* input)
{
    // TODO: cleanup, make SSE and AVX versions as well
    for (int i = 0; i < 2; ++i)
    {
        for (int j = IIR_FILTER_SIZE - 1; j > 0; --j)
        {
            x[j] = x[j - 1];
            y[j] = y[j - 1];
        }

        // lowpass Chebyshev filter generated with:
        // http://www-users.cs.york.ac.uk/~fisher/mkfilter/
        //
        // filtertype   =   Chebyshev
        // passtype     =   Lowpass
        // ripple   =   -0.5
        // order    =   10
        // samplerate   =   88200
        // corner1  =   20000
        x[0] = input[i] / 7.290426784e+03f;
        y[0] = x[0] + x[10] +
               10.0f * (x[9] + x[1]) +
               45.0f * (x[8] + x[2]) +
               120.0f * (x[7] + x[3]) +
               210.0f * (x[6] + x[4]) +
               252.0f * x[5] +
               ( -0.1979824657f * y[10]) +
               (  1.2383760478f * y[9]) +
               ( -4.1332027761f * y[8]) +
               (  9.3602145643f * y[7]) +
               (-15.7280128999f * y[6]) +
               ( 20.3867628124f * y[5]) +
               (-20.6672437866f * y[4]) +
               ( 16.2722769093f * y[3]) +
               ( -9.6890187121f * y[2]) +
               (  4.0173721372f * y[1]);
    }

    return y[0];
}

mwInterpolator::mwInterpolator()
{
    m_FilterPhases = m_FilterSize = 0;
    m_FilterPhasesF = 0.0f;
    data = nullptr;
}

mwInterpolator::~mwInterpolator()
{
    if (data)
    {
        for (int i = 0; i <= m_FilterPhases; i++)
            free_simd(data[i]);

        free(data);
        data = nullptr;
    }
}


int mwInterpolator::Setup(int filterSize, int filterPhases)
{
    data = (float**)malloc(sizeof(void*) * (filterPhases + 1));

    for (int i = 0; i <= filterPhases; i++)
        data[i] = (float*)malloc_simd(sizeof(float) * filterSize);

    m_FilterPhases = filterPhases;
    m_FilterSize = filterSize;
    m_FilterPhasesF = (float)filterPhases;

    /*
        Interpolator cutoff frequency offset. Increased to flatten passband
    */
    const double cutoff = 1.1f;

    for (int j = 0; j < filterPhases; j++)
    {
        double phase = (float)j / (float)filterPhases;

        for (int i = 0; i < filterSize; i++)
        {
            double x = (float)(i - filterSize / 2);
            double y = Math::Sinc((x + phase) * cutoff);
            double window = Math::BlackmanHarris(((float)i + phase) / (float)filterSize);
            data[filterPhases - 1 - j][i] = (float)(y * window);
        }
    }


    double phase = -1.0 / (float)filterPhases;
    for (int i = 0; i < filterSize; i++)
    {
        double x = (float)(i - filterSize / 2);
        double y = Math::Sinc((x + phase) * cutoff);
        double window = Math::BlackmanHarris(((float)i + phase) / (float)filterSize);
        data[filterPhases][i] = (float)(y * window);
    }

    return 0;
}



mwWaveTable::mwWaveTable()
{
    m_ppData = nullptr;
    m_MipsNum = 0;
    m_RootSize = 0;
}

mwWaveTable::~mwWaveTable()
{
    Release();
}

void mwWaveTable::Release()
{
    if (m_ppData)
    {
        for (int i = 0; i < m_MipsNum; ++i)
        {
            free(m_ppData[i]);
            m_ppData[i] = nullptr;
        }

        free(m_ppData);
        m_ppData = nullptr;
    }

    m_MipsNum = 0;
    m_RootSize = 0;
}

/*
    order:
    0 -> 1 sample, 1 -> 2 samples, 2 -> 4 samples, ...
*/
int mwWaveTable::LoadData(float* pData, int order, const mwInterpolator* pInterpolator)
{
    int inSamples = 1 << order; //input samples
    m_MipsNum = order + 2;
    m_RootSize = inSamples * 2; // additional level for upsampled signal

    m_ppData = (float**)malloc(sizeof(float*) * m_MipsNum);

    // calculate FFT of original signal
    float* pBase = (float*)malloc(sizeof(float) * inSamples * 2);
    for (int i = 0; i < inSamples; i++)
    {
        pBase[2 * i] = pData[i] / inSamples; // real
        pBase[2 * i + 1] = 0.0f; // imaginary
    }
    //printf("--- MIPMAP %i ---\n", 1);
    //PrintFFT(pBase, inSamples);
    Math::FFT(pBase, inSamples);


    // upsample by IFFT
    float* pHelper = (float*)malloc(sizeof(float) * inSamples * 4);
    memcpy(pHelper, pBase, sizeof(float) * inSamples);
    memcpy(pHelper + 3 * inSamples, pBase + inSamples, sizeof(float) * inSamples);
    memset(pHelper + inSamples, 0, 2 * sizeof(float) * inSamples);
    Math::IFFT(pHelper, 2 * inSamples);
    //printf("\n--- MIPMAP %i ---\n", 0);
    //PrintFFT(pHelper, 2*inSamples);

    m_ppData[0] = (float*)malloc(sizeof(float) * (m_RootSize + pInterpolator->m_FilterSize));
    for (int i = 0; i < m_RootSize; i++)
        m_ppData[0][i] = pHelper[2 * i];
    memcpy(&m_ppData[0][m_RootSize], &m_ppData[0][0],
           sizeof(float) * (pInterpolator->m_FilterSize)); // fold


    m_ppData[1] = (float*)malloc(sizeof(float) * (m_RootSize / 2 + pInterpolator->m_FilterSize));
    memcpy(m_ppData[1], pData, sizeof(float) * m_RootSize / 2);
    memcpy(&m_ppData[1][m_RootSize / 2], &m_ppData[1][0],
           sizeof(float) * (pInterpolator->m_FilterSize)); // fold

    int samples = inSamples / 2;
    for (int i = 2; i < m_MipsNum; i++)
    {
        m_ppData[i] = (float*)malloc(sizeof(float) * (samples + pInterpolator->m_FilterSize));

        memcpy(pHelper, pBase, sizeof(float) * samples);
        memcpy(pHelper + samples, pBase + 2 * inSamples - samples, sizeof(float) * samples);
        //pHelper[samples] *= 1 << i;
        Math::IFFT(pHelper, samples);

        for (int j = 0; j < samples; j++)
            m_ppData[i][j] = pHelper[2 * j];

        for (int j = 0; j < pInterpolator->m_FilterSize; ++j) // fold
            m_ppData[i][j + samples] = m_ppData[i][j];

        samples /= 2;
    }

    free(pHelper);
    free(pBase);
    return 0;
}

/*
    Interpolate sample from specified mipmap level.
    Simple FPU version.
    * ratio - sampling ratio: 1.0f - normal, 2.0f - 2x less harmonics, etc...
    * phase - sampling point [0...1.0f)
*/
float mwWaveTable::Sample_FPU(int mipmap, float phase, const mwInterpolator* pInterpolator) const
{
    unsigned int mipmapSize = m_RootSize >> mipmap;
    phase *= (float)mipmapSize;
    unsigned int id = (unsigned int)phase;
    const float* pSrc = m_ppData[mipmap];

    float base_phase = phase - floor(phase);
    float phase_scaled = base_phase * (float)(
                             pInterpolator->m_FilterPhases); // scale phase to 0..INTERPOLATOR_FILTER_PHASES-1
    float phase_factor = phase_scaled - floor(phase_scaled); // phase interpolation coefficient

    // select base phase indicies
    int phase_sel = (int)phase_scaled;

    // FIR filter
    float sum = 0.0;
    for (int i = 0; i < pInterpolator->m_FilterSize; ++i)
    {
        float x0 = pSrc[id + i];
        float f0 = pInterpolator->data[phase_sel][i];
        float f1 = pInterpolator->data[phase_sel + 1][i];
        float f = f0 - phase_factor * (f0 - f1); // lineary interpolate filter coefficients
        sum += x0 * f;
    }

    return sum;
}

void mwWaveTable::Synth_FPU(size_t samplesNum, const float* pFreq, mwWaveSynthContext* pCtx,
                            float* pOutput)
{
    for (size_t i = 0; i < samplesNum; i++)
    {
        float freq = pFreq[0];
        float phaseA = pCtx->m_Phase + freq * 0.5f;
        if (phaseA > 1.0f) phaseA -= 1.0f;
        float phaseB = phaseA + freq * 0.5f;
        if (phaseB > 1.0f) phaseB -= 1.0f;

        float ratio = freq * (float)(m_RootSize);
        int mipmap = Math::log2_int(ratio);
        if (mipmap > m_MipsNum) mipmap = m_MipsNum;
        if (mipmap < 0) mipmap = 0;

        pCtx->m_Phase = phaseB;

        float samples[2];
        samples[0] = Sample_FPU(mipmap, phaseA, pCtx->pInterpolator);
        samples[1] = Sample_FPU(mipmap, phaseB, pCtx->pInterpolator);

        pOutput[i] = pCtx->Downsample(samples);
    }
}
