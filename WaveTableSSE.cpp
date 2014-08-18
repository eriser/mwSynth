#include "stdafx.hpp"
#include "Wavetable.hpp"

#include <xmmintrin.h>

// sum all XMM register horizontally
__forceinline float _mm_hsum(__m128 x)
{
    const __m128 t = _mm_add_ps(x, _mm_movehl_ps(x, x));
    const __m128 sum = _mm_add_ss(t, _mm_shuffle_ps(t, t, 1));
    return _mm_cvtss_f32(sum);
}

float mwWaveTable::Sample_SSE(float ratio, float phase, const mwInterpolator* pInterpolator) const
{
    int mipmap = Math::log2_int(ratio);
    if (mipmap >= m_MipsNum) return 0.0f;
    if (mipmap < 0) mipmap = 0;

    const float* pSrc = m_ppData[mipmap];
    float** pFilter = pInterpolator->data;

    unsigned int mipmapSize = m_RootSize >> mipmap;
    phase *= (float)mipmapSize;
    unsigned int id = (unsigned int)phase;

    __m128 phase_v = _mm_set1_ps(phase);
    __m128 base_phase = _mm_sub_ps(phase_v, _mm_floor_ps(phase_v));

    // scale phase to 0..INTERPOLATOR_FILTER_PHASES-1
    __m128 phase_scaled = _mm_mul_ps(base_phase, _mm_set1_ps(pInterpolator->m_FilterPhasesF));

    // select base phase indicies
    int phase_sel = _mm_cvttss_si32(phase_scaled); // float -> int casting

    // phase interpolation coefficient
    __m128 phase_factor = _mm_sub_ps(phase_scaled, _mm_floor_ps(phase_scaled));


    // FIR filter
    __m128 sum_v = _mm_setzero_ps();
    for (int i = 0; i < pInterpolator->m_FilterSize; i += 4)
    {
        __m128 f0 = _mm_load_ps(pFilter[phase_sel] + i);
        __m128 f1 = _mm_load_ps(pFilter[phase_sel + 1] + i);
        __m128 x = _mm_loadu_ps(pSrc + (id + i));

        /// sum_v += x * (f0 - phase_factor * (f0-f1))
        sum_v = _mm_add_ps(sum_v, _mm_mul_ps(_mm_sub_ps(f0, _mm_mul_ps(phase_factor, _mm_sub_ps(f0, f1))),
                                             x));
    }

    return _mm_hsum(sum_v);
}

void mwWaveTable::Synth_SSE(size_t samplesNum, const float* pFreq, mwWaveSynthContext* pCtx,
                            float* pOutput)
{
    float phases[4];
    float samples[4];

    phases[3] = pCtx->m_Phase;
    for (size_t i = 0; i < samplesNum; i += 2)
    {
        // TODO: optimize
        phases[0] = phases[3] + pFreq[i] * 0.5f;
        phases[1] = phases[0] + pFreq[i] * 0.5f;
        phases[2] = phases[1] + pFreq[i + 1] * 0.5f;
        phases[3] = phases[2] + pFreq[i + 1] * 0.5f;
        if (phases[0] > 1.0f) phases[0] -= 1.0f;
        if (phases[1] > 1.0f) phases[1] -= 1.0f;
        if (phases[2] > 1.0f) phases[2] -= 1.0f;
        if (phases[3] > 1.0f) phases[3] -= 1.0f;

        samples[0] = Sample_SSE(pFreq[i] * (float)(m_RootSize), phases[0], pCtx->pInterpolator);
        samples[1] = Sample_SSE(pFreq[i] * (float)(m_RootSize), phases[1], pCtx->pInterpolator);
        samples[2] = Sample_SSE(pFreq[i + 1] * (float)(m_RootSize), phases[2], pCtx->pInterpolator);
        samples[3] = Sample_SSE(pFreq[i + 1] * (float)(m_RootSize), phases[3], pCtx->pInterpolator);

        pOutput[i] =     pCtx->Downsample(samples);
        pOutput[i + 1] = pCtx->Downsample(samples + 2);
    }
    pCtx->m_Phase = phases[3];
}