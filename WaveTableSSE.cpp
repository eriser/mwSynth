#include "stdafx.hpp"
#include "Wavetable.hpp"
#include "Math.hpp"
#include <xmmintrin.h>
#include <smmintrin.h> // for SSE4

__m128 mwWaveTable::Sample_SSE(int mipmap, __m128 phase, const mwInterpolator* pInterpolator) const
{
    const float* pSrc;
    float** pFilter = pInterpolator->data;
    int mipmapSize = m_RootSize >> mipmap;

    _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);

    __m128 phase_v = _mm_mul_ps(phase, _mm_set1_ps((float)mipmapSize));
    __m128i id = _mm_cvtps_epi32(phase_v);
    __m128 base_phase = _mm_sub_ps(phase_v, _mm_floor_ps(phase_v));

    // scale phase to 0..INTERPOLATOR_FILTER_PHASES-1
    __m128 phase_scaled = _mm_mul_ps(base_phase, _mm_set1_ps(pInterpolator->m_FilterPhasesF));

    // select base phase indicies
    __m128i phase_sel = _mm_cvtps_epi32(phase_scaled);

    // phase interpolation coefficient
    __m128 phase_factor = _mm_sub_ps(phase_scaled, _mm_floor_ps(phase_scaled));

    __m128 sums[4];
    for (int k = 0; k < 4; ++k)
    {
        pSrc = m_ppData[mipmap] + id.m128i_i32[k];
        __m128 sum = _mm_setzero_ps();
        __m128 factor = _mm_set_ps1(phase_factor.m128_f32[k]);
        int j = phase_sel.m128i_i32[k];
        for (int i = 0; i < pInterpolator->m_FilterSize; i += 4)
        {
            __m128 f0 = _mm_load_ps(pFilter[j] + i);
            __m128 f1 = _mm_load_ps(pFilter[j + 1] + i);
            __m128 x = _mm_loadu_ps(pSrc + i);
            __m128 term = _mm_sub_ps(f0, _mm_mul_ps(factor, _mm_sub_ps(f0, f1)));
            sum = _mm_add_ps(sum, _mm_mul_ps(term, x));
        }
        sums[k] = sum;
    }

    _MM_TRANSPOSE4_PS(sums[0], sums[1], sums[2], sums[3]);
    return _mm_add_ps(_mm_add_ps(sums[0], sums[1]), _mm_add_ps(sums[2], sums[3]));
}

void mwWaveTable::Synth_SSE(size_t samplesNum, const float* pFreq, mwWaveSynthContext* pCtx,
                            float* pOutput)
{
    __m128 samples;
    float phases[4];

    __m128 ones = _mm_set1_ps(1.0f);
    
    for (size_t i = 0; i < samplesNum; i += 2)
    {
        samples = _mm_setzero_ps();

        for (size_t j = 0; j < pCtx->phases.size(); ++j)
        {
            // TODO: optimize
            phases[0] = pCtx->phases[j] + pFreq[i] * 0.5f;
            phases[1] = pCtx->phases[j] + pFreq[i];
            phases[2] = phases[1] + pFreq[i + 1] * 0.5f;
            phases[3] = phases[1] + pFreq[i + 1];
            if (phases[0] > 1.0f) phases[0] -= 1.0f;
            if (phases[1] > 1.0f) phases[1] -= 1.0f;
            if (phases[2] > 1.0f) phases[2] -= 1.0f;
            if (phases[3] > 1.0f) phases[3] -= 1.0f;
            pCtx->phases[j] = phases[3];

            float ratio = pFreq[i] * (float)(m_RootSize);
            int mipmap = Math::log2_int(ratio);
            if (mipmap > m_MipsNum) mipmap = m_MipsNum;
            if (mipmap < 0) mipmap = 0;

            __m128 phases_v = _mm_set_ps(phases[3], phases[2], phases[1], phases[0]);
            samples = _mm_add_ps(samples, Sample_SSE(mipmap, phases_v, pCtx->pInterpolator));
        }

        pOutput[i] =     pCtx->Downsample(samples.m128_f32);
        pOutput[i + 1] = pCtx->Downsample(samples.m128_f32 + 2);
    }
}
