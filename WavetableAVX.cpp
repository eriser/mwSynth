#include "stdafx.hpp"
#include "Wavetable.hpp"

// sum YMM register horizontally
__forceinline float _mm256_hsum(__m256 x)
{
    // hiQuad = ( x7, x6, x5, x4 )
    const __m128 hiQuad = _mm256_extractf128_ps(x, 1);
    // loQuad = ( x3, x2, x1, x0 )
    const __m128 loQuad = _mm256_castps256_ps128(x);
    // sumQuad = ( x3 + x7, x2 + x6, x1 + x5, x0 + x4 )
    const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);
    // loDual = ( -, -, x1 + x5, x0 + x4 )
    const __m128 loDual = sumQuad;
    // hiDual = ( -, -, x3 + x7, x2 + x6 )
    const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);
    // sumDual = ( -, -, x1 + x3 + x5 + x7, x0 + x2 + x4 + x6 )
    const __m128 sumDual = _mm_add_ps(loDual, hiDual);
    // lo = ( -, -, -, x0 + x2 + x4 + x6 )
    const __m128 lo = sumDual;
    // hi = ( -, -, -, x1 + x3 + x5 + x7 )
    const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);
    // sum = ( -, -, -, x0 + x1 + x2 + x3 + x4 + x5 + x6 + x7 )
    const __m128 sum = _mm_add_ss(lo, hi);
    return _mm_cvtss_f32(sum);
}

float mwWaveTable::Sample_AVX(float ratio, float phase, const mwInterpolator* pInterpolator) const
{
    int mipmap = Math::log2_int(ratio);
    if (mipmap >= m_MipsNum) return 0.0f;
    if (mipmap < 0) mipmap = 0;

    const float* pSrc = m_ppData[mipmap];
    float** pFilter = pInterpolator->data;

    unsigned int mipmapSize = m_RootSize >> mipmap;
    phase *= (float)mipmapSize;
    unsigned int id = (unsigned int)phase;

    __m256 phase_v = _mm256_set1_ps(phase);
    __m256 base_phase = _mm256_sub_ps(phase_v, _mm256_floor_ps(phase_v));
    __m256 phase_scaled = _mm256_mul_ps(base_phase,
                                        _mm256_set1_ps(pInterpolator->m_FilterPhasesF)); // scale phase to 0..INTERPOLATOR_FILTER_PHASES-1
    __m256 phase_factor = _mm256_sub_ps(phase_scaled,
                                        _mm256_floor_ps(phase_scaled)); // phase interpolation coefficient

    // select base phase indicies
    int phase_sel = _mm256_cvttps_epi32(phase_scaled).m256i_i32[0]; // float -> int casting


    // FIR filter
    __m256 sum_v = _mm256_setzero_ps();
    for (int i = 0; i < pInterpolator->m_FilterSize; i += 8)
    {
        __m256 f0 = _mm256_load_ps(pFilter[phase_sel] + i);
        __m256 f1 = _mm256_load_ps(pFilter[phase_sel + 1] + i);

        //__m256 f = _mm256_sub_ps(f0, _mm256_mul_ps(phase_factor_v, _mm256_sub_ps(f0, f1)));
        __m256 f = _mm256_fmadd_ps(phase_factor, _mm256_sub_ps(f1, f0), f0);

        __m256 x = _mm256_loadu_ps(pSrc + (id + i));

        //sum_v = _mm256_add_ps(sum_v, _mm256_mul_ps(f, x));
        sum_v = _mm256_fmadd_ps(f, x, sum_v);
    }

    return _mm256_hsum(sum_v);
}
