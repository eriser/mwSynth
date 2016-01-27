#include "stdafx.hpp"
#include "Wavetable.hpp"
#include "Math.hpp"

#include <immintrin.h> // for AVX

namespace mvSynth {

// sum YMM register horizontally
inline float _mm256_hsum(__m256 x)
{
    const __m128 hiQuad = _mm256_extractf128_ps(x, 1);
    const __m128 loQuad = _mm256_castps256_ps128(x);
    const __m128 sumQuad = _mm_add_ps(loQuad, hiQuad);
    const __m128 loDual = sumQuad;
    const __m128 hiDual = _mm_movehl_ps(sumQuad, sumQuad);
    const __m128 sumDual = _mm_add_ps(loDual, hiDual);
    const __m128 lo = sumDual;
    const __m128 hi = _mm_shuffle_ps(sumDual, sumDual, 0x1);
    const __m128 sum = _mm_add_ss(lo, hi);
    return _mm_cvtss_f32(sum);
}

float WaveTable::Sample_AVX(int mipmap, float phase, const Interpolator& interpolator) const
{
    const float* waveData = mData[mipmap];
    const float const* const* filterData = interpolator.data;

    int mipmapSize = mRootSize >> mipmap;
    phase *= (float)mipmapSize;
    unsigned int id = (unsigned int)phase;

    __m256 phase_v = _mm256_set1_ps(phase);
    __m256 base_phase = _mm256_sub_ps(phase_v, _mm256_floor_ps(phase_v));
    // scale phase to 0..INTERPOLATOR_FILTER_PHASES-1
    __m256 phase_scaled = _mm256_mul_ps(base_phase,
                                        _mm256_set1_ps(interpolator.mFilterPhasesF));
    // phase interpolation coefficient
    __m256 phase_factor = _mm256_sub_ps(phase_scaled, _mm256_floor_ps(phase_scaled));

    // select base phase indicies
    // _mm256_extractf128_si256(_mm256_cvttps_epi32(phase_scaled));
    int phase_sel = 0;

    // FIR filter
    __m256 sum_v = _mm256_setzero_ps();
    for (int i = 0; i < interpolator.mFilterSize; i += 8)
    {
        __m256 f0 = _mm256_load_ps(filterData[phase_sel] + i);
        __m256 f1 = _mm256_load_ps(filterData[phase_sel + 1] + i);

        __m256 f = _mm256_sub_ps(f0, _mm256_mul_ps(phase_factor, _mm256_sub_ps(f0, f1)));
        //__m256 f = _mm256_fmadd_ps(phase_factor, _mm256_sub_ps(f1, f0), f0);

        __m256 x = _mm256_loadu_ps(waveData + (id + i));

        sum_v = _mm256_add_ps(sum_v, _mm256_mul_ps(f, x));
        // sum_v = _mm256_fmadd_ps(f, x, sum_v);
    }

    return _mm256_hsum(sum_v);
}

void WaveTable::Synth_AVX(size_t samplesNum, const float* freqBuff, WaveTableContext& ctx,
                          const Interpolator& interpolator, float* output) const
{
    __m128 samples;
    float phases[4];

    __m128 ones = _mm_set1_ps(1.0f);

    for (size_t i = 0; i < samplesNum; i += 4)
    {
        samples = _mm_setzero_ps();

        for (size_t j = 0; j < ctx.mPhases.size(); ++j)
        {
            // TODO: optimize
            phases[0] = ctx.mPhases[j] + freqBuff[i] * 0.5f;
            phases[1] = ctx.mPhases[j] + freqBuff[i];
            phases[2] = phases[1] + freqBuff[i + 1] * 0.5f;
            phases[3] = phases[1] + freqBuff[i + 1];
            phases[4] = phases[3] + freqBuff[i + 2] * 0.5f;
            phases[5] = phases[3] + freqBuff[i + 2];
            phases[6] = phases[5] + freqBuff[i + 3] * 0.5f;
            phases[7] = phases[5] + freqBuff[i + 3];

            if (phases[0] > 1.0f) phases[0] -= 1.0f;
            if (phases[1] > 1.0f) phases[1] -= 1.0f;
            if (phases[2] > 1.0f) phases[2] -= 1.0f;
            if (phases[3] > 1.0f) phases[3] -= 1.0f;
            if (phases[4] > 1.0f) phases[4] -= 1.0f;
            if (phases[5] > 1.0f) phases[5] -= 1.0f;
            if (phases[6] > 1.0f) phases[6] -= 1.0f;
            if (phases[7] > 1.0f) phases[7] -= 1.0f;
            ctx.mPhases[j] = phases[7];

            float ratio = freqBuff[i] * (float)(mRootSize);
            int mipmap = Math::log2_int(ratio);
            if (mipmap > mMipsNum)
                mipmap = mMipsNum;
            if (mipmap < 0)
                mipmap = 0;

            __m128 phases_v = _mm_set_ps(phases[3], phases[2], phases[1], phases[0]);
            samples = _mm_add_ps(samples, Sample_SSE(mipmap, phases_v, interpolator));
        }

        output[i] = ctx.Downsample(samples.m128_f32);
        output[i + 1] = ctx.Downsample(samples.m128_f32 + 2);
    }
}

} // namespace mvSynth