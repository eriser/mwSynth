#include "stdafx.hpp"
#include "Wavetable.hpp"
#include "Math.hpp"

#include <xmmintrin.h>
#include <smmintrin.h> // for SSE4

namespace mvSynth {

float WaveTableContext::Downsample_SSE(float* input)
{
    for (int i = 0; i < 2; ++i)
    {
        for (int j = IIR_FILTER_SIZE - 1; j > 0; --j)
        {
            mX[j] = mX[j - 1];
            mY[j] = mY[j - 1];
        }

        mX[0] = input[i];

        __m128 sumA, sumB;
        sumA = _mm_mul_ps(_mm_load_ps(a), _mm_load_ps(mX));
        sumB = _mm_mul_ps(_mm_load_ps(b), _mm_load_ps(mY));
        sumA = _mm_add_ps(sumA, _mm_mul_ps(_mm_load_ps(a + 4), _mm_load_ps(mX + 4)));
        sumB = _mm_add_ps(sumB, _mm_mul_ps(_mm_load_ps(b + 4), _mm_load_ps(mY + 4)));
        sumA = _mm_add_ps(sumA, _mm_mul_ps(_mm_load_ps(a + 8), _mm_load_ps(mX + 8)));
        sumB = _mm_add_ps(sumB, _mm_mul_ps(_mm_load_ps(b + 8), _mm_load_ps(mY + 8)));
        __m128 sum = _mm_sub_ps(sumA, sumB);

        // horizontal sum
        mY[0] = (sum.m128_f32[0] + sum.m128_f32[1]) + (sum.m128_f32[2] + sum.m128_f32[3]);
    }

    return mY[0];
}

__m128 WaveTable::Sample_SSE(int mipmap, __m128 phase, const Interpolator& interpolator) const
{
    const float* src;
    const float const* const* filterData = interpolator.data;
    int mipmapSize = mRootSize >> mipmap;

    __m128 phase_v = _mm_mul_ps(phase, _mm_set1_ps((float)mipmapSize));
    __m128i id = _mm_cvtps_epi32(phase_v);
    __m128 base_phase = _mm_sub_ps(phase_v, _mm_floor_ps(phase_v));

    // scale phase to 0..INTERPOLATOR_FILTER_PHASES-1
    __m128 phase_scaled = _mm_mul_ps(base_phase, _mm_set1_ps(interpolator.mFilterPhasesF));

    // select base phase indicies
    __m128i phase_sel = _mm_cvtps_epi32(phase_scaled);

    // phase interpolation coefficient
    __m128 phase_factor = _mm_sub_ps(phase_scaled, _mm_floor_ps(phase_scaled));

    __m128 sums[4];
    for (int k = 0; k < 4; ++k)
    {
        src = mData[mipmap] + id.m128i_i32[k];
        __m128 sum = _mm_setzero_ps();
        __m128 factor = _mm_set_ps1(phase_factor.m128_f32[k]);
        int j = phase_sel.m128i_i32[k];
        for (int i = 0; i < interpolator.mFilterSize; i += 4)
        {
            __m128 f0 = _mm_load_ps(filterData[j] + i);
            __m128 f1 = _mm_load_ps(filterData[j + 1] + i);
            __m128 x = _mm_loadu_ps(src + i);
            // lineary interpolate filter coefficients
            __m128 term = _mm_sub_ps(f0, _mm_mul_ps(factor, _mm_sub_ps(f0, f1)));
            // FIR filtering
            sum = _mm_add_ps(sum, _mm_mul_ps(term, x));
        }
        sums[k] = sum;
    }

    _MM_TRANSPOSE4_PS(sums[0], sums[1], sums[2], sums[3]);
    return _mm_add_ps(_mm_add_ps(sums[0], sums[1]), _mm_add_ps(sums[2], sums[3]));
}

void WaveTable::Synth_SSE(size_t samplesNum, const float* freqBuff, WaveTableContext& ctx,
                          const Interpolator& interpolator, float* output) const
{
    __m128 samples;
    __m128 phases;

    const __m128 ones = _mm_set1_ps(1.0f);
    const __m128 freq0_factor = _mm_set_ps(1.0f, 1.0f, 1.0f, 0.5f);
    const __m128 freq1_factor = _mm_set_ps(1.0f, 0.5f, 0.0f, 0.0f);

    _MM_SET_ROUNDING_MODE(_MM_ROUND_DOWN);

    for (size_t i = 0; i < samplesNum; i += 2)
    {
        samples = _mm_setzero_ps();

        for (size_t j = 0; j < ctx.mPhases.size(); ++j)
        {
            // calculate phases, equal to:
            // phases[0] = ctx[j] + freqBuff[i] * 0.5f;
            // phases[1] = ctx[j] + freqBuff[i];
            // phases[2] = ctx[j] + freqBuff[i] + freqBuff[i + 1] * 0.5f;
            // phases[3] = ctx[j] + freqBuff[i] + freqBuff[i + 1];
            phases = _mm_set1_ps(ctx.mPhases[j]);
            phases = _mm_add_ps(phases, _mm_mul_ps(freq0_factor, _mm_set1_ps(freqBuff[i])));
            phases = _mm_add_ps(phases, _mm_mul_ps(freq1_factor, _mm_set1_ps(freqBuff[i + 1])));

            // remove 1.0 if value exceedes 1.0
            __m128 compareMask = _mm_cmpgt_ps(phases, ones);
            __m128 subValue = _mm_and_ps(ones, compareMask);
            phases = _mm_sub_ps(phases, subValue);

            ctx.mPhases[j] = phases.m128_f32[3];

            float ratio = freqBuff[i] * (float)(mRootSize);
            int mipmap = Math::log2_int(ratio);
            if (mipmap > mMipsNum)
                mipmap = mMipsNum;
            if (mipmap < 0)
                mipmap = 0;

            samples = _mm_add_ps(samples, Sample_SSE(mipmap, phases, interpolator));
        }

        output[i] =     ctx.Downsample_SSE(samples.m128_f32);
        output[i + 1] = ctx.Downsample_SSE(samples.m128_f32 + 2);
    }
}

} // namespace mvSynth