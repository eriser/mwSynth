#include "stdafx.hpp"
#include "Wavetable.hpp"
#include "Math.hpp"

namespace mvSynth {


// Elliptic IIR filter coefficients generated with Octave:
// ellip(11, 0.5, 100.0, 0.475)

const float WaveTableContext::a[] =
{
    0.00237215228669364f,
    0.01071603665454480f,
    0.03035791370490471f,
    0.06070430305500466f,
    0.09407806181688412f,
    0.11617732288099189f,
    0.11617732288099196f,
    0.09407806181688420f,
    0.06070430305500469f,
    0.03035791370490472f,
    0.01071603665454481f,
    0.00237215228669364
};

const float WaveTableContext::b[] =
{
    0.0f,
    -3.225641030483127f,
    7.937266825410543,
    -13.361604195910033f,
    18.160728660560746f,
    -19.561560159271554f,
    17.264056989018833f,
    -12.296805901712231f,
    6.989483345163727f,
    -3.043659827605777f,
    0.932448506264840f,
    -0.165901630637920f,
};

WaveTableContext::WaveTableContext()
{
    mPhases.push_back(0.0f);
    Reset();
}

void WaveTableContext::Reset()
{
    for (int i = 0; i < IIR_FILTER_SIZE; ++i)
        mX[i] = mY[i] = 0.0f;
}

void WaveTableContext::Init(size_t numVoices, float* newPhases)
{
    mPhases.clear();
    for (int i = 0; i < numVoices; ++i)
        mPhases.push_back(newPhases[i]);
}

float WaveTableContext::Downsample(float* input)
{
    for (int i = 0; i < 2; ++i)
    {
        for (int j = IIR_FILTER_SIZE - 1; j > 0; --j)
        {
            mX[j] = mX[j - 1];
            mY[j] = mY[j - 1];
        }

        mX[0] = input[i];
        float sum = 0.0f;
        for (int j = 0; j < IIR_FILTER_SIZE; ++j)
            sum += a[j] * mX[j];
        for (int j = 0; j < IIR_FILTER_SIZE; ++j)
            sum -= b[j] * mY[j];
        mY[0] = sum;
    }

    return mY[0];
}



WaveTable::WaveTable()
{
    mData = nullptr;
    mMipsNum = 0;
    mRootSize = 0;
}

WaveTable::~WaveTable()
{
    Release();
}

void WaveTable::Release()
{
    if (mData)
    {
        for (int i = 0; i < mMipsNum; ++i)
        {
            free(mData[i]);
            mData[i] = nullptr;
        }

        free(mData);
        mData = nullptr;
    }

    mMipsNum = 0;
    mRootSize = 0;
}

int WaveTable::LoadData(float* data, int order, const Interpolator& interpolator)
{
    int inSamples = 1 << order; // calculate number of input samples
    mMipsNum = order + 2;
    mRootSize = inSamples * 2; // additional level for upsampled signal

    mData = (float**)malloc(sizeof(float*) * mMipsNum);

    // calculate FFT of the original signal
    float* rootMipmap = (float*)malloc(sizeof(float) * inSamples * 2);
    for (int i = 0; i < inSamples; i++)
    {
        rootMipmap[2 * i] = data[i] / inSamples; // real part
        rootMipmap[2 * i + 1] = 0.0f; // imaginary part
    }
    Math::FFT(rootMipmap, inSamples);

    // upsample by IFFT
    float* tmpData = (float*)malloc(sizeof(float) * inSamples * 4);
    memcpy(tmpData, rootMipmap, sizeof(float) * inSamples);
    memcpy(tmpData + 3 * inSamples, rootMipmap + inSamples, sizeof(float) * inSamples);
    memset(tmpData + inSamples, 0, 2 * sizeof(float) * inSamples);
    Math::IFFT(tmpData, 2 * inSamples);

    // first mipmap - 2x upsampled input data
    mData[0] = (float*)malloc(sizeof(float) * (mRootSize + interpolator.mFilterSize));
    for (int i = 0; i < mRootSize; i++)
        mData[0][i] = tmpData[2 * i]; // take only real part (imaginary should be zero)
    memcpy(&mData[0][mRootSize], &mData[0][0], sizeof(float) * (interpolator.mFilterSize)); // folding

    // second mipmap - input data
    mData[1] = (float*)malloc(sizeof(float) * (mRootSize / 2 + interpolator.mFilterSize));
    memcpy(mData[1], data, sizeof(float) * mRootSize / 2);
    memcpy(&mData[1][mRootSize / 2], &mData[1][0],
           sizeof(float) * (interpolator.mFilterSize)); // folding

    // generate rest of mipmaps
    int samples = inSamples / 2;
    for (int i = 2; i < mMipsNum; i++)
    {
        mData[i] = (float*)malloc(sizeof(float) * (samples + interpolator.mFilterSize));
        memcpy(tmpData, rootMipmap, sizeof(float) * samples);
        memcpy(tmpData + samples, rootMipmap + 2 * inSamples - samples, sizeof(float) * samples);
        Math::IFFT(tmpData, samples);

        for (int j = 0; j < samples; j++)
            mData[i][j] = tmpData[2 * j]; // take only real part

        for (int j = 0; j < interpolator.mFilterSize; ++j) // folding
            mData[i][j + samples] = mData[i][j];

        samples /= 2;
    }

    free(tmpData);
    free(rootMipmap);
    return 0;
}

/*
    Interpolate sample from specified mipmap level.
    Simple FPU version.
    * ratio - sampling ratio: 1.0f - normal, 2.0f - 2x less harmonics, etc...
    * phase - sampling point [0...1.0f)
*/
float WaveTable::Sample_FPU(int mipmap, float phase, const Interpolator& interpolator) const
{
    unsigned int mipmapSize = mRootSize >> mipmap;
    phase *= (float)mipmapSize;
    unsigned int id = (unsigned int)phase;
    const float* pSrc = mData[mipmap];

    float base_phase = phase - floor(phase);

    // scale phase to 0..INTERPOLATOR_FILTER_PHASES-1
    float phase_scaled = base_phase * (float)(interpolator.mFilterPhases);

    // phase interpolation coefficient
    float phase_factor = phase_scaled - floor(phase_scaled);

    // select base phase indicies
    int phase_sel = (int)phase_scaled;

    // interpolation (via FIR filter)
    float sum = 0.0;
    for (int i = 0; i < interpolator.mFilterSize; ++i)
    {
        float x0 = pSrc[id + i];
        float f0 = interpolator.data[phase_sel][i];
        float f1 = interpolator.data[phase_sel + 1][i];
        float f = f0 - phase_factor * (f0 - f1); // lineary interpolate filter coefficients
        sum += x0 * f;
    }

    return sum;
}

void WaveTable::Synth_FPU(size_t samplesNum, const float* freqBuff, WaveTableContext& ctx,
                          const Interpolator& interpolator, float* output) const
{
    // iterato through samples
    for (size_t i = 0; i < samplesNum; i++)
    {
        float freq = freqBuff[i];

        // cap frequency to half sample rate
        if (freq > 0.5f)
            freq = 0.5f;

        float samples[2] = {0.0f, 0.0f};

        // iterate through subvoices
        for (size_t j = 0; j < ctx.mPhases.size(); ++j)
        {
            // calculate phases
            float phaseA = ctx.mPhases[j] + freq * 0.5f;
            float phaseB = ctx.mPhases[j] + freq;
            ctx.mPhases[j] = phaseB;

            if (phaseA > 1.0f)
                phaseA -= 1.0f;
            if (phaseB > 1.0f)
                phaseB -= 1.0f;

            float ratio = freq * (float)(mRootSize);

            // mipmap selection
            // TODO: smooth mipmap blending when ratio is near 1.0
            int mipmap = Math::log2_int(ratio);
            if (mipmap > mMipsNum)
                mipmap = mMipsNum;
            if (mipmap < 0)
                mipmap = 0;

            samples[0] += Sample_FPU(mipmap, phaseA, interpolator);
            samples[1] += Sample_FPU(mipmap, phaseB, interpolator);
        }

        output[i] = ctx.Downsample(samples);
    }
}

} // namespace mvSynth