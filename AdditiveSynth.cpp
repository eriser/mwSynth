#include "stdafx.hpp"
#include "Wavetable.hpp"

int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    std::vector<float> samples;
    samples.reserve(1000000);

    LARGE_INTEGER start, stop, frequency;
    QueryPerformanceFrequency(&frequency);

    mwInterpolator interp;
    interp.Setup(16, 64);

    mwWaveSynthContext ctx;
    ctx.pInterpolator = &interp;

#define WAVE_TABLE_SIZE_POW 11
#define WAVE_TABLE_SIZE (1<<WAVE_TABLE_SIZE_POW)
    float saw[WAVE_TABLE_SIZE];
    for (int i = 0; i < WAVE_TABLE_SIZE; i++)
        //saw[i] = (i<WAVE_TABLE_SIZE/2) ? -0.5f : 0.5f;
        saw[i] = 1.2f * (-0.5f + (float)i / (float)WAVE_TABLE_SIZE);

    mwWaveTable wt;
    wt.LoadData(saw, WAVE_TABLE_SIZE_POW, &interp);

    QueryPerformanceCounter(&start);
    {
        double freq = 0.0005;
        while (freq < 0.5)
        {
            float freqBuffer[2];
            float samplesBuffer[2];
            freqBuffer[0] = (float)(freq *= 1.000001);
            freqBuffer[1] = (float)(freq *= 1.000001);
            wt.Synth_SSE(2, freqBuffer, &ctx, samplesBuffer);

            samples.push_back(samplesBuffer[0]);
            samples.push_back(samplesBuffer[1]);
        }
    }
    QueryPerformanceCounter(&stop);

    float t = (float)(stop.QuadPart - start.QuadPart) / (float)frequency.QuadPart;
    printf("Total time             = %.3f\n", t);
    printf("Megasamples generated  = %.2f\n", (float)samples.size() * 0.000001f);
    printf("Megasamples per second = %.2f\n", (float)samples.size() / t * 0.000001f);
    printf("Max generators num     = %.1f\n", (float)samples.size() / t / 44100.0f);

    FILE* pFile;
    fopen_s(&pFile, "a.wav", "wb");
    fwrite(samples.data(), samples.size() * sizeof(float), 1, pFile);
    fclose(pFile);

    system("pause");
    return 0;
}
