#include "stdafx.hpp"
#include "Wavetable.hpp"
#include <chrono>

using namespace mvSynth;

// #define DISABLE_FILE_WRITE


int main(int argc, char* argv[])
{
    (void)argc;
    (void)argv;

    std::vector<float> samples;
    samples.reserve(1000000);

    WaveTableContext ctx;
    Interpolator interpolator;
    interpolator.Setup(16, 64);

#define WAVE_TABLE_SIZE_POW 11
#define WAVE_TABLE_SIZE (1<<WAVE_TABLE_SIZE_POW)
    float saw[WAVE_TABLE_SIZE];
    for (int i = 0; i < WAVE_TABLE_SIZE; i++)
        //saw[i] = (i<WAVE_TABLE_SIZE/2) ? -0.5f : 0.5f;
        saw[i] = 1.2f * (-0.5f + (float)i / (float)WAVE_TABLE_SIZE);

    std::cout << "Wavetable size: " << WAVE_TABLE_SIZE << std::endl;

    WaveTable wt;
    auto load_start = std::chrono::system_clock::now();
    wt.LoadData(saw, WAVE_TABLE_SIZE_POW, interpolator);
    auto load_end = std::chrono::system_clock::now();
    auto load_elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(load_end - load_start);
    std::cout << "Wavetable generation: " << (double)load_elapsed.count() * 1.0e-6 << " ms" << std::endl;

    auto start = std::chrono::system_clock::now();
    {
        double freq = 0.0005;
        while (freq < 0.5)
        {
            float freqBuffer[2];
            float samplesBuffer[2];
            freqBuffer[0] = (float)(freq *= 1.00001);
            freqBuffer[1] = (float)(freq *= 1.00001);
            wt.Synth_SSE(2, freqBuffer, ctx, interpolator, samplesBuffer);

            samples.push_back(samplesBuffer[0]);
            samples.push_back(samplesBuffer[1]);
        }
    }
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
    double t = (double)elapsed.count() * 1.0e-9;
    std::cout << "Elapsed time: " << t << std::endl;
    std::cout << "Megasamples generated: " << (float)samples.size() * 1.0e-6 << std::endl;
    std::cout << "Megasamples per second = " << (float)samples.size() / t * 1.0e-6 << std::endl;
    std::cout << "Max generators num     = " << (float)samples.size() / t / 44100.0f << std::endl;

#ifndef DISABLE_FILE_WRITE
    FILE* pFile = fopen("a.wav", "wb");
    fwrite(samples.data(), samples.size() * sizeof(float), 1, pFile);
    fclose(pFile);
#endif

    system("pause");
    return 0;
}
