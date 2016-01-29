#include "stdafx.hpp"
#include "Interpolator.hpp"
#include "Math.hpp"

namespace mvSynth {

Interpolator::Interpolator()
{
    mFilterPhases = mFilterSize = 0;
    mFilterPhasesF = 0.0f;
    data = nullptr;
}

Interpolator::~Interpolator()
{
    if (data)
    {
        for (int i = 0; i <= mFilterPhases; i++)
            _aligned_free(data[i]);

        free(data);
        data = nullptr;
    }
}


int Interpolator::Setup(int filterSize, int filterPhases)
{
    data = (float**)malloc(sizeof(void*) * (filterPhases + 1));

    for (int i = 0; i <= filterPhases; i++)
        data[i] = (float*)_aligned_malloc(sizeof(float) * filterSize, 16);

    mFilterPhases = filterPhases;
    mFilterSize = filterSize;
    mFilterPhasesF = (float)filterPhases;

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
            double y = Sinc((x + phase) * cutoff);
            double window = BlackmanHarris(((float)i + phase) / (float)filterSize);
            data[filterPhases - 1 - j][i] = (float)(y * window);
        }
    }


    double phase = -1.0 / (float)filterPhases;
    for (int i = 0; i < filterSize; i++)
    {
        double x = (float)(i - filterSize / 2);
        double y = Sinc((x + phase) * cutoff);
        double window = BlackmanHarris(((float)i + phase) / (float)filterSize);
        data[filterPhases][i] = (float)(y * window);
    }

    return 0;
}

} // namespace mvSynth