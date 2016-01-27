#pragma once

namespace mvSynth {

class Interpolator final
{
    friend class WaveTable;

    float** data;
    int mFilterSize;
    int mFilterPhases;
    float mFilterPhasesF;

public:

    Interpolator();
    ~Interpolator();
    int Setup(int filterSize, int filterPhases);
};

} // namespace mvSynth