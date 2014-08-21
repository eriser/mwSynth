#pragma once

namespace Math {

extern double PI;

double BlackmanHarris(double x);
double Sinc(double phase);

void FFT(float* data, unsigned long nn);
void IFFT(float* data, unsigned long nn);


inline int log2_int(float x)
{
    unsigned int ix = (unsigned int&)x;
    unsigned int exp = (ix >> 23) & 0xFF;
    return int(exp) - 127;
}

} // namespace mwMath
