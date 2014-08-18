#pragma once

#define M_2PI (2.0*3.14159265359)
#define M_PI (3.14159265359)

namespace Math {

double BlackmanHarris(double x);
double Sinc(double phase);

void FFT(float* data, unsigned long nn);
void IFFT(float* data, unsigned long nn);

__forceinline int log2_int(float x)
{
    unsigned int ix = (unsigned int&)x;
    unsigned int exp = (ix >> 23) & 0xFF;
    return int(exp) - 127;
}

} // namespace mwMath
