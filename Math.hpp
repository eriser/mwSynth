#pragma once

namespace mvSynth {

#define PI 3.14159265358979

/**
 * Blackman-Harris window function.
 */
double BlackmanHarris(double x);

/**
 * Calculate normalized sinc function.
 * sinc(x) = sin(pi * x) / (pi * x)
 */
double Sinc(double phase);

/**
 * Calculate Fast Fourier Transform.
 * @param data Input and output samples (even - real part, odd - imaginary part).
 * @param nn   Number of samples.
 */
void FFT(float* data, unsigned long nn);

/**
 * Calculate inverse Fast Fourier Transform.
 * @param data Input and output samples (even - real part, odd - imaginary part).
 * @param nn   Number of samples.
 */
void IFFT(float* data, unsigned long nn);

/**
 * Fast logarithm base 2.
 * Result is rounded down to the nearest integer.
 */
__forceinline int log2_int(float x)
{
    unsigned int ix = (unsigned int&)x;
    unsigned int exp = (ix >> 23) & 0xFF;
    return int(exp) - 127;
}

/**
 * Fast logarithm base 2.
 */
__forceinline float fast_log2(float val)
{
    int* const exp_ptr = (int*)(&val);
    int x = *exp_ptr;
    const int log_2 = ((x >> 23)) - 128;
    x &= ~(255 << 23);
    x += 127 << 23;
    *exp_ptr = x;
    val = ((-1.0f / 3.0f) * val + 2.0f) * val - (2.0f / 3.0f);
    //val = ((0.16404256f * val - 1.09886529f) * val + 3.1482979f) * val - 1.2134752f;
    return (val + log_2);
}

} // namespace mvSynth
