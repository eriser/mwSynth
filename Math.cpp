#include "stdafx.hpp"
#include "Math.hpp"

namespace Math {

double PI = 3.14159265358979;

double BlackmanHarris(double x)
{
    return 0.35875 - 0.48828 * cos(2.0 * PI * x) + 0.14128 * cos(4.0 * PI * x) -
           0.01168 * cos(6.0 * PI * x); //Blackman - Harris
}

double Sinc(double phase)
{
    return ((abs(phase) < 1e-2) ?
            (1.0 - PI * PI * phase * phase) :    // Taylor expansion at phase = 0
            (sin(PI * phase) / (PI * phase)));
}

void FFT(float* data, unsigned long nn)
{
    unsigned long n, mmax, m, j, istep, i;
    float wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;

    // reverse-binary reindexing
    n = nn << 1;
    j = 1;
    for (i = 1; i < n; i += 2)
    {
        if (j > i)
        {
            std::swap(data[j - 1], data[i - 1]);
            std::swap(data[j], data[i]);
        }
        m = nn;
        while (m >= 2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    };

    // here begins the Danielson-Lanczos section
    mmax = 2;
    while (n > mmax)
    {
        istep = mmax << 1;
        theta = -2.0f * (float)PI / (float)mmax;
        wtemp = sin(0.5f * theta);
        wpr = -2.0f * wtemp * wtemp;
        wpi = sin(theta);
        wr = 1.0;
        wi = 0.0;
        for (m = 1; m < mmax; m += 2)
        {
            for (i = m; i <= n; i += istep)
            {
                j = i + mmax;
                tempr = wr * data[j - 1] - wi * data[j];
                tempi = wr * data[j] + wi * data[j - 1];

                data[j - 1] = data[i - 1] - tempr;
                data[j] = data[i] - tempi;
                data[i - 1] += tempr;
                data[i] += tempi;
            }
            wtemp = wr;
            wr += wr * wpr - wi * wpi;
            wi += wi * wpr + wtemp * wpi;
        }
        mmax = istep;
    }
}

void IFFT(float* data, unsigned long nn)
{
    // conjugate
    for (unsigned long i = 1; i < 2 * nn; i += 2)
        data[i] = -data[i];

    FFT(data, nn);

    // conjugate
    for (unsigned long i = 1; i < 2 * nn; i += 2)
        data[i] = -data[i];

    //scale
    //float scale = 1.0 / (float)nn;
    //for (int i = 0; i<2*nn; i++)
    //  data[i] *= scale;
}

}