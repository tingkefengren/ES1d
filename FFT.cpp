//An fft code from https://tfetimes.com/c-fast-fourier-transform/

#include <complex>
#include <iostream>
#include <valarray>
#include "FFT.h"
#define _USE_MATH_DEFINES

typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

// Cooley–Tukey FFT (in-place)
void FFT::fft(CArray& x)
{
	const size_t N = x.size();
	if (N <= 1) return;

	// divide
	CArray even = x[std::slice(0, N/2, 2)];
	CArray  odd = x[std::slice(1, N/2, 2)];

	// conquer
	fft(even);
	fft(odd);

	// combine
	for (size_t i = 0; i < N/2; ++i)
	{
		Complex t = std::polar(1.0, -2 * M_PI * i / N) * odd[i];
		x[i    ] = even[i] + t;
		x[i+N/2] = even[i] - t;
	}
}

// inverse fft (in-place)
void FFT::ifft(CArray& x)
{
	// conjugate the complex numbers
	x= x.apply(std::conj);

	// forward fft
	fft( x );

	// conjugate the complex numbers again
	x = x.apply(std::conj);

	// scale the numbers
	x /= x.size();
}

