#ifndef CLASSFFT_H
#define CLASSFFT_H

#include<complex>
#include<valarray>
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray; 

class FFT
{
public:
	void fft(CArray& x);
	void ifft(CArray& x);
};
#endif
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


