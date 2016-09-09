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


//对类中的成员及成员函数参数的选取尚不了解，对类成员，成员函数的关系及外界参数与它们的关系尚不了解。
//在类中，public成员函数将作为外部的接口，而public成员函数可以调用private中的成员，类中的成员函数也可以相互调用。
//所以private中的成员的选取，可以是除调用类时输入与输出的中间变量，如果没有想要写的private成员，private一行也可以不写。
