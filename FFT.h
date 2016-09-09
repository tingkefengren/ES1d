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


//�����еĳ�Ա����Ա����������ѡȡ�в��˽⣬�����Ա����Ա�����Ĺ�ϵ�������������ǵĹ�ϵ�в��˽⡣
//�����У�public��Ա��������Ϊ�ⲿ�Ľӿڣ���public��Ա�������Ե���private�еĳ�Ա�����еĳ�Ա����Ҳ�����໥���á�
//����private�еĳ�Ա��ѡȡ�������ǳ�������ʱ������������м���������û����Ҫд��private��Ա��privateһ��Ҳ���Բ�д��
