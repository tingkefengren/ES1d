#ifndef CLASSFFT_H
#define CLASSFFT_H

class FFT
{
public:
	void changeOrder(double (**rhok),double (**rhoki),int n,int t);
	void FFT_1D(double (**rhoj),double (**rhoji),double (**rhok),double (**rhoki), int len,int t);
	void IFFT_1D(double (**phik),double (**phiki),double (**phij),double (**phiji),int len,int t);

};
#endif


//�����еĳ�Ա����Ա����������ѡȡ�в��˽⣬�����Ա����Ա�����Ĺ�ϵ�������������ǵĹ�ϵ�в��˽⡣
//�����У�public��Ա��������Ϊ�ⲿ�Ľӿڣ���public��Ա�������Ե���private�еĳ�Ա�����еĳ�Ա����Ҳ�����໥���á�
//����private�еĳ�Ա��ѡȡ�������ǳ�������ʱ������������м���������û����Ҫд��private��Ա��privateһ��Ҳ���Բ�д��