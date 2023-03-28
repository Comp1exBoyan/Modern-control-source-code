// __2020_8_3.cpp : It is a executable cpp file 
// Implenment of fast Fourier transform in cpp
#define _CRT_SECURE_NO_WARNINGS  
#include <iostream>
#include <string>
#include <cmath>

#define FFT_Size 8
#define FFT_PI 3.1415
#define POW log(FFT_Size) / log(2)

namespace FFTToolbox {
	class ComplexNum {

	private:
		double Real;
		double Image;
	public:
		ComplexNum(double r = 0.0, double v = 0.0)
			: Real(r), Image(v) {};

		ComplexNum operator + (ComplexNum& Num)
		{
			ComplexNum v1;
			v1.Real = this->Real + Num.Real;
			v1.Image = this->Image + Num.Image;
			return v1;

		}
		ComplexNum operator - (ComplexNum& Num)
		{
			ComplexNum v2;
			v2.Real = this->Real - Num.Real;
			v2.Image = this->Image - Num.Image;
			return v2;

		}
		ComplexNum operator * (ComplexNum& Num)
		{
			ComplexNum v3;
			v3.Real = this->Real * Num.Real - this->Image * Num.Image;
			v3.Image = this->Real * Num.Image + this->Image * Num.Real;
			return v3;

		}
		friend std::ostream& operator<< (std::ostream& out, const ComplexNum& c);
		friend void InitWn(ComplexNum* W, int n);
		friend void InitArray(ComplexNum* W, int n);

	};
	std::ostream& operator<< (std::ostream& out, const ComplexNum& c) {
		out << c.Real << " " << c.Image << std::endl;
		return out;
	}
	 void InitWn(ComplexNum* W, int n)
	{
		for (int counter = 0; counter < n; counter++)
		{
			W[counter].Real = cos(2 * (FFT_PI / n )* counter);
			W[counter].Image = -1 * sin(2 * (FFT_PI / n) * counter);
			 std::cout << W[counter] << std::endl;
		}


	}
	 void InitArray(ComplexNum* W, int n)
	 {
		 for (int counter = 0; counter < n; counter++)
		 {
			 W[counter].Real = cos(counter+1);
			 W[counter].Image = sin(counter + 1.2);
			 std::cout << W[counter] << std::endl; 
		 }
	 }
	 void FFTCalculator(ComplexNum* IN,ComplexNum *W )
	{
		int L = 0;
		ComplexNum ResultTop, ResultBottom, ResultTmp;
		long int Temp[10] = { 0 };
		long int Out[FFT_Size] = { 0 };
		long int q = 0;
		for (int i = 0; i < FFT_Size; i++)
		{
			int Reg = i;
			for (int j = 0; j < POW; j++)
			{
				Temp[j] = Reg % 2;
				Reg /= 2;
				q = Temp[j] * pow(2, POW - j - 1);
				Out[i] += q;

			}
		}
		
		for (int m = 0; m < POW; m++)
		{
			L = 1 << m;
			for (int j = 0; j < FFT_Size; j++)
			{
				for (int k = 0; k < L; k++)
				{
					ResultTmp = IN[Out[j] + k + L] * W[FFT_Size * k / 2 / L];
					ResultTop = IN[Out[j] + k] + ResultTmp;
					ResultBottom = IN[Out[j] + k] - ResultTmp;
					IN[Out[j] + k] = ResultTop;
					IN[Out[j] + k + L] = ResultBottom;

				 }



			}

		}
		
	}

}


int main()
{
	FFTToolbox::ComplexNum c1[8] = { 0 }, c2[8] = {0};
	FFTToolbox::InitWn(c2,FFT_Size);
	//FFTToolbox::InitArray(c1, FFT_Size);
	//FFTToolbox::FFTCalculator(c1, c2);
	/*
	  complex operation test
	double real, image; 
	std::cin >> real >> image;
	ComplexNum c1(real, image);
	std::cin >> real >> image;
	ComplexNum c2(real, image);
	ComplexNum c3 = c1 + c2;
	std :: cout << c3;
	*******************
          passed         */


}
