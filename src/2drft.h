#ifndef RFT2D_H
#define RFT2D_H

#include <random>
#include <vector>
#include <fftw3.h>
namespace mynamespace {
class rft2d
{
	public:
		rft2d( int N1_=32, int N2_=32, double L1=3.2, double L2=3.2, double sigma2_x=1.0, double half_ilength2_x=10.0);
		~rft2d();
		void run();
		void get_last_configuration(double * T)
		{
			for (k1 = 0; k1 < N1; k1++)
			{
				for (k2 = 0; k2 < N2; k2++)
				{
					T[k1*N1+k2] = B[k1*N1+k2][0]*(2.0*(k2%2)-1)*(2.0*(k1%2)-1);
				}
			}
		};

	private:
		//----------One time grid init----------------------------
		int N2_low, N2_high, N1_low, N1_high, k1, k2, half_N1, half_N2;
		int const N1, N2;
		int index, c_index;
		double temp;
		bool ** exp_arg_clip;
		double ** exp_tab1;
		double ** exp_tab2;
		int ** index_tab;
		int ** c_index_tab;
		
		//---------One time gaussian random variable list-------------------
		std::mt19937 generator;
		double * gaussian_list;
		int const Nmax;
		
		//---------FFTW3--------------------------------------------
		fftw_complex *A, *B;
		double ** C;
		
		fftw_plan plan;
};
}

#endif
