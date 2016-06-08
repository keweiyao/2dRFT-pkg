#include "2drft.h"
#include <iostream>
#include <cmath> 

mynamespace::rft2d::rft2d(int const N1_, int const N2_, double L1, double L2, double Var_Phi, double lx)
: 	N1(N1_), N2(N2_), Nmax(10000)
{
	if (N1%2 !=0 || N2%2 != 0)
	{
		std::cout << "division must be even number." << std::endl;
		return;
	}
	else {half_N1 = N1/2; half_N2 = N2/2;}
	double dV = L1*L2/N1/N2;
	double sigma2_p_D = sqrt(Var_Phi*M_PI*2.0*lx*lx/L1/L2);	
	double N_half_ilength2_p = -M_PI*M_PI*2.0*lx*lx;
	double sqrt2 = sqrt(2.0);
	double iDelta2 = 1.0;;
	//-----------Gaussian-table--------------------------
	std::cout << "Set up gaussian random table, size = " << Nmax << std::endl;
	std::normal_distribution<double> norm(0.0, std::sqrt(iDelta2));
	gaussian_list = new double [Nmax];
	for(int i=0;i<Nmax;i++)
	{
		gaussian_list[i] = norm(generator);
	}
	//---------------Grid-setup--------------------------------
	std::cout << "One time grids (index, argument) set up" << std::endl;
	double p1, p2, arg;
	exp_arg_clip = new bool * [half_N2+1];
	exp_tab1 = new double * [half_N2+1];
	exp_tab2 = new double * [half_N2+1];
	index_tab = new int * [half_N2+1];
	c_index_tab = new int * [half_N2+1];
	N2_low = half_N2+2; N2_high = 0; N1_low = N1+1; N1_high = 0;
	for (k2 = 0; k2 <= half_N2; k2++)
	{
		exp_arg_clip[k2] = new bool [N1];
		exp_tab1[k2] = new double [N1];
		exp_tab2[k2] = new double [N1];
		index_tab[k2] = new int [N1];
		c_index_tab[k2] = new int [N1];
		p2 = (k2 - half_N2)/L2;
		for (k1 = 0; k1 < N1; k1++)
		{
			p1 = (k1 - half_N1)/L1;
			arg = (p1*p1+p2*p2)*N_half_ilength2_p*0.5;
			exp_arg_clip[k2][k1] = true;//(arg >= -8.0);
			if (exp_arg_clip[k2][k1])
			{
				if (k2 < N2_low) N2_low = k2;
				if (k2 > N2_high) N2_high = k2;
				if (k1 < N1_low) N1_low = k1;
				if (k1 > N1_high) N1_high = k1;
			}
			exp_tab1[k2][k1] =  sigma2_p_D*std::exp(arg);
			exp_tab2[k2][k1] =  exp_tab1[k2][k1]/sqrt2;
			index_tab[k2][k1] = k1*N1+k2;
			c_index_tab[k2][k1] = ((N1-k1)%N1)*N1 + ((N2-k2)%N2);
		}
	}
	std::cout << "cut grid: [" << N2_low << ", " << N2_high << "] x [" << N1_low << ", " << N1_high << "]" << std::endl;
	//----------------FFTW-------------------------------------
	A = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N1*N2);
	B = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N1*N2);
	std::cout << "creating FFT plan" << std::endl;
	plan = fftw_plan_dft_2d(N1, N2, A, B, FFTW_BACKWARD, FFTW_MEASURE);
	std::cout << "FFT plan finished" << std::endl;
	for (k2 = 0; k2 < N2; k2++)
	{		
		for (k1 = 0; k1 < N1; k1++)
		{
			A[k1*N1+k2][0] = 0.0;
			A[k1*N1+k2][1] = 0.0;
		}
	}
	//---------------HDF5------------------------------
}

mynamespace::rft2d::~rft2d()
{
	fftw_destroy_plan(plan);
}

void mynamespace::rft2d::run()
{ 
		for (k2 = N2_low; k2 <= N2_high; k2++)
		{
			for (k1 = N1_low; k1 <= N1_high; k1++)
			{
				index = index_tab[k2][k1];
				if ((k1 == 0 || k1 == half_N1) && (k2 == 0 || k2 == half_N2))
				{
					if(exp_arg_clip[k2][k1])
					{
						A[index][0] = exp_tab1[k2][k1]*gaussian_list[generator()%Nmax];
						A[index][1] = 0.0;
					}
				}
				else
				{
					c_index = c_index_tab[k2][k1];
					if(exp_arg_clip[k2][k1])
					{
				  		temp = exp_tab2[k2][k1];
						A[index][0] = temp*gaussian_list[generator()%Nmax];
						A[index][1] = temp*gaussian_list[generator()%Nmax];
						A[c_index][0] = A[index][0];
						A[c_index][1] = -A[index][1];
					}
				}
			}
		}
		fftw_execute(plan);
}


