cdef extern from "2drft.h" namespace "mynamespace":
	cdef cppclass rft2d:
		rft2d( int N1, int N2, double L1, double L2, double Var_Phi, double lx)
		void run()
		void get_last_configuration(double *T)
