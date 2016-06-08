cimport cpy_2drft
cimport numpy as np
import numpy as np
import sys

cdef class py_2drft:
	cdef cpy_2drft.rft2d * thisptr
	cdef int N1_, N2_
	def __cinit__(self, int N1, int N2, double L1, double L2, double Var, double hiL2):
		if (N1%2 != 0) or (N1%2 != 0):
			print "Please use even grids. Stopped!"
			sys.exit(-1)
			return
		self.N1_ = N1
		self.N2_ = N2
		self.thisptr = new cpy_2drft.rft2d(N1, N2, L1, L2, Var, hiL2)
	def rft_run(self):
		self.thisptr.run()
		
		cdef np.ndarray[np.double_t, ndim=1, mode="c"] XT
		XT = np.ascontiguousarray(np.zeros(self.N1_*self.N2_,), dtype=np.double)

		self.thisptr.get_last_configuration(<double*> XT.data)

		return XT.reshape([self.N1_, self.N2_])
