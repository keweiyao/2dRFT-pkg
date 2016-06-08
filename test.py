import matplotlib.pyplot as plt
import cpy_2drft as rft
import numpy as np
import scipy.signal as sg


N = 200
L = 80
Var_Phi = 1.0
lx = 2.0

f = rft.py_2drft(N, N, L, L, Var_Phi, lx)
x = f.rft_run()
plt.imshow(x, extent = [0, L, 0, L])


mean = np.mean(x)
std = np.std(x)
phi4 = np.mean(x**4)
phi6 = np.mean(x**6)
print "<phi>/std = ", mean/std
print "std = ", std
print "<phi^4>/std^4", phi4/std**4 
print "<phi^6>/std^6", phi6/std**6 
plt.colorbar()
plt.show()

