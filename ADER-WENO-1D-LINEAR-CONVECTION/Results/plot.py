import numpy as np 
import matplotlib.pyplot as plt 

x, U = np.loadtxt("sol.csv", delimiter =',', skiprows=1, unpack = True)
x_exact, U_exact = np.loadtxt("sol_exact.csv", delimiter =',', skiprows=1, unpack = True)

plt.figure(1)
plt.plot(x_exact, U_exact, linewidth=1, color='black', label='Exact')
plt.scatter(x, U, color='blue', s=5, label='Numerical')
plt.xlabel(r'$x$')
plt.ylabel(r'$u$')
plt.grid()
plt.legend()
plt.title('Linear Advection')
plt.savefig('U.pdf')

