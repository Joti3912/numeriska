import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


def rhsODE(t,y): 
    return np.sin(t) - y

def solve(fun, tstop, y0, h):
    t = 0
    y = y0
    for i in [0, tspan]:
        k1 = rhsODE(t,y)
        k2 = rhsODE(t+h, y + h*k1)


y0 = [0.0]
tspan = (0.0, 10.0)
ts = np.arange(0.0, 10.0, 0.01)


plt.plot(result.t, result.y[0], label="y(t)")
plt.xlabel('t')
plt.xlabel('y')
plt.legend()
plt.show()
