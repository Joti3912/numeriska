import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

def rhsODE(t,y): 
    return np.sin(t) - y

#Euler Implicit
def solve(fun, tstop, y0, h):
    steps = round(tstop/h)
    tvec = np.linspace(0, tstop, steps)
    i = 0
    y = np.zeros((len(tvec)))
    y[i] = y0
    for t in tvec[0:len(tvec)-1]: 
        g = lambda y_new: y_new - y[i] - h*fun(t+h, y_new)
        y[i+1] = fsolve(g, y[i])
        i += 1
    return tvec, y


ts, y = solve(rhsODE, 10.0, 0, 0.01)
plt.plot(ts, y, label="y(t)")
plt.title('Uppgift 4, Implicit Euler')
plt.xlabel('t')
plt.xlabel('y')
plt.legend()
plt.show()
