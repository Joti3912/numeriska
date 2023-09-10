import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


def rhsODE(t,y): 
    return np.sin(t) - y

def solve(fun, tstop, y0, h):
    steps = round(tstop/h)
    tvec = np.linspace(0, tstop, steps)
    i = 0
    y = np.zeros((len(tvec), 1))
    y[i] = y0
    for t in tvec[0:len(tvec)-1]: #OBS! kan ej använda range då rhsODE beror på t
        k1 = rhsODE(t,y[i])
        k2 = rhsODE(t+h, y[i] + h*k1)
        y[i+1, :] = y[i] + (h/2)*(k1 + k2)
        i += 1
    return tvec, y


ts, y = solve(rhsODE, 10.0, [0], 0.01)
plt.plot(ts, y[:,0], label="y(t)")
plt.xlabel('t')
plt.xlabel('y')
plt.legend()
plt.show()
