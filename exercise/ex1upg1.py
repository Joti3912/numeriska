import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp



def rhsODE(t,y): 
    return np.sin(t) - y

y0 = [0.0]
tspan = (0.0, 10.0)
ts = np.arange(0.0, 10.0, 0.01)
result = solve_ivp(rhsODE, tspan, y0, t_eval=ts) 
plt.plot(result.t, result.y[0], label="y(t)")
plt.xlabel('t')
plt.xlabel('y')
plt.legend()
plt.show()

