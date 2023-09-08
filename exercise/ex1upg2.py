import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

#u1 = y(t) -> u1' = y'(t)
#u2 = y'(t) -> u2' = y''(t)
def rhsODE(t, u):
    u1, u2 = u
    du1 = u2
    return [u2, -3*du1 - 2*u1]

u0 = [2,1]
tspan = [0,10]
t_eval = np.arange(0, 10, 0.1)
result = solve_ivp(rhsODE, tspan, u0, t_eval = t_eval)
plt.plot(result.t, result.y[0], label="y(t)")
plt.plot(result.t, result.y[1], label="y'(t)")
plt.xlabel('t')
plt.xlabel('y')
plt.legend()
plt.show()