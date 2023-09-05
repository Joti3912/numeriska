import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

m = 1
l = 1 
g = 9.81
ml2 = m*l*l
#theta1 = np.pi/10
#theta2 = np.pi/10
#p1 = 0
#p2 = 0

def dtheta1(theta1, theta2, p1, p2):
    taljare = 2*p1 - 3*np.cos(theta1 - theta2)*p2
    namnare = 16 - 9*(np.cos(theta1 - theta2))**2
    return (6/ml2)*(taljare /  namnare) 

def dtheta2(theta1, theta2, p1, p2):
    taljare = 8*p2 - 3*np.cos(theta1 - theta2)*p1
    namnare = 16 - 9*(np.cos(theta1 - theta2))**2
    return (6/ml2)*(taljare / namnare) 

def dp1(theta1, theta2, p1, p2):
    x = dtheta1(theta1, theta2, p1, p2)
    y = dtheta2(theta1, theta2, p1, p2)
    negml2 = (-(1/2))*ml2
    return negml2*(x*y*np.sin(theta1 - theta2) + 3*(g/l)*np.sin(theta1))

def dp2(theta1, theta2, p1, p2):
    x = dtheta1(theta1, theta2, p1, p2)
    y = dtheta2(theta1, theta2, p1, p2)
    negml2 = (-(1/2))*ml2
    return negml2*(-x*y*np.sin(theta1 - theta2) + (g/l)*np.sin(theta2))

def f(t, y):
    theta1, theta2, p1, p2 = y
    dtheta1_dt = dtheta1(theta1, theta2, p1, p2)
    dtheta2_dt = dtheta2(theta1, theta2, p1, p2)
    dp1_dt = dp1(theta1, theta2, p1, p2)
    dp2_dt = dp2(theta1, theta2, p1, p2)
    return [dtheta1_dt, dtheta2_dt, dp1_dt, dp2_dt]

    
result = solve_ivp(f, (0,10), (np.pi/10, np.pi/10, 0, 0))
print(result)
plt.plot(result.t, result.y[0], label="theta1(t)")
plt.plot(result.t, result.y[1], label="theta2(t)")
plt.plot(result.t, result.y[2], label="p1(t)")
plt.plot(result.t, result.y[3], label="p2(t)")
plt.show()