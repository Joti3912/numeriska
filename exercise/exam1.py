import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

#   -------UPGIFT 1--------

def upg1ode(t,y):
    return 1.8*y*(1 - (y/10000))

upg1tspan = [0,8]
upg1t_eval = np.arange(0, 8, 0.1)
result1 = solve_ivp(upg1ode, upg1tspan, [10], t_eval = upg1t_eval)
plt.figure("Upgift 1")
plt.title("Upgift 1")
plt.plot(result1.t, result1.y[0], label="y(t)")
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()



#   -------UPGIFT 2--------

def runge_kutta(fun, tstop, y0, h):
    steps = round((tstop)/h)               
    tvec = np.linspace(0, tstop, steps)     
    y = np.zeros((len(tvec)))
    i = 0
    y[i] = y0
    for t in range(len(tvec)-1):
        k1 = fun(t,y[i])
        k2 = fun(t + (h/2), y[i]+(h/2)*k1)
        k3 = fun(t+(h/2),y[i]+(h/2)*k2)
        k4 = fun(t+h,y[i]+h*k3)
        y[i+1] = y[i] + (h/6)*(k1+2*k2+2*k3+k4)
        i+=1
    return tvec, y

ts_upg2, y_upg2 = runge_kutta(upg1ode, 8, 10, 0.1)
plt.figure("Upgift 2")
plt.title("Upgift 2")
plt.plot(ts_upg2, y_upg2, label="y(t)")
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()


#   --------UPGIFT 3--------
#u1 = y(t) -> u1' = y'(t)
#u2 = y'(t) -> u2' = y''(t)

def up3ode(t, u):
    u1, u2 = u
    du1 = u2
    return [u2, ((-0.1/2)*du1 - 9.81*np.sin(u1))]


up3_u0 = [0,3]
upg3tspan = [0,10]
upg3t_eval = np.arange(0, 10, 0.1)
result3 = solve_ivp(up3ode, upg3tspan, up3_u0, t_eval = upg3t_eval)
plt.figure("Upgift 3")
plt.title("Upgift 3")
plt.plot(result3.t, result3.y[0], label="y(t)")
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()


#  ------UPGIFT 4------- 
def upg4ode(t,y): 
    return (np.e)**(t*(np.sin(y)))

#trapets: y_i+1 = y_i + h/2[f(t_i, y_i) + f(t_i+1, y_i+1)]
def trapets(fun, tstop, y0, h):
    steps = round(tstop/h)
    tvec = np.linspace(0, tstop, steps)
    i = 0
    y = np.zeros((len(tvec)))
    y[i] = y0
    for t in tvec[0:len(tvec)-1]: 
        g = lambda y_new: y_new - y[i] - (h/2)*(fun(t, y[i]) + fun(t+h, y_new))
        y[i+1] = fsolve(g, y[i])[0]
        i += 1
    return tvec, y


ts_upg4, y_upg4 = trapets(upg4ode, 3, 0, 0.01)
plt.figure("Upgift 4")
plt.title('Uppgift 4, Trapets')
plt.plot(ts_upg4, y_upg4, label="y(t)")
plt.xlabel('t')
plt.ylabel('y')
plt.legend()
plt.show()