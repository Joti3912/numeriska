import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.animation as anim

h_ivp = 0.01
h_euler = 0.001
h_runge = 0.01
m = 1
l = 1 
g = 9.81
ml2 = m*l*l

#vinkelhastighet theta1 (övre pendel)
def dtheta1(theta1, theta2, p1, p2):
    taljare = 2*p1 - 3*np.cos(theta1 - theta2)*p2
    namnare = 16 - 9*(np.cos(theta1 - theta2))**2
    return (6/ml2)*(taljare /  namnare) 
#vinkelhastighet theta2 (nedre pendel)
def dtheta2(theta1, theta2, p1, p2):
    taljare = 8*p2 - 3*np.cos(theta1 - theta2)*p1
    namnare = 16 - 9*(np.cos(theta1 - theta2))**2
    return (6/ml2)*(taljare / namnare) 

#förändring av rörelsemängden p1
def dp1(theta1, theta2, p1, p2):
    x = dtheta1(theta1, theta2, p1, p2)
    y = dtheta2(theta1, theta2, p1, p2)
    negml2 = (-(1/2))*ml2
    return negml2*(x*y*np.sin(theta1 - theta2) + 3*(g/l)*np.sin(theta1))

#förändring av rörelsemängden p2
def dp2(theta1, theta2, p1, p2):
    x = dtheta1(theta1, theta2, p1, p2)
    y = dtheta2(theta1, theta2, p1, p2)
    negml2 = (-(1/2))*ml2
    return negml2*(-x*y*np.sin(theta1 - theta2) + (g/l)*np.sin(theta2))

#system av ODE:er som beskriver dubbelpendeln
def ode(t, y):
    theta1, theta2, p1, p2 = y
    dtheta1_dt = dtheta1(theta1, theta2, p1, p2)
    dtheta2_dt = dtheta2(theta1, theta2, p1, p2)
    dp1_dt = dp1(theta1, theta2, p1, p2)
    dp2_dt = dp2(theta1, theta2, p1, p2)
    return np.array([dtheta1_dt, dtheta2_dt, dp1_dt, dp2_dt])

t_eval = np.arange(0.0, 10.0, h_ivp)
result1 = solve_ivp(ode, (0,10), (np.pi/10, np.pi/10, 0, 0), t_eval=t_eval)
plt.figure("Figure 1: solve_ivp, h = {}" .format(h_ivp))
plt.title("Solve_ivp, h = {}" .format(h_ivp))
plt.plot(result1.t, result1.y[0], label="theta1(t)")
plt.plot(result1.t, result1.y[1], label="theta2(t)")
plt.xlabel('t')
plt.ylabel('theta')
plt.legend()
plt.show()

#Beräknar y för varje steg dt inom tspan med eulers metod på funktion f(t,y).
def euler(f, tspan, y0, dt):
    steps = round((tspan)/dt)               #antalet tidpunkter
    tvec = np.linspace(0, tspan, steps)     #array med alla tidpunkter
    y = np.zeros((len(tvec),len(y0)))       #tidpunkter x variabler
    i = 0
    y[i, :] = y0                            #set start values   
    for t in range(len(tvec)-1):            #obs: ODE:n beror ej av t så vi kan använda range
        y[i+1,:] = y[i,:]+dt*f(t,y[i,:])
        i += 1
    return tvec , y

t_euler, y_euler = euler(ode, 10, [np.pi/10, np.pi/10, 0, 0], h_euler)

plt.figure("Figure 2: Euler method, h = {}" .format(h_euler))
plt.title("Euler method, h = {}" .format(h_euler))
plt.plot(t_euler, y_euler[:,0], label="theta1(t)")
plt.plot(t_euler, y_euler[:,1], label='theta2(t)')  
plt.xlabel('t')
plt.ylabel('theta')
plt.legend()
plt.show()

#Beräknar y för varje steg dt inom tspan med klassisk runge-kutta (4-steg) på funktionen f(t,y).
def runge_kutta(f, tspan, y0, dt):
    steps = round((tspan)/dt)               #antalet tidpunkter
    tvec = np.linspace(0, tspan, steps)     #array med alla tidpunkter
    y = np.zeros((len(tvec),len(y0)))       #tidpunkter x variabler
    i = 0
    y[i, :] = y0  
    h2 = dt/2
    h6 = dt/6
    for t in range(len(tvec)-1):
        k1 = f(t,y[i,:])
        k2 = f(t + h2, y[i,:]+h2*k1)
        k3 = f(t+h2,y[i,:]+h2*k2)
        k4 = f(t+dt,y[i,:]+dt*k3)
        y[i+1,:] = y[i,:] + h6*(k1+2*k2+2*k3+k4)
        i+=1
    return tvec,y

 
t_runge, y_runge = runge_kutta(ode, 20, [np.pi/2.2, np.pi/2.2, 0, 0], h_runge)
plt.figure("Figure 3: Runge-kutta method, h = {}" .format(h_runge))
plt.title("Runge-kutta, h = {}" .format(h_runge))
plt.plot(t_runge, y_runge[:,0], label="theta1(t)")
plt.plot(t_runge, y_runge[:,1], label='theta2(t)')
plt.xlabel('t')
plt.ylabel('theta')
plt.legend()
plt.show()

# Omvandlar vinklar till x,y koordinater för respektive pendel (längd från global variabel l) 
def theta_to_xy(theta1, theta2, l):
    return  (l*np.sin(theta1),
            -l*np.cos(theta1),
            l*np.sin(theta1) + l*np.sin(theta2),
            -l*np.cos(theta1) - l*np.cos(theta2))

x1, y1, x2, y2 = theta_to_xy(y_runge[:, 0], y_runge[:, 1], l)

def animate(i):
    ln1.set_data([0, x1[2*i], x2[2*i]], [0, y1[2*i], y2[2*i]])  

plt.figure("Figure 4, animated double pendulum with runge-kutta method")
fig, ax = plt.subplots(1,1, figsize=(5,5))
ax.set_facecolor('w')
ax.get_xaxis().set_ticks([])
ax.get_yaxis().set_ticks([])
ln1, = plt.plot([],[], 'ro--', lw=3, markersize=8)
ax.set_ylim(-3,3)
ax.set_xlim(-3,3)
plt.title("theta1, theta2, p1, p2 = pi/2.2, pi/2.2, 0, 0")
ani = anim.FuncAnimation(fig, animate, frames=round(len(t_runge)/2), interval=1)
ani.save('dubpen.gif', writer='pillow', fps=(50))