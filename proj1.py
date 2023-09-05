import numpy as np
from scipy.integrate import solve_ivp

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
    namnare = 16 - 9*(np.cos(theta1 - theta2))^2
    return (6/ml2)*(taljare /  namnare) 

def dtheta2(theta1, theta2, p1, p2):
    taljare = 8*p2 - 3*np.cos(theta1 - theta2)*p1
    namnare = 16 - 9*(np.cos(theta1 - theta2))^2
    return (6/ml2)*(taljare / namnare) 

def dp1(theta1, theta2, p1, p2):
    dtheta1 = dtheta1(theta1, theta2, p1, p2)
    dtheta2 = dtheta2(theta1, theta2, p1, p2)
    negml2 = (-(1/2))*ml2
    return negml2*(dtheta1*dtheta2*np.sin(theta1 - theta2) + 3*(g/l)*np.sin(theta1))

def dp2(theta1, theta2, p1, p2):
    dtheta1 = dtheta1(theta1, theta2, p1, p2)
    dtheta2 = dtheta2(theta1, theta2, p1, p2)
    negml2 = (-(1/2))*ml2
    return negml2*(-dtheta1*dtheta2*np.sin(theta1 - theta2) + 3*(g/l)*np.sin(theta2))

#skriv delta 