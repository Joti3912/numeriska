import numpy as np
import scipy as sp

m = 1
l = 1 
g = 9.81
ml2 = m*l*l
#theta1 = np.pi/10
#theta2 = np.pi/10
#p1 = 0
#p2 = 0

def dtheta1(t, theta1, theta2, p1, p2):
    taljare = 2*p1 - 3*np.cos(theta1 - theta2)*p2
    namnare = 16 - 9*(np.cos(theta1 - theta2))^2
    return (6/ml2)*(taljare /  namnare) 

def dtheta2(t, theta1, theta2, p1, p2):
    taljare = 8*p2 - 3*np.cos(theta1 - theta2)*p1
    namnare = 16 - 9*(np.cos(theta1 - theta2))^2
    return (6/ml2)*(taljare /  namnare) 

#skriv delta 