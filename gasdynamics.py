#!/usr/bin/env python3
import numpy as np
import scipy.optimize as opt
import time
import warnings
# might just move opt into the prandtlMeyerM cause it might be the only place its used but that would also slow things down in the middle instead of at the beginning
r2d = 180 / np.pi
d2r = np.pi / 180
def prandtlMeyerAng(M, gamma=1.4):
    # i have a hatred for try except statements cause id rather know whats happening and use if statements specifically but i have no way of knowing if it will be iterable so this is how i have to do it i guess
    try:
        iter(M)
        M = np.array(M)
        # if any(M < 1):  print(M[M < 1])#;  M[M < 1] = 1
    except:
        pass
        # if M < 1:  M = 1
    c = np.sqrt((gamma+1) / (gamma-1))
    return c * np.arctan(np.sqrt(M**2 - 1) / c)*r2d - np.arccos(1 / M)*r2d

def prandtlMeyerM(nu, gamma=1.4):
    # solve for M where prandtlMeyerAng(M) = nu based on the PMM function below
    # as far as i can tell it works fairly well
    # turn M into an array so PMM doesnt die if its a list
    nu = np.array(nu)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        M = opt.fsolve(lambda M: prandtlMeyerAng(M) - nu, PMM(nu, gamma))
    return M if len(M) > 1 else M[0]    # this should be fine cause i think fsolve always outputs an array

def PMM(nu, gamma=1.4):
    # PMM outputs an approcimation of the prandtl-Meyer angle
    # its a bunch of bs
    # first I took the tangent of both sides of the PM function so i could use the sum of angles identity and turn it into soemthing
    # then i noticed that the the top could be simplified and one of the terms could be ommited without a lot of error
    # then i saw that the result was a function with a low value untill an assymtote
    # then i saw that solving for the location of the assymtote as a function of c resulted in a radical function
    # then i linearized it and found the fit, given by a
    # then i found an approximat function for the assymtote as b*g/(a*gamma(a-c)) where b is sqrt(M^2-1)
    # lastly i solved this for M, getting results that are too low in front of the assymtote and too high behind
    c = np.sqrt((gamma+1)/(gamma-1))
    a = 1.9127285932696898/(c-2) + 2.489
    M =  np.sqrt((a*(1 - a/(a + c*gamma*np.tan(d2r * nu))))**2 + 1)
    return M + gamma/5.6 * (1 - np.exp(-2*(M-1)) - (2*M/a) / (1 + np.exp(-a/2*(M-a))))

def areaRatio(M, gamma=1.4):
    return (2/(gamma+1) * (1 + (gamma-1)/2 * M**2))**((gamma+1)/2/(gamma-1)) / M

def dbmBeta(M, delta, gamma=1.4):
    betaEq = [1, -((M**2 + 2) / M**2 + gamma * np.sin(d2r*delta)**2), ((2*M**2 + 1) / M**4 + ((gamma+1)**2 / 4 + (gamma-1) / M**2) * np.sin(d2r*delta)**2), -np.cos(d2r*delta)**2 / M**4]
    Bs = np.roots(betaEq)
    B = np.arcsin(np.sqrt(Bs[Bs > min(Bs)]))*r2d
    return min(B), max(B)

def dbmDelta(M, Beta, gamma=1.4):
    return np.arctan2(M**2 * np.sin(d2r*2*Beta) - 2/np.tan(d2r*Beta), M**2 * (gamma + np.cos(d2r*2*Beta)) + 2)*r2d

def dbmMach(delta, Beta, gamma=1.4):
    return np.sqrt(2 * (np.cot(Beta*d2r) + np.tan(delta*d2r)) / (np.sin(2*Beta*d2r) - np.tan(delta*d2r) * (gamma + cos(2*Beta*d2r))))

def p2Byp1(M, gamma=1.4):
    return (2*gamma * M**2 - (gamma - 1)) / (gamma + 1)

def pByp0(M, gamma=1.4):
    return (1 + (gamma-1)/2 * M**2) ** (-gamma / (gamma - 1))

def normShockM2(M1, gamma=1.4):
    return np.sqrt(((gamma - 1) * M1**2 + 2) / (2*gamma * M1**2 - (gamma - 1)))

def obliqueShock(M1, p1, delta, gamma=1.4):
    # calculate Beta from delta and M1, then normal component of M1 from beta
    Beta = dbmBeta(M1, delta, gamma)[0]
    M1n = M1 * np.sin(d2r*Beta)
    # claculate M2 from M1n, Beta, and delta
    M2 = normShockM2(M1n, gamma) / np.sin(d2r*(Beta - delta))
    # calculate p2 from M1n and p1
    p2 = p2Byp1(M1n, gamma) * p1
    return M2, p2

