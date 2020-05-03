#!/usr/bin/env python3
import numpy as np
import scipy.optimize as opt
# might just move opt into the prandtlMeyerM cause it might be the only place its used but that would also slow things down in the middle instead of at the beginning
r2d = 180 / np.pi
d2r = np.pi / 180
def prandtlMeyerAng(M, gamma=1.4):
    # i have a hatred for try except statements cause id rather know whats happening and use if statements specifically but i have no way of knowing if it will be iterable so this is how i have to do it i guess
    try:
        iter(M)
        M = np.array(M)
    except:
        pass
    return np.sqrt((gamma+1) / (gamma-1)) * np.arctan(np.sqrt((gamma-1) / (gamma+1) * (M**2 - 1)))*r2d - np.arccos(1 / M)*r2d

def prandtlMeyerM(nu, gamma=1.4):
    # solve for M where prandtlMeyerAng(M) = nu based on the PMM function below
    # as far as i can tell it works fairly well
    # turn M into an array so PMM doesnt die if its a list
    nu = np.array(nu)
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
    return np.sqrt((a*(1 - a/(a + c*gamma*np.tan(d2r * nu))))**2 + 1)


