import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import gasdynamics as gas

gamma = 1.4
d2r = np.pi / 180
r2d = 180 / np.pi 
Me = 5
throatHeight = 1
throatRad = throatHeight / 2
numInit = 3
numFullBounce = 1

numTotalLen = numInit * (numFullBounce*(numInit + 1) + 1)
nuE = gas.prandtlMeyerAng(Me, gamma)
# TODO figure out if this is actually what i want for the max angle
thetaE = nuE / 2
print(thetaE)
delta = (np.arange(1, numInit+1) * thetaE / numInit).tolist() + [0] * (numTotalLen - numInit)
# nu at the throat is 0 because Mt = 0
M = gas.prandtlMeyerM(delta, gamma)
states = pd.DataFrame({'M':M, 'delta':delta})

states.loc[:numInit, 'nu'] = gas.prandtlMeyerAng(states.loc[:numInit, 'M'], gamma)
states.loc[:numInit, 'mu'] = np.arcsin(1 / states.loc[:numInit, 'M'])*r2d
states.loc[:numInit, 'L'] = states.loc[:numInit, 'nu'] - states.loc[:numInit, 'delta']
states.loc[:numInit, 'R'] = states.loc[:numInit, 'nu'] + states.loc[:numInit, 'delta']
states.loc[:numInit, 'x'] = throatRad * np.sin(states['delta'][:numInit]*d2r)
states.loc[:numInit, 'y'] = (throatHeight/2 + throatRad) - throatRad * np.cos(states['delta'][:numInit]*d2r)

print(states)


