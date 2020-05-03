import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import gasdynamics as gas
import time

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
# it seems that the solver will try numbers < 1 and then it throws warnings
# but if you let it recieve the nan itll handle it
thetaE = nuE / 2
delta = (np.arange(1, numInit+1) * thetaE / numInit).tolist() + [0] * (numTotalLen - numInit)
# nu at the throat is 0 because Mt = 0
# TODO sometimes i get warnings here, check it out
M = np.concatenate((gas.prandtlMeyerM(delta[:numInit], gamma), np.zeros(numTotalLen - numInit)))
states = pd.DataFrame({'M':M, 'delta':delta})

states.loc[:numInit-1, 'nu'] = gas.prandtlMeyerAng(states.loc[:numInit-1, 'M'], gamma)
states.loc[:numInit-1, 'mu'] = np.arcsin(1 / states.loc[:numInit-1, 'M'])*r2d
states.loc[:numInit-1, 'L'] = states.loc[:numInit-1, 'nu'] - states.loc[:numInit-1, 'delta']
states.loc[:numInit-1, 'R'] = states.loc[:numInit-1, 'nu'] + states.loc[:numInit-1, 'delta']
states.loc[:numInit-1, 'x'] = throatRad * np.sin(states.loc[:numInit, 'delta']*d2r)
states.loc[:numInit-1, 'y'] = (throatHeight/2 + throatRad) - throatRad * np.cos(states.loc[:numInit, 'delta']*d2r)

print(states.__sizeof__())
print(states.loc[:numInit-1, :], len(states))

fig, ax = plt.subplots()
ax.plot(states.x, states.y, '.')
ax.plot([0, throatRad], [1, 1], '.')
fig.savefig('init.png')
plt.close(fig)


