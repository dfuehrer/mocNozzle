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
numInit = 16
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
# left and right characteristic curves
states.loc[:numInit-1, 'L'] = states.loc[:numInit-1, 'nu'] - states.loc[:numInit-1, 'delta']
states.loc[:numInit-1, 'R'] = states.loc[:numInit-1, 'nu'] + states.loc[:numInit-1, 'delta']
states.loc[:numInit-1, 'x'] = throatRad * np.sin(states.loc[:numInit-1, 'delta']*d2r)
states.loc[:numInit-1, 'y'] = (throatHeight/2 + throatRad) - throatRad * np.cos(states.loc[:numInit-1, 'delta']*d2r)
# lLine and rLine group all the characteristic lines together by the index of the point on the circular section of the nozzle
# the top wall will have the id of -1 and centerline of -2
states.loc[:numInit-1, 'lLine'] = -1
states.loc[:numInit-1, 'rLine'] = states.loc[:numInit-1, :].index

print(states.__sizeof__())
print(states.loc[:numInit-1, :], len(states))
print(states, len(states))

fig, ax = plt.subplots()
ax.plot(states.x, states.y, '.')
ax.plot([0, throatRad], [1, 1], '.')
fig.savefig('init.png')
plt.close(fig)


######## Start of copying from last time ####### 
# remember the goal is to go untill the wall makes the mach angle inf weak i guess
# also remember we start at the top going right and last time was opposite

i = numInit
for bounce in range(numFullBounce):
    for curLine in range(numInit):
        # calc centerline bounce
        r = states.loc[i-numInit, :]
        states.loc[i, 'lLine'] = curLine
        states.loc[i, 'rLine'] = -2
        states.loc[i, 'delta'] = 0     # the mid line will have a horizontal flow
        states.loc[i, 'nu'] = r['R'] - states['delta'][i]
        states.loc[i, 'R']  = r['R']
        states.loc[i, 'L']  = states['nu'][i] - states['delta'][i]
        states.loc[i, 'M']  = gas.prandtlMeyerM(states['nu'][i])
        states.loc[i, 'mu'] = np.arcsin(1 / states['M'][i])*r2d
        states.loc[i, 'x']  = r['x']  - r['y'] / np.tan(d2r*(r['delta'] - r['mu']))
        states.loc[i, 'y']  = 0     # centerline defined at 0
        i += 1

        fig, ax = plt.subplots()
        ax.plot(states.x, states.y, '.')
        ax.plot([0, throatRad], [1, 1], '.')
        # fig.savefig('init.png')
        plt.show()
        plt.close(fig)

        for intLine in np.arange(curLine+1, curLine+numInit) % numInit:
            # calc for each line intersected in order
            # i did this stupidly but im going with it
            # the previous point is on this line so its an R line
            # the other L line is numInit back
            l = states.loc[i-1, :]
            r = states.loc[i-numInit, :]
            #print(r, l)
            states.loc[i, 'lLine'] = curLine
            states.loc[i, 'rLine'] = intLine
            states.loc[i, 'R']  = r['R']
            states.loc[i, 'L']  = l['L']
            states.loc[i, 'nu'] = (r['R'] + l['L']) / 2
            states.loc[i, 'delta'] = (r['R'] - l['L']) / 2
            states.loc[i, 'M']  = gas.prandtlMeyerM(states['nu'][i])
            states.loc[i, 'mu'] = np.arcsin(1 / states['M'][i])*r2d
            states.loc[i, 'x']  = (l['y'] - r['y'] + r['x'] * np.tan(d2r*(r['delta']-r['mu'])) - l['x'] * np.tan(d2r*(l['delta']+l['mu']))) / (np.tan(d2r*(r['delta']-r['mu'])) - np.tan(d2r*(l['delta']+l['mu'])))
            states.loc[i, 'y']  = r['y'] + (states['x'][i] - r['x']) * np.tan(d2r*(r['delta']-r['mu']))
            i += 1
            #print(stuffs)

        # calc top wall
        l = states.loc[i-1, :]
        states.loc[i, 'lLine'] = -1
        states.loc[i, 'rLine'] = curLine
        # this is the state of the last point on the wall
        prevWall = states.loc[i-numInit-1, :]
        # we know that the right running characteristic intersects the final mach line
        # therefore the right characteristic is the same as for that final mach line
        # that mach line last hits the centerline where the delta = 0 and then doesnt change
        # so R = finalCharNu = prandtlMeyerAng(Me)
        # so delta = R - nu = finalCharNu - nu
        # where nu is given by all the properties being constant
        finalCharNu = gas.prandtlMeyerAng(Me)
        states.loc[i, 'delta'] = finalCharNu - l['nu']
        states.loc[i, 'nu'] = l['nu']
        states.loc[i, 'R']  = states['nu'][i] + states['delta'][i]
        states.loc[i, 'L']  = states['nu'][i] - states['delta'][i]
        states.loc[i, 'M']  = gas.prandtlMeyerM(states['nu'][i])
        states.loc[i, 'mu'] = np.arcsin(1 / states['M'][i])*r2d
        states.loc[i, 'x']  = (l['y'] - prevWall['y'] + prevWall['x'] * np.tan(d2r*prevWall['delta']) - l['x'] * np.tan(d2r*(l['delta'] + l['mu']))) / (np.tan(d2r*prevWall['delta']) - np.tan(d2r*(l['delta'] + l['mu'])))
        states.loc[i, 'y']  = prevWall['y'] + (states['x'][i] - prevWall['x']) * np.tan(d2r*prevWall['delta'])
        i += 1
        #print(stuffs)

print(states.__sizeof__())
print(states.loc[:numInit-1, :], len(states))
print(states, len(states))


fig, ax = plt.subplots()
ax.plot(states.x, states.y, '.')
ax.plot([0, throatRad], [1, 1], '.')
for i in range(numInit):
    ax.plot(states.query(f'(lLine == {i}) or (rLine == {i})')['x'], states.query(f'(lLine == {i}) or (rLine == {i})')['y'], 'c-')
ax.plot(states.query(f'(lLine == {-2}) or (rLine == {-2})')['x'], states.query(f'(lLine == {-2}) or (rLine == {-2})')['y'], 'k--')
ax.plot(states.query(f'(lLine == {-1}) or (rLine == {-1})')['x'], states.query(f'(lLine == {-1}) or (rLine == {-1})')['y'], 'k')
ax.grid()
ax.set_aspect('equal')
fig.savefig('final.png')
plt.show()
plt.close(fig)
