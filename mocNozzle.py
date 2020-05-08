import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
# import seaborn as sns;  sns.set()
import gasdynamics as gas
import time

gamma = 1.4
d2r = np.pi / 180
r2d = 180 / np.pi 
Me = 5
throatHeight = 1
throatRad = throatHeight / 2
numInit = 256
numFullBounce = 1

numTotalLen = numInit * (numFullBounce*(numInit + 1) + 1)
nuE = gas.prandtlMeyerAng(Me, gamma)
# TODO figure out if this is actually what i want for the max angle
thetaE = nuE / 2
dxE = throatRad*np.sin(d2r*thetaE)
dyE = throatRad*(1 - np.cos(d2r*thetaE)) + throatHeight/2
print(dxE, dyE)
ntS = lambda thetaS: r2d*np.arcsin((dxE - (throatHeight/2 + dyE) * np.sqrt(gas.prandtlMeyerM(thetaS)**2 -1)) / throatRad)
# thetaS = thetaE / numInit
newThetaS = abs(np.arcsin(dxE/throatRad * (.5 - dyE/throatHeight))*r2d)
thetaS = newThetaS/2
# newThetaS = ntS(thetaS)
tmp = ntS(newThetaS)
print(thetaE, thetaS, newThetaS)
tS = thetaS
while (newThetaS < tS) or (abs(newThetaS - tS) > .01):
    old = tS
    tS = newThetaS
    tmp2 = tmp
    tmp = ntS(tS) - tS
    newThetaS = tS - (tmp*(tS-old)) / (tmp - tmp2)
    if newThetaS < 0:   newThetaS = tS/2
    if not np.isfinite(newThetaS):  newThetaS = tS
    print(tS, newThetaS, tmp, tmp2)
thetaS = newThetaS
delta = (np.arange(numInit) * (thetaE - thetaS) / (numInit-1) + thetaS).tolist() + [0] * (numTotalLen - numInit)
# nu at the throat is 0 because Mt = 1
# TODO sometimes i get warnings here, check it out
# it seems that the solver will try numbers < 1 and then it throws warnings
# but if you let it recieve the nan itll handle it
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
print(states['R'][numInit-1])

# print(states.__sizeof__())
# print(states.loc[:numInit-1, :], len(states))

# fig, ax = plt.subplots()
# ax.plot(states.x, states.y, '.')
# ax.plot([0, throatRad], [1, 1], '.')
# fig.savefig('init.png')
# plt.close(fig)


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
        finalCharNu = states['R'][numInit-1]
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
print(states, len(states))

# print(states.query(f'(lLine == {numInit-1}) or (rLine == {numInit-1})')['R'])

extra = states.loc[[0, numInit, states.query(f'(lLine == {-2}) or (rLine == {-2})').index[-1]], :].reset_index(drop=True)
extra.loc[:1, 'M'] = 1
extra.loc[:1, 'x'] = 0
extra.loc[0, 'y'] = throatHeight/2
extra.loc[1, 'y'] = 0
extra.loc[2, 'x'] = states['x'].max()
extra.loc[1:2, 'lLine'] = -3
extra.loc[0, 'rLine'] = -3
print(extra)

states = extra.loc[:1, :].append(states.append(extra.loc[2, :], ignore_index=True), ignore_index=True)

fig, ax = plt.subplots(2, 1)
sc = ax[0].scatter(states['x'], states['y'], c=states['M'], alpha=.95, zorder=1, cmap=mpl.cm.gnuplot, marker='.')
print(states)
for i in range(numInit):
    ax[0].plot(states.query(f'(lLine == {i}) or (rLine == {i})')['x'], states.query(f'(lLine == {i}) or (rLine == {i})')['y'], '-', LineWidth=.5, Color='grey', zorder=0)
ax[0].plot(states.query(f'rLine == {-2}')['x'], states.query(f'rLine == {-2}')['y'], 'k--', LineWidth=1, zorder=1)
ax[0].plot(states.query(f'lLine == {-1}')['x'], states.query(f'lLine == {-1}')['y'], 'k', LineWidth=1, zorder=1)
fig.colorbar(sc, ax=ax[0], shrink=.6, aspect=10, label='M', ticks=np.arange(1, int(states['M'].max()+10**-10) + 1))
ax[0].grid()
ax[0].set_aspect('equal')
ax[0].set_xlim(0, max(states['x']))
ax[0].set_ylim(0, max(states['y']))
ax[0].set_xlabel('x/t')
ax[0].set_ylabel('y/t')


triang = mpl.tri.Triangulation(states['x'], states['y'])

con = ax[1].tripcolor(triang, states['M'], shading='gouraud', cmap=mpl.cm.gnuplot)
# for i in range(numInit):
#     ax[1].plot(states.query(f'(lLine == {i}) or (rLine == {i})')['x'], states.query(f'(lLine == {i}) or (rLine == {i})')['y'], '-', LineWidth=.5, Color='grey')
ax[1].plot(states.query(f'rLine == {-2}')['x'], states.query(f'rLine == {-2}')['y'], 'k--', LineWidth=1, zorder=1)
ax[1].plot(states.query(f'lLine == {-1}')['x'], states.query(f'lLine == {-1}')['y'], 'k', LineWidth=1, zorder=1)
fig.colorbar(con, ax=ax[1], shrink=.6, aspect=10, label='M', ticks=np.arange(1, int(states['M'].max()+10**-10) + 1))
ax[1].grid()
ax[1].set_aspect('equal')
ax[1].set_xlim(0, max(states['x']))
ax[1].set_ylim(0, max(states['y']))
ax[1].set_xlabel('x/t')
ax[1].set_ylabel('y/t')

fig.suptitle(f'M = {Me}, Rc = {throatRad}, nwaves = {numInit}')

fig.savefig('final.png', dpi=1200)
plt.show()
plt.close(fig)
