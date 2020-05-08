import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
# import seaborn as sns;  sns.set()
import gasdynamics as gas
import time
import sys
import argparse as ap

d2r = np.pi / 180
r2d = 180 / np.pi 

# this variable just exists cause theres a couple places i want to have default values but i want to be able to change them from one place
defaultVals = { 'gamma' : 1.4,
                'Me' : 5,
                'throatHeight' : 1,
                'throatRad' : 1 / 2,
                'numInit' : 64,
                'numFullBounce' : 1,
                'show' : True,
                'name' : 'mocTunnel'
                }

class moc:
    def __init__(self, throatRad=defaultVals['throatRad'], throatHeight=defaultVals['throatHeight'], gamma=defaultVals['gamma'], Me=defaultVals['Me'], numInit=defaultVals['numInit'], numFullBounce=defaultVals['numFullBounce'], name=defaultVals['name'], show=defaultVals['show']):
        self.gamma = gamma
        self.Me = Me
        self.throatHeight = throatHeight
        self.throatRad = throatRad
        self.numInit = numInit
        self.numFullBounce = numFullBounce
        self.name = name
        self.show = show

        self.numTotalLen = numInit * (numFullBounce*(numInit + 1) + 1)

    def run(self):
        for num in self.numInit:
            self.calculateTunnel(num)
            print(states.__sizeof__())
            print(states, len(states))
            self.plot()
        if self.show:
            plt.show()
            plt.close('all')

    def plot(self):
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

    def calculateTunnel(self, numInit):
        self.calcFirstPoint(numInit)
        # initialize the state
        delta = (np.arange(numInit) * (self.thetaE - self.thetaS) / (numInit-1) + self.thetaS).tolist() + [0] * (self.numTotalLen - numInit)
        # nu at the throat is 0 because Mt = 1
        # TODO sometimes i get warnings here, check it out
        # it seems that the solver will try numbers < 1 and then it throws warnings
        # but if you let it recieve the nan itll handle it
        M = np.concatenate((gas.prandtlMeyerM(delta[:numInit], self.gamma), np.zeros(self.numTotalLen - numInit)))
        self.states = pd.DataFrame({'M':M, 'delta':delta})

        self.states.loc[:numInit-1, 'nu'] = gas.prandtlMeyerAng(self.states.loc[:numInit-1, 'M'], self.gamma)
        self.states.loc[:numInit-1, 'mu'] = np.arcsin(1 / self.states.loc[:numInit-1, 'M'])*r2d
        # left and right characteristic curves
        self.states.loc[:numInit-1, 'L'] = self.states.loc[:numInit-1, 'nu'] - self.states.loc[:numInit-1, 'delta']
        self.states.loc[:numInit-1, 'R'] = self.states.loc[:numInit-1, 'nu'] + self.states.loc[:numInit-1, 'delta']
        self.states.loc[:numInit-1, 'x'] = self.throatRad * np.sin(self.states.loc[:numInit-1, 'delta']*d2r)
        self.states.loc[:numInit-1, 'y'] = (self.throatHeight/2 + self.throatRad) - self.throatRad * np.cos(self.states.loc[:numInit-1, 'delta']*d2r)
        # lLine and rLine group all the characteristic lines together by the index of the point on the circular section of the nozzle
        # the top wall will have the id of -1 and centerline of -2
        self.states.loc[:numInit-1, 'lLine'] = -1
        self.states.loc[:numInit-1, 'rLine'] = self.states.loc[:numInit-1, :].index
        print(self.states['R'][numInit-1])
         
        # print(self.states.__sizeof__())
        # print(self.states.loc[:numInit-1, :], len(self.states))

        # fig, ax = plt.subplots()
        # ax.plot(self.states.x, self.states.y, '.')
        # ax.plot([0, throatRad], [1, 1], '.')
        # fig.savefig('init.png')
        # plt.close(fig)
        
        i = numInit
        for bounce in range(numFullBounce):
            for curLine in range(numInit):
                i = self.calcCurve(numInit, i, curLine)


    def calcCurve(self, numInit, i, curLine):
        # calc centerline bounce
        r = self.states.loc[i-numInit, :]
        self.states.loc[i, 'lLine'] = curLine
        self.states.loc[i, 'rLine'] = -2
        self.states.loc[i, 'delta'] = 0     # the mid line will have a horizontal flow
        self.states.loc[i, 'nu'] = r['R'] - self.states['delta'][i]
        self.states.loc[i, 'R']  = r['R']
        self.states.loc[i, 'L']  = self.states['nu'][i] - self.states['delta'][i]
        self.states.loc[i, 'M']  = gas.prandtlMeyerM(self.states['nu'][i])
        self.states.loc[i, 'mu'] = np.arcsin(1 / self.states['M'][i])*r2d
        self.states.loc[i, 'x']  = r['x']  - r['y'] / np.tan(d2r*(r['delta'] - r['mu']))
        self.states.loc[i, 'y']  = 0     # centerline defined at 0
        i += 1

        for intLine in np.arange(curLine+1, curLine+numInit) % numInit:
            # calc for each line intersected in order
            # i did this stupidly but im going with it
            # the previous point is on this line so its an R line
            # the other L line is numInit back
            l = self.states.loc[i-1, :]
            r = self.states.loc[i-numInit, :]
            #print(r, l)
            self.states.loc[i, 'lLine'] = curLine
            self.states.loc[i, 'rLine'] = intLine
            self.states.loc[i, 'R']  = r['R']
            self.states.loc[i, 'L']  = l['L']
            self.states.loc[i, 'nu'] = (r['R'] + l['L']) / 2
            self.states.loc[i, 'delta'] = (r['R'] - l['L']) / 2
            self.states.loc[i, 'M']  = gas.prandtlMeyerM(self.states['nu'][i])
            self.states.loc[i, 'mu'] = np.arcsin(1 / self.states['M'][i])*r2d
            self.states.loc[i, 'x']  = (l['y'] - r['y'] + r['x'] * np.tan(d2r*(r['delta']-r['mu'])) - l['x'] * np.tan(d2r*(l['delta']+l['mu']))) / (np.tan(d2r*(r['delta']-r['mu'])) - np.tan(d2r*(l['delta']+l['mu'])))
            self.states.loc[i, 'y']  = r['y'] + (self.states['x'][i] - r['x']) * np.tan(d2r*(r['delta']-r['mu']))
            i += 1
            #print(stuffs)

        # calc top wall
        l = self.states.loc[i-1, :]
        self.states.loc[i, 'lLine'] = -1
        self.states.loc[i, 'rLine'] = curLine
        # this is the state of the last point on the wall
        prevWall = self.states.loc[i-numInit-1, :]
        # we know that the right running characteristic intersects the final mach line
        # therefore the right characteristic is the same as for that final mach line
        # that mach line last hits the centerline where the delta = 0 and then doesnt change
        # so R = finalCharNu = prandtlMeyerAng(Me)
        # so delta = R - nu = finalCharNu - nu
        # where nu is given by all the properties being constant
        finalCharNu = self.states['R'][numInit-1]
        self.states.loc[i, 'delta'] = finalCharNu - l['nu']
        self.states.loc[i, 'nu'] = l['nu']
        self.states.loc[i, 'R']  = self.states['nu'][i] + self.states['delta'][i]
        self.states.loc[i, 'L']  = self.states['nu'][i] - self.states['delta'][i]
        self.states.loc[i, 'M']  = gas.prandtlMeyerM(self.states['nu'][i])
        self.states.loc[i, 'mu'] = np.arcsin(1 / self.states['M'][i])*r2d
        self.states.loc[i, 'x']  = (l['y'] - prevWall['y'] + prevWall['x'] * np.tan(d2r*prevWall['delta']) - l['x'] * np.tan(d2r*(l['delta'] + l['mu']))) / (np.tan(d2r*prevWall['delta']) - np.tan(d2r*(l['delta'] + l['mu'])))
        self.states.loc[i, 'y']  = prevWall['y'] + (self.states['x'][i] - prevWall['x']) * np.tan(d2r*prevWall['delta'])
        i += 1
        #print(stuffs)

        return i

    def calcFirstPoint(self, numInit):
        nuE = gas.prandtlMeyerAng(self.Me, self.gamma)
        # TODO figure out if this is actually what i want for the max angle
        self.thetaE = nuE / 2
        dxE = self.throatRad*np.sin(d2r*self.thetaE)
        dyE = self.throatRad*(1 - np.cos(d2r*self.thetaE)) + self.throatHeight/2
        print(dxE, dyE)
        ntS = lambda thetaS: r2d*np.arcsin((dxE - (self.throatHeight/2 + dyE) * np.sqrt(gas.prandtlMeyerM(thetaS)**2 -1)) / self.throatRad)
        thetaS = self.thetaE / numInit
        newThetaS = abs(np.arcsin(dxE/self.throatRad * (.5 - dyE/self.throatHeight))*r2d)
        # newThetaS = ntS(thetaS)
        tmp = ntS(newThetaS)
        # print(self.thetaE, thetaS, newThetaS)
        tS = thetaS
        while (newThetaS < tS) or (abs(newThetaS - tS) > .01):
            old = tS
            tS = newThetaS
            tmp2 = tmp
            tmp = ntS(tS) - tS
            newThetaS = tS - (tmp*(tS-old)) / (tmp - tmp2)
            if newThetaS < 0:   newThetaS = tS/2
            if not np.isfinite(newThetaS):  newThetaS = tS
            # print(tS, newThetaS, tmp, tmp2)
        self.thetaS = newThetaS
        return self.thetaS


def parse(clargs=sys.argv[1:]):
    # TODO write a real help message so that its usable even without other instructions
    parser = ap.ArgumentParser()
    parser.add_argument('--throatRad',     '-r', default=defaultVals['throatRad'],     type=float,          help='radius of the throat (as a fraction of the throat height')
    parser.add_argument('--throatHeight',  '-H', default=defaultVals['throatHeight'],  type=float,          help='throat height (leave as 1)')
    parser.add_argument('--gamma',         '-g', default=defaultVals['gamma'],         type=float,          help='ratio of specific heats of the gas')
    parser.add_argument('--Mach',          '-M', default=defaultVals['Me'],            type=float,          help='design Mach number for the wind tunnel')
    parser.add_argument('--numInit',       '-i', default=defaultVals['numInit'],       type=int, nargs='+', help='list of number of initial points to make up the circular throat area. ex: 8,16,32,64')
    parser.add_argument('--numFullBounce', '-b', default=defaultVals['numFullBounce'], type=int,            help='number of times the characteristic curves should bounce off the centerline (not working, dont use)')
    parser.add_argument('--name',          '-n', default=defaultVals['name'],          type=str,            help='name of the run (for image output)')
    parser.add_argument('--show',          '-s', default=defaultVals['show'],          action='store_true', help='show interactive plots at the end')

    args = parser.parse_args()
    # TODO actually do some error checking on the values inputted
    return args

def main():
    args = parse()
    tunnel = moc(throatRad=args.throatRad, throatHeight=args.throatHeight, gamma=args.gamma, Me=args.Mach, numInit=args.numInit, numFullBounce=args.numFullBounce, name=args.name, show=args.show)
    tunnel.run()

if __name__ == '__main__':
    main()
