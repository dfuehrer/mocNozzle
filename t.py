import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib as mpl
import gasdynamics as gas
import time
import sys
import argparse as ap

d2r = np.pi / 180
r2d = 180 / np.pi 

# this variable just exists cause theres a couple places i want to have default values but i want to be able to change them from one place
defaultVals = { 'gamma' : 1.4,
                'R'  : 287,
                'Me' : 5,
                'throatHeight' : 1,
                'throatRad' : 1 / 2,
                'numInit' : [64],
                'numFullBounce' : 1,
                'show' : True,
                'name' : 'mocNozzle',
                'outputDir' : './',
                # TODO make this a reasonable number of things to plot
                'convergeNums' : (2 ** np.arange(3, 9)).tolist(),
                'plotConvergence' : True,
                'plotCenterline'  : True,
                'calcOperation'   : True,
                'gridOverContour' : False
                }

class moc:
    def __init__(self, throatRad=defaultVals['throatRad'], throatHeight=defaultVals['throatHeight'], gamma=defaultVals['gamma'], Me=defaultVals['Me'], numInit=defaultVals['numInit'], numFullBounce=defaultVals['numFullBounce'], name=defaultVals['name'], outputDir=defaultVals['outputDir'], show=defaultVals['show'], gridOverContour=defaultVals['gridOverContour'], plotConvergence=defaultVals['plotConvergence'], plotCenterline=defaultVals['plotCenterline'], calcOperation=defaultVals['calcOperation']):
        self.gamma = gamma
        self.Me = Me
        self.throatHeight = throatHeight
        self.throatRad = throatRad * self.throatHeight
        self.numInit = numInit
        self.numFullBounce = numFullBounce
        self.name = name
        self.outputDir = outputDir
        self.show = show
        self.gridOverContour = gridOverContour
        self.plotConvergence = plotConvergence
        self.plotCenterline  = plotCenterline
        self.calcOperation   = calcOperation
        self.convergeNums    = defaultVals['convergeNums'] if self.plotConvergence else []

        self.theoAreaRat = gas.areaRatio(self.Me, self.gamma)


    def run(self):
        prt = time.time()
        self.ARat = pd.DataFrame({'num':np.unique(self.convergeNums + self.numInit)})
        ind = 0
        for num in self.ARat['num']:
            ptime = time.time()
            self.calculateNozzle(num)
            ntime = time.time()
            print(num, 'point nozzle calc took :', ntime - ptime)
            # print(self.states.__sizeof__())
            # print(self.states, len(self.states))
            # self.ARat.loc[self.ARat.query(f'num == {num}').index, 'rat'] = self.states['y'].max() * 2
            # self.ARat.query(f'num == {num}').loc['rat'] = self.states['y'].max() * 2
            self.ARat.loc[ind, 'rat'] = self.states['y'].max() * 2
            ind += 1
            if num in self.numInit:
                print(num, 'point nozzle calc took :', ntime - ptime)
                ppt = time.time()
                self.plot(num)
                print('plot time:', time.time() - ppt)
                self.states = self.states.assign(pRat = lambda df: (1 + (self.gamma-1)/2 * df.M**2) ** (-self.gamma/(self.gamma-1)),
                                                 TRat = lambda df: (1 + (self.gamma-1)/2 * df.M**2) ** -1)
                if self.calcOperation:
                    pcot = time.time()
                    self.calcOperatingVals()
                    print('calcOperatingVals time: ', time.time() - pcot)
                if self.plotCenterline:
                    ppct = time.time()
                    self.plotCenter()
                    print('calcOperatingVals time: ', time.time() - ppct)

        if self.plotConvergence:
            ppCt = time.time()
            self.ARat['err'] = abs(1 - self.ARat['rat'] / self.theoAreaRat)
            self.plotConverge()
            print('calcOperatingVals time: ', time.time() - ppCt)
        print('Total time: ', time.time() - prt, ' [s]')
        if self.show:
            plt.show()
            plt.close('all')

    def plot(self, numInit):
        split = self.calcExtraPoints(numInit)  # split is the starting index of the points to go at the end of states
        self.states = self.extra.loc[:split-1, :].append(self.states.append(self.extra.loc[split:, :], ignore_index=True), ignore_index=True)

        fig, ax = plt.subplots(2, 1, figsize=(12, 6))
        fig.subplots_adjust(left=0.1, right=1, bottom=0, top=1, hspace=0)
        # TODO try plotting without alpha and see if its faster
        sc = ax[0].scatter(self.states['x'], self.states['y'], c=self.states['M'], alpha=.95, zorder=1, cmap=mpl.cm.gnuplot, marker='.', linewidth=1)
        pmlt = time.time()
        for i in range(numInit):
            ax[0].plot(self.states.query(f'(lLine == {i}) or (rLine == {i})')['x'], self.states.query(f'(lLine == {i}) or (rLine == {i})')['y'], '-', LineWidth=.5, Color='grey', zorder=0)
        print('gen grid lines time:', time.time() - pmlt)
        ax[0].plot(self.states.query(f'rLine == {-2}')['x'], self.states.query(f'rLine == {-2}')['y'], 'k--', LineWidth=1, zorder=1)
        ax[0].plot(self.states.query(f'lLine == {-1}')['x'], self.states.query(f'lLine == {-1}')['y'], 'k', LineWidth=1, zorder=1)
        fig.colorbar(sc, ax=ax[0], shrink=.6, aspect=10, label='M', ticks=np.arange(1, int(self.states['M'].max()+10**-10) + 1))
        ax[0].grid()
        ax[0].set_aspect('equal')
        ax[0].set_xlim(0, max(self.states['x']))
        ax[0].set_ylim(0, max(self.states['y']))
        ax[0].set_xlabel('x/t')
        ax[0].set_ylabel('y/t')


        triang = mpl.tri.Triangulation(self.states['x'], self.states['y'])

        con = ax[1].tripcolor(triang, self.states['M'], shading='gouraud', cmap=mpl.cm.gnuplot)
        if self.gridOverContour:
            for i in range(numInit):
                ax[1].plot(self.states.query(f'(lLine == {i}) or (rLine == {i})')['x'], self.states.query(f'(lLine == {i}) or (rLine == {i})')['y'], '-', LineWidth=.5, Color='grey')
        ax[1].plot(self.states.query(f'rLine == {-2}')['x'], self.states.query(f'rLine == {-2}')['y'], 'k--', LineWidth=1, zorder=1)
        ax[1].plot(self.states.query(f'lLine == {-1}')['x'], self.states.query(f'lLine == {-1}')['y'], 'k', LineWidth=1, zorder=1)
        fig.colorbar(con, ax=ax[1], shrink=.6, aspect=10, label='M', ticks=np.arange(1, int(round(self.states['M'].max())) + 1))
        ax[1].grid()
        ax[1].set_aspect('equal')
        ax[1].set_xlim(0, max(self.states['x']))
        ax[1].set_ylim(0, max(self.states['y']))
        ax[1].set_xlabel('x/t')
        ax[1].set_ylabel('y/t')

        fig.suptitle(f'M = {self.Me}, Rc = {self.throatRad}, nwaves = {numInit}')

        psft = time.time()
        fig.savefig(f'{self.name}_nozzle_M{self.Me}_num{numInit}.png', dpi=1200)
        print('savefig time:', time.time() - pmlt)

    def plotConverge(self):
        fig, ax = plt.subplots()
        ax.plot(self.ARat['num'], self.ARat['err'], '.-')
        ax.set_xlabel('Number of Initial points')
        ax.set_ylabel('Relative Error')
        ax.set_xlim([0, int(self.ARat['num'].max() / 10 + 1) * 10])
        ax.set_ylim([0, int(self.ARat['err'].max() * 2 + 1) / 2])
        ax.grid()
        fig.suptitle('Relative Error in Area Ratio vs Number of Initial Points')
        fig.savefig(self.name + '_converge.png')
        # plt.show()
        # plt.close(fig)

    def plotCenter(self):
        fig, ax1 = plt.subplots()

        lc = 'green'
        rc = 'purple'
        ax1.plot(self.states.query(f'rLine == {-2}')['x'], self.states.query(f'rLine == {-2}')['M'], 'g', label='Mach')
        ax1.set_ylabel('Mach Number', color=lc)
        ax1.tick_params(axis='y', labelcolor=lc)

        ax2 = ax1.twinx()
        ax2.plot(self.states.query(f'rLine == {-2}')['x'], self.states.query(f'rLine == {-2}')['pRat'], 'r', label='pressure ratio')
        ax2.plot(self.states.query(f'rLine == {-2}')['x'], self.states.query(f'rLine == {-2}')['TRat'], 'b', label='Temp ratio')
        ax2.set_ylabel('Total Pressure and Temperature Ratios', color=rc)
        ax2.tick_params(axis='y', labelcolor=rc)

        ax1.set_xlabel('X / Throat Height')
        # ax.set_xlim([0, int(self.ARat['num'].max() / 10 + 1) * 10])
        # ax.set_ylim([0, int(self.ARat['err'].max() / 10 + 1) * 10])
        ax1.grid()
        fig.legend(loc=[.65, .7])
        fig.subplots_adjust(top=0.58)
        fig.suptitle('Mach and Total Pressure and Temperature Ratios vs X')
        fig.tight_layout()
        fig.savefig(self.name + '_MpT.png')
        # plt.show()
        # plt.close(fig)

    def calcInit(self, numInit):
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
         
        # print(self.states.__sizeof__())
        # print(self.states.loc[:numInit-1, :], len(self.states))

        # fig, ax = plt.subplots()
        # ax.plot(self.states.x, self.states.y, '.')
        # ax.plot([0, throatRad], [1, 1], '.')
        # fig.savefig('init.png')
        # plt.close(fig)

    def calculateNozzle(self, numInit):
        self.numTotalLen = numInit * (self.numFullBounce*(numInit + 1) + 1)
        i, cur = self.calcFirstPoint(numInit)
        # self.calcInit(numInit)
        
        # i = numInit
        for bounce in range(self.numFullBounce):
            for curLine in range(cur, numInit):
                i = self.calcCurve(numInit, i, curLine)
            cur = 0

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

    def calcOperatingVals(self):
        tl = self.states['x'].max() - self.states.query('(rLine == -2) and (lLine != -3)')['x'].iloc[-1]
        print(f'Triangle length: {tl:.4f} times the throat height')
        Pcr3 = self.states.iloc[-1].pRat
        print(f'Pcr3 = {Pcr3:.4g}')
        minT0 = 55 / self.states['TRat'].iloc[-1]
        print(f'Minimum total Temperature: {minT0:.4f} K')
        minp0 = .1 / Pcr3
        print(f'Minimum total pressure: {minp0:.4f} [psia]')
        maxp00 = 500 * 6894.76
        minp0 *= 6894.76
        p00 = np.linspace(minp0, maxp00, 1000)
        # T0 = minT0 * (p00 / minp0) ** ((self.gamma-1)/self.gamma)
        eH = .35
        w = .5
        A = eH * w
        areaRat = self.states['y'].max() * 2
        # mdot = A / areaRat * p00 * np.sqrt(self.gamma / defaultVals['R'] / T0) * (2/(self.gamma+1)) ** ((self.gamma+1)/2/(self.gamma-1))
        mdot = A / areaRat * p00 * np.sqrt(self.gamma / defaultVals['R'] / minT0) * (2/(self.gamma+1)) ** ((self.gamma+1)/2/(self.gamma-1))
        mfig, mAx = plt.subplots()
        mAx.plot(p00, mdot)
        mAx.set_xlabel('Operating Total Pressure [Pa]')
        mAx.set_ylabel('Mass Flow Rate [kg/s]')
        mAx.grid()
        mfig.suptitle('Mass Flow Rate vs Operating Total Pressure')
        mfig.savefig(self.name + '_mdot.png')
        # plt.show()
        # plt.close(mfig)

        T00 = 300
        V = 40
        c = lambda p00: self.gamma * np.sqrt(self.gamma * defaultVals['R'] * T00) / V * p00**((1-self.gamma)/2/self.gamma) * A / areaRat * (2/(self.gamma+1)) ** ((self.gamma+1)/2/(self.gamma-1))
        t = 2*self.gamma / (c(p00)*(self.gamma-1)) * (minp0**((1-self.gamma)/2/self.gamma) - p00**((1-self.gamma)/2/self.gamma))
        # t = np.arange(0, 71)
        # p0 = ((self.gamma-1)*c(maxp00)*t / (2*self.gamma) + maxp00**((1-self.gamma)/2/self.gamma)) ** (2*self.gamma/(1-self.gamma))
        tfig, tax = plt.subplots()
        tax.plot(p00, t)
        # ax.plot(t, p0)
        tax.set_xlabel('Operating Total Pressure [Pa]')
        tax.set_ylabel('Total Operation Time [s]')
        tax.grid()
        tfig.suptitle('Total Opteration Time vs Operating Total Pressure')
        tfig.savefig(self.name + '_time.png')
        # plt.show()
        # plt.close(tfig)

        T = minT0 * self.states.iloc[-1]['TRat']
        mu = 1.716e-5 * (T / 273.15)**1.5 * (273.15+110.4) / (T + 110.4)
        rho = p00 * self.states.iloc[-1]['pRat'] / defaultVals['R'] / T
        U = np.sqrt(self.gamma * defaultVals['R'] * T) * self.states.iloc[-1]['M']
        Reprime = rho * U / mu
        Refig, Reax = plt.subplots()
        Reax.plot(p00, Reprime)
        Reax.set_xlabel('Operating Total Pressure [Pa]')
        Reax.set_ylabel('Free-Stream Unit Reynolds Number')
        Reax.grid()
        Refig.suptitle('Unit Reynolds Number vs Operating Total Pressure')
        Refig.savefig(self.name + '_Re.png')
        # plt.show()
        # plt.close(Refig)

        tH = eH / areaRat
        nl = tH * self.states['x'].max()
        print(f'Nozzle length is {nl:.4f} [m]')

        lowArea = gas.areaRatio(0.05, self.gamma) * tH*w
        print(f'The area required for M <= 0.05 is {lowArea:.4g} [m^2] ({lowArea*1000000:.4g} [mm^2])')

        ang = 7
        wLen = .5
        wWidth = .2
        topLen = wLen / np.cos(d2r*ang)
        topVect = np.array([-np.sin(d2r*ang), np.cos(d2r*ang)])
        botVect = np.array([0, -1])
        topP = gas.obliqueShock(self.states['M'].iloc[-1], maxp00 * self.states['pRat'].iloc[-1], ang, gamma=self.gamma)[1]
        botP = gas.obliqueShock(self.states['M'].iloc[-1], maxp00 * self.states['pRat'].iloc[-1], 0,   gamma=self.gamma)[1]
        F = -wWidth * (topLen * topP * topVect + wLen * botP * botVect)
        print(f'The force on the wedge is {abs(F[0]):.4f} [N] to the {"right" if F[0] >= 0 else "left"} and {abs(F[1]):.4f} [N] {"upward" if F[1] >= 0 else "downward"}')


    def calcFirstPoint(self, numInit):
        # TODO what i really want is
        # 1 for this to be more efficient so i dont have to wait for this part to calculate and 
        # 2 for this to use the calcCurve function to calculate the curve for the first line and find the right starting angle to make the first point on the wall within about 2 steps out of the last point of the circle
        # so i should TODO test if using the fsolve function is faster than my newtons or secant method
        nuE = gas.prandtlMeyerAng(self.Me, self.gamma)
        # TODO figure out if this is actually what i want for the max angle
        self.thetaE = nuE / 2
        dxE = self.throatRad*np.sin(d2r*self.thetaE)


        def firstPntErr(ts):
            self.thetaS = ts
            self.calcInit(numInit)
            i = self.calcCurve(numInit, numInit, 0)
            i = self.calcCurve(numInit, i, 1)
            # print(self.states.loc[:numInit*3+1, ['x', 'y', 'M', 'delta', 'nu', 'mu']])
            x1 = self.states.loc[3*numInit+1, 'x']
            # print('x= ', self.states.loc[2*numInit, 'x'], x1, dxE, (self.states.loc[2*numInit, 'x'] - x1) / (x1 - dxE) - .5)
            return (self.states.loc[2*numInit, 'x'] - x1) / (x1 - dxE) + .45, i
        

        ptime = time.time()
        tol = .05
        newThetaS = self.thetaE / numInit/2
        err, i = firstPntErr(newThetaS)
        if abs(err) > tol:
            tS = self.thetaE / numInit * 2/3
            print(tS, newThetaS)
            tmp = firstPntErr(tS)[0]
            while abs(err) > tol:
                old = tS
                tS = newThetaS
                tmp2 = tmp
                tmp = err
                newThetaS = tS - (tmp*(tS-old)) / (tmp - tmp2)
                if newThetaS < 0:   newThetaS = tS/2
                if not np.isfinite(newThetaS):  newThetaS = tS
                err, i = firstPntErr(newThetaS)
                print(tS, newThetaS, err, tmp)
        print('first point time: ', time.time() - ptime)
        self.thetaS = newThetaS
        # print(self.thetaS)

        

        return i, 2

    def calcExtraPoints(self, numInit):
        # returns the starteing index of the points that should go at the end of states
        self.extra = self.states.loc[[0, numInit, self.states.query(f'(lLine == {-2}) or (rLine == {-2})').index[-1]], :].reset_index(drop=True)
        self.extra.loc[:1, 'M'] = 1
        self.extra.loc[:1, 'x'] = 0
        self.extra.loc[0, 'y'] = self.throatHeight/2
        self.extra.loc[1, 'y'] = 0
        self.extra.loc[2, 'x'] = self.states['x'].max()
        self.extra.loc[1:2, 'lLine'] = -3
        self.extra.loc[0, 'rLine'] = -3
        # print(self.extra)
        return 2


def parse(clargs=sys.argv[1:]):
    storeTF = lambda choice: 'store_false' if defaultVals[choice] else 'store_true'
    # TODO write a real help message so that its usable even without other instructions
    parser = ap.ArgumentParser()
    parser.add_argument('--throatRad',         '-r', default=defaultVals['throatRad'],     type=float,          help='radius of the throat (as a fraction of the throat height')
    parser.add_argument('--throatHeight',      '-H', default=defaultVals['throatHeight'],  type=float,          help='throat height (leave as 1)')
    parser.add_argument('--gamma',             '-g', default=defaultVals['gamma'],         type=float,          help='ratio of specific heats of the gas')
    parser.add_argument('--Mach',              '-M', default=defaultVals['Me'],            type=float,          help='design Mach number for the wind nozzle')
    parser.add_argument('--numInit',           '-i', default=defaultVals['numInit'],       type=int, nargs='+', help='list of number of initial points to make up the circular throat area. ex: 8 16 32 64')
    parser.add_argument('--numFullBounce',     '-b', default=defaultVals['numFullBounce'], type=int,            help='number of times the characteristic curves should bounce off the centerline (not implemented, do not use)')
    parser.add_argument('--name',              '-n', default=defaultVals['name'],          type=str,            help='name of the run (for image output)')
    parser.add_argument('--outputDir',         '-o', default=defaultVals['outputDir'],     type=str,            help='name of the directory to output the plot images to')
    parser.add_argument('--noshow',            '-s', action=storeTF('show'),                                    help='toggle off showing interactive plots at the end')
    parser.add_argument('--gridOverContour',   '-G', action=storeTF('gridOverContour'),                         help='show characteristic curves over contour plot')
    parser.add_argument('--noplotConvergence', '-C', action=storeTF('plotConvergence'),                         help='toggle off plotting the convergence WRT number of initial waves')
    parser.add_argument('--noplotCenterline',  '-c', action=storeTF('plotCenterline'),                          help='toggle off plotting Mach and total pressure and temperature ratios on the centerline')
    parser.add_argument('--nocalcOperation',   '-O', action=storeTF('calcOperation'),                           help='toggle off calculation the answers to the opterations questions for part 3')

    args = parser.parse_args()
    # TODO actually do some error checking on the values inputted
    return args

def main():
    args = parse()
    nozzle = moc(throatRad=args.throatRad, throatHeight=args.throatHeight, gamma=args.gamma, Me=args.Mach, numInit=args.numInit, numFullBounce=args.numFullBounce, name=args.name, outputDir=args.outputDir, show=args.noshow, gridOverContour=args.gridOverContour, plotConvergence=args.noplotConvergence, plotCenterline=args.noplotCenterline, calcOperation=args.nocalcOperation)
    nozzle.run()

if __name__ == '__main__':
    main()
