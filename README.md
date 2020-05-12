## mocNozzle
 * Calculate wind tunnel nozzle contour
 * Calculate some operation specifications

## Usage
 * Download folder
 * Make sure you have python and pip installed
 * Run packages.bat file to install the necessary python libraries with pip
   * Depends on numpy, pandas, scipy, matplotlib 
 * Run mocNozzle.py from command line with desired arguments
```
mocNozzle.py [-h] [--throatRad THROATRAD] [--throatHeight THROATHEIGHT] [--gamma GAMMA] [--Mach MACH]
                    [--numInit NUMINIT [NUMINIT ...]] [--numFullBounce NUMFULLBOUNCE] [--name NAME]
                    [--outputDir OUTPUTDIR] [--save] [--noshow] [--gridOverContour] [--noplotConvergence]
                    [--noplotCenterline] [--nocalcOperation]
```
### Arguments
 * `-h, --help`            
   * show this help message and exit
 * `--throatRad THROATRAD, -r THROATRAD`
   * set radius of the throat (as a fraction of the throat height
 * `--throatHeight THROATHEIGHT, -H THROATHEIGHT`
   * set throat height (leave as 1)
 * `--gamma GAMMA, -g GAMMA`
   * set ratio of specific heats of the gas
 * `--Mach MACH, -M MACH`  
   * set design Mach number for the wind nozzle
 * `--numInit NUMINIT [NUMINIT ...], -i NUMINIT [NUMINIT ...]`
   * set list of number of initial points to make up the circular throat area. ex: 8 16 32 64
 * `--numFullBounce NUMFULLBOUNCE, -b NUMFULLBOUNCE`
   * set number of times the characteristic curves should bounce off the centerline (not implemented, do not use)
 * `--name NAME, -n NAME`  
   * set name of the run (for image output)
 * `--outputDir OUTPUTDIR, -o OUTPUTDIR`
   * specify name of the directory to output the plot images to
 * `--save, -S`            
   * save the data to an pickle file (h5 requires another dependency)
 * `--noshow, -s`          
   * toggle off showing interactive plots at the end
 * `--gridOverContour, -G`
   * show characteristic curves over contour plot
 * `--noplotConvergence, -C`
   * toggle off plotting the convergence WRT number of initial waves
 * `--noplotCenterline, -c`
   * toggle off plotting Mach and total pressure and temperature ratios on the centerline
 * `--nocalcOperation, -O`
   * toggle off calculation the answers to the operations questions for part 3
### Defaults
 By default, mocNozzle.py runs with
 * gamma = 1.4
 * R  = 287
 * Me = 5
 * throatHeight = 1
 * throatRad = 1 / 2
 * numInit = [256]
 * numFullBounce = 1
 * save = False
 * show = True
 * name = 'mocNozzle'
 * outputDir = './'
 * convergeNums = [8 16 32 64 128 256]
 * plotConvergence = True
 * plotCenterline  = True
 * calcOperation   = True
 * gridOverContour = False
### Example
 * generate all plots and output without interactive plots
```
python .\mocNozzle.py -s
```
 * run with interactive plots, without convergence, with 64 starting points
```
python .\mocNozzle.py -Ci 64
```

 ![final graph](/plots/mocNozzle_nozzle_M5_num512_rot.png)
