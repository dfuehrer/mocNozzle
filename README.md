## mocNozzle
 * calculate wind tunnel nozzle contour
 * calculate some operation specifications

## Usage
 * Download folder
 * run packages.bat file to install the necessary python libraries with pip
   * depends on numpy, pandas, scipy, matplotlib 
 * run mocNozzle.py from command line with desired arguments
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
   * radius of the throat (as a fraction of the throat height
 * `--throatHeight THROATHEIGHT, -H THROATHEIGHT`
   * throat height (leave as 1)
 * `--gamma GAMMA, -g GAMMA`
   * ratio of specific heats of the gas
 * `--Mach MACH, -M MACH`  
   * design Mach number for the wind nozzle
 * `--numInit NUMINIT [NUMINIT ...], -i NUMINIT [NUMINIT ...]`
   * list of number of initial points to make up the circular throat area. ex: 8 16 32 64
 * `--numFullBounce NUMFULLBOUNCE, -b NUMFULLBOUNCE`
   * number of times the characteristic curves should bounce off the centerline (not implemented, do not use)
 * `--name NAME, -n NAME`  
   * name of the run (for image output)
 * `--outputDir OUTPUTDIR, -o OUTPUTDIR`
   * name of the directory to output the plot images to
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

 ![final graph](/plots/mocNozzle_nozzle_M5_num512_rot.png)
