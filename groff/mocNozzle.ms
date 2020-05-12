.nr FM 1
.PP
.am PP
.nh
..
.de label
.DS C
\\$1
.DE
..
.EQ
delim $$
.EN
.TL
Method of Characteristics Nozzle Design
.AU
Daniel Fuehrer

.NH
Introduction
.PP
The method of characteristics is a simple and effective method for calculating flow characteristics for a supersonic flow, such as modeling the flow within a supersonic nozzle.
This method can be used to find the contour of such a nozzle by specifying that the Mach waves coming off the tunnel wall must be infinitely weak.
This paper focuses on a program that calculates the contour for a nozzle based on the given Mach number using the method of characteristics.

.NH
Methods
.NH 2
Method of Characteristics Nozzle
.PP
The mocNozzle.py script has a class with methods to calculate the necessary wind tunnel data.
It first calculates the flow properties in the nozzle, using them to find its contour.
The nozzle starts at the throat as a section of a circle.
Points on the circle are evenly spaced out up to $delta sub {roman max} = nu sub {roman E}$, with the first point offset from the center of the throat.
The script starts by choosing the first point such that where its path ends on the wall is roughly halfway between the points on the wall on either side, shown in the image below.
.PSPIC ps/zi3.eps
.label "Figure 1.  First Point Choice"
For each point, the script calculates the Mach number, direction of the flow, Prandtl-Meyer angle, Mach angle, left and right running characteristic values, and $x$ and $y$ coordinates.
At each point in the flow, the left and right running characteristics are the same as the preceding points, so the flow direction and Prandtl-Meyer angles are calculated from the equations below.
.EQ (1)
nu = {R + L} over 2
.EN
.EQ (2)
delta = {R - L} over 2
.EN
From $nu$, the Mach number can be found from numerically solving the Prandtl-Meyer equation, and $mu$ can be found from the Mach number with $mu = sin ( 1 smallover M )$.
Lastly the coordinates of the new point are found from the intersection of the Mach lines at each point using
.EQ (3)
x = {y sub L - y sub R + x sub R tan ( delta sub R - mu sub R ) - x sub R tan ( delta sub L + mu sub L )}
over 
{ tan ( delta sub R - mu sub R ) - tan ( delta sub L + mu sub L )}
.EN
.EQ (4)
y = y sub R + ( x - x sub R ) tan ( delta sub R - mu sub R )
.EN
At the centerline, the flow direction is known to be horizontal ($delta = 0$), and the other values follow from there. 

.PP
At the top wall, the shape has to cancel out the incoming characteristic curves.
Since the curves bouncing off the wall will intersect with the final curve of the design Mach number with horizontal flow:
.EQ (5)
R sub {roman wall} = L sub {roman final} = nu sub {roman final}
.EN
From this, the flow angle at the wall can be obtained from $delta = R - nu$, and the other flow characteristics can be found the same way as for other points.
.NH 2
Nozzle Contour
.NH 3
Convergence
.PP
The error in the calculation can be simply measured by comparing the end height to the theoretical area ratio for a nozzle of this design Mach number.
.EQ (6)
A over {A sup *} =
1 over M left ( 2 over { gamma + 1 } left ( 1 + { gamma - 1 } over 2 M sup 2 right ) right ) sup {{ gamma + 1 } over {2( gamma - 1 )}} =
1 over 5 left ( 2 over { 1.4 + 1 } left ( 1 + { 1.4 - 1 } over 2 5 sup 2 right ) right ) sup {{ 1.4 + 1 } over {2( 1.4 - 1 )}} =
25
.EN
To test convergence, the nozzle was calculated for different numbers of initial points, from the $3 sup {roman rd}$ to the $9 sup {roman th}$ power of 2.
Then the relative error in area ratio was plotted versus the number of initial points.
.NH 3
Nozzle Contour
.PP
The contour of the nozzle was plotted with two sub-plots.
The top plot shows the points in their $x$ and $y$ coordinates, bounded by the black line of the nozzle wall.
During the nozzle calculation, an ID was stored for the left and right curves at each point.
Then the characteristic curves were made by plotting the lines connecting points with the same IDs.
Similarly, the centerline and wall had IDs to identify them for plotting.
.PP
The bottom plot is a contour map of the flow within the nozzle.
To plot a contour of the whole nozzle rather than solely the calculated points, two points were added at the top and bottom of the throat and one point at the end was added on the centerline.
The colors in both plots vary based on the Mach number at the point.
.NH 3
Centerline Distributions
.PP
After calculating the nozzle, the total pressure and total temperature ratios were calculated at each point.
Based on the IDs mentioned above, the pressure ratios, temperature ratios, and Mach numbers were selected to just the points on the centerline.
These values were then plotted against their $x$ position all in one plot with separate y-axes to make the values of the ratios clearly visible compared to Mach number.
One benefit of plotting them this way is that the relationship between the Mach number and ratios is clear.
.PP
The $3 sup {roman rd}$ critical pressure ratio $P sub {roman cr3}$ is given by the total pressure ratio at the nozzle exit.
.NH 3
Triangle Length
.PP
The triangle length was calculated as the difference between the maximum $x$ value and the $x$ value of the last point on the centerline.
.NH 2
Tunnel Characteristics
.NH 3
Minimum Total Temperature
.PP
The minimum total temperature was calculated as the temperature of liquidation divided by the minimum temperature ratio.
This gives liquidation at the exit of the nozzle.
.NH 3
Minimum Total Pressure
.PP
The minimum total pressure is given by the vacuum pressure of the test section divided by $P sub {roman cr3}$.
.NH 3
Mass Flow Rate
.PP
The mass flow rate is given by:
.EQ (7)
m dot =
A left ( A sup \(** over A right )
p sub 0 over sqrt T sub 0
sqrt { gamma over R }
left ( 2 over { gamma + 1 } right ) sup 
{{ gamma + 1 } over { 2 ( gamma - 1 ) }}
.EN
Plugging in the given height, width, and $T sub {0 roman min}$ gives an expression for mass flow rate in terms of $p sub 0$.
This relationship was then plotted.
.NH 3
Operation Time
.PP
The mass flow rate is given by:
.EQ (8)
m dot =
{dm} over {dt} =
V left ( { d rho sub 0 } over {dt} right ) =
V over R left ( {d ( p sub 0 smallover T sub 0 ) } over {dt} right ) =
A left ( A sup \(** over A right )
p sub 0 over sqrt T sub 0
sqrt { gamma over R }
left ( 2 over { gamma + 1 } right ) sup 
{{ gamma + 1 } over { 2 ( gamma - 1 ) }}
.EN
From isentropic flows:
.EQ (9)
T sub 0 =
T sub {0 0} left ( p sub 0 over p sub {0 0} right ) sup {{ gamma - 1 } over gamma }
.EN
Solving for $p sub 0 smallover T sub 0$
.EQ (10)
p sub 0 over T sub 0 =
{{p sub {0 0}} sup {{ gamma - 1 } smallover gamma } {p sub 0} sup { 1 smallover gamma }} over T sub {0 0}
.EN
Substitute back into mass flow rate
.EQ (11)
V over R left ( {d ( p sub 0 smallover T sub 0 )} over {dt} right ) =
V over R left ( {p sub {0 0}} sup {{ gamma - 1 } smallover gamma } over T sub {0 0} right ) ^ left ( {d ( {p sub 0} sup { 1 smallover gamma } )} over {dt} right ) =
V over R left ( {p sub {0 0}} sup {{ gamma - 1 } smallover gamma } over T sub {0 0} right ) ^ left ( 1 smallover gamma {p sub 0} sup {{ 1 - gamma } smallover gamma } right ) ^ left ( {d p sub 0} over {dt} right )
.EN
Solve for ${d p sub 0} smallover {dt}$
.EQ (12)
{d p sub 0} over {dt} =
{ gamma sqrt { gamma R T sub {0 0} }}
over V
{p sub {0 0}} sup {{ 1 - gamma } smallover gamma }
A left ( A sup \(** over A right ) ^
left ( 2 over { gamma + 1 } right ) sup 
{{ gamma + 1 } over { 2 ( gamma - 1 ) }}

left ( {p sub 0} sup {{ 3 gamma - 1 } over gamma } right ) 
.EN
Let
.EQ (13)
c =
{ gamma sqrt { gamma R T sub {0 0} }}
over V
{p sub {0 0}} sup {{ 1 - gamma } smallover gamma }
A left ( A sup \(** over A right ) ^
left ( 2 over { gamma + 1 } right ) sup 
{{ gamma + 1 } over { 2 ( gamma - 1 ) }}
.EN
.EQ (14)
{d p sub 0} over {dt} =
-c ^
left ( {p sub 0} sup {{ 3 gamma - 1 } over gamma } right ) 
.EN
Separate and integrate.
Find $C sub 0$ from $p sub 0 ( 0 ) = p sub {0 0}$.
.EQ (15)
left ( { 2 gamma } over { 1 - gamma } right ) 
{p sub 0} sup {{1 - gamma } over { 2 gamma }} =
-ct + C sub 0 = 
-ct + 
left ( { 2 gamma } over { 1 - gamma } right ) ^
{p sub {0 0}} sup {{1 - gamma } over { 2 gamma }}
.EN
Solving for $t$ at $p sub {0 roman min}$
.EQ (16)
t =
left ( { 2 gamma } over { c ( gamma - 1 )} right ) ^
left ( {p sub {0 roman min}} sup {{1 - gamma } over { 2 gamma }} - {p sub {0 0}} sup {{1 - gamma } over { 2 gamma }} right )
.EN
Using this expression, the operation time $t$ was plotted against the operation total pressure $p sub 0$.
.NH 3
Unit Reynolds Number
.PP
The unit Reynolds number is given by the equation below, using Sutherland's law to approximate viscosity with respect to the temperature and isentropic relations for $rho$ and $T$.
.EQ (17)
Re prime =
{ rho U } over mu =
{ rho U } 
over {
mu sub {roman ref} left ( T over T sub {roman ref} right ) sup {3/2} ^
left ( { T sub {roman ref} + S } over {  T + S } right )
} =
{ ( p smallover p sub 0 ) p sub 0 U } 
over {
mu sub {roman ref} R
{( {( T smallover T sub 0 ) T sub {0 roman min}} ) sup {5/2}}
over 
T sub {roman ref} sup {3/2} ^
left ( { T sub {roman ref} + S }
over
{  ( T smallover T sub 0 ) T sub {0 roman min} + S } right )
}
.EN
where $mu sub {roman ref} = 1.716 E -5 ~[ {roman Pa} fat cdot {roman s} ]$, $T sub {roman ref} = 273.15~ roman K$, and $S = 110.4~ roman K$.
A plot of unit Reynolds number versus operating total pressure was made using the above expression.
.NH 3
Tunnel Geometry
.PP
The length of the nozzle was found by multiplying the maximum $x$ value by the throat height.
The throat height was found by dividing the exit height by the area ratio, since the tunnel has a constant width.
.PP
The area ratio needed for a $M <= 0.05$ flow was found using the area ratio equation (7) with that Mach number
Multiplying the area ratio by the width times the throat height gives the area needed.
.NH 3
Model Forces
.PP
The forces on the wedge were found first by calculating the areas of the sides using the given angle, length, and width.
Then normal vectors to the sides were defined based on the sines and cosines of the angle.
After that, the angle was used to find the angle of the oblique shock wave on the top and bottom.
From there, the normal components of the flow to the oblique shocks were used to find the pressures on the top and bottom of the wedge with (18) from normal shock relations.
.EQ (18)
p sub 2 over p sub 1 = {2 gamma M sup 2 - ( gamma - 1)} over { gamma + 1 }
.EN
Lastly, the negatives of the pressures were multiplied by the areas and normal vectors to get the total force on each side, which were added together to get the total force on the wedge.

.NH
Results
.NH 2
Nozzle Contour
.NH 3
Convergence
.PP
Below is the convergence plot.
Qualitatively, the relative error seems to be inversely proportional to the number of initial points.
The error at 512 points is 1.07%, compared to the error of 129.6% with 8 initial points.
.PSPIC ./ps/mocNozzle_converge.eps
.label "Figure 2.  Convergence Plot"
.NH 3
Nozzle Contour
.PP
Below is the plot of the nozzle contour.
As discussed above, the final height is very close to the theoretical 12.5.
The nozzle is quite long, about 75 times the length of the throat height.
Because of this, the semi-circular portion of the nozzle at the throat is indistinguishable.
The Mach contour quickly increases along the x-axis, with the Mach number reaching 3 around 10% of the total length and reaching 5 around 20% of the length.

.KS
.PSPIC ./ps/mocNozzle_nozzle_M5_num512_rot.eps 7.5 9.5
.label "Figure 3.  Nozzle Contour"
.KE

.NH 3
Centerline Distributions
.PP
The centerline distributions of Mach number, total temperature and total pressure are shown in Figure 4 below.
This plot shows how quickly the Mach number rises from 1 to 5, where it stays steady to the end of the length.
Since the total temperature and total pressure ratios are dependent on the Mach number, they in turn fall quickly from their initial values of 0.8333 and 0.5283 to 0.1667 and 0.00189, respectively.
As indicated above, $P sub {cr3} = 0.00189$.
The relationship between the Mach number and total temperature and pressure ratios is obvious with the ratios sharply plateauing as the Mach number does.
From this chart, it is clear that the total pressure ratio drops more quickly than the total temperature ratio.
.PSPIC ./ps/mocNozzle_MpT.eps
.label "Figure 4. $M$, $p smallover p sub 0$, $T smallover T sub 0$ Centerline Distributions"
.NH 3
Triangle Length
.PP
The length of the triangle of uniform flow is 61.8934 times the throat height, which can be seen in the centerline plot above.
If the nozzle had an exit height of 350 [mm] as in next section, the triangle of uniform flow would be about 866.5 [mm] long.
.NH 2
Tunnel Characteristics
.NH 3
Minimum Total Temperature
.PP
The minimum total temperature found was 330 K; total temperatures less than this would result in liquefaction inside the nozzle.
.NH 3
Minimum Total Pressure
.PP
The minimum total pressure for the wind tunnel is 52.9090 [psia].
If the total pressure was lower than this, then the flow would be over-expanded, causing shocks to form at the exit, and the nozzle would not be started.
.NH 3
Mass Flow Rate
.PP
The mass flow rate for a nozzle 500 [mm] wide with an exit height of 350 [mm] operating at the minimum total temperature is plotted below.
The relationship is linear with the operating total pressure, as mandated by the driving expression, ranging from 5.62 to 53.12 $[ roman { kg smallover s } ]$.
.PSPIC ./ps/mocNozzle_mdot.eps
.label "Figure 5.  Mass Flow Rate"
.NH 3
Operation Time
.PP
The possible operation time is plotted against the initial total pressure in the Figure 6.
The relationship is nonlinear, with increases in initial total pressure bringing diminishing returns in possible operation time.
The operation time given by this graph is less than it would be if the total temperature was assumed to remain constant.
The maximum operation time for this tunnel is 54.37 [s].
.PP
It should be noted that the total temperature of starts at 300 K, which is less than the minimum total temperature of 330 K, and decreases from there.
This indicates that these calculated results are inaccurate since they do not account for liquefaction.
.PSPIC ./ps/mocNozzle_time.eps
.label "Figure 6.  Possible Operation Time per Initial Reservoir Pressure"
.NH 3
Unit Reynolds Number
.PP
The unit Reynolds number in the free-stream of the wind tunnel is shown in the plot in Figure 7 against the total pressure.
This plot is linear since the driving equation described above (17) is linear with total pressure.
The unit Reynolds number in the free-stream of this tunnel ranges from 9029941 to 85334671 [1/m].
.PSPIC ./ps/mocNozzle_Re.eps
.label "Figure 7.  Unit Reynolds Number vs Total Pressure"
.NH 3
Tunnel Geometry
.PP
This nozzle is 1.0494 [m] long.
For the Mach number in the settling chamber to be \(<= 0.05, it would have to have a cross-sectional area of 0.08028 $[ roman m sup 2 ]$ or 80280 $[ roman {mm} sup 2 ]$.
.NH 3
Model Forces
.PP
The forces on the wedge would be 179.5354 [N] to the right and 810.6305 [N] downward if the tunnel was operating at full total pressure.

.NH
Usage
.PP
The mocNozzle.py script should be called from the command line to give it the desired arguments.
The --help flag gives a list of the flags with their descriptions.
By default, the script does all the calculations, defaulting to plotting the nozzle with 256 starting points.
Flags can be used to change the number of initial points, Mach number, throat radius, or turn off calculation off specific calculations such as the convergence test.
The README.md document contains more information on the script's usage.

.NH
Summary and Conclusions
.PP
This script calculates the contour of a wind tunnel nozzle based on the given parameters and generates the contour plot describing the flow inside it.
The script also calculates flow characteristics useful for understanding its behavior in a wind tunnel.
It generates plots and prints out data describing the tunnel operation, geometry, and limits.

.SH
Appendix
.PP
There are many functions used to organize the calculations in this script.
There are separate functions for generating the specific plots and wind tunnel characteristics.
The gasdynamics.py file contains a collection of function for the basic relationships commonly used.
.PP
calcInit() calculates the positions of the initial points and the flow state values at those points.
.PP
calcFirstPoint() decides the position of the first point by finding where its curve ends on the top wall and changing its starting angle to make that halfway between the points on either side.
It in turn calls calcInit and calculates the first 2 curves.
.PP
calcCurve() calculates the states at the points along characteristic curve starting from the top wall, bouncing off the centerline, and ending with a point on the top wall.
.PP
calculateNozzle() calls calcFirstPoint and then loops through each curve calling calcCurve untill the nozzle is nozzle is finished.
