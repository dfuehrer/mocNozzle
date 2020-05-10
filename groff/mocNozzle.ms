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


.NH
Methods
.NH 2
Method of Characteristics Nozzle
.PP
The mocNozzle.py script has a class with methods to calculate the necessary wind tunnel data.
It first calculates the flow properties in the nozzle, using them to find its contour.
The nozzle starts at the throat as a section of a circle.
Points on the circle are evenly spaced out up to $delta sub {roman max}$, with the first point offset from the center of the throat.  The script starts by choosing the first point such that where its path ends on the wall is halfway between the points on the wall on either side.
For each point, the script calculates the Mach number, direction of the flow, Prandtl-Meyer angle, Mach angle, left and right running characteristic values, and $x$ and $y$ coordinates using the equations below.
At each point, the left and right running characteristics are the same as the preceding points.
Then the flow direction and Prandtl-Meyer angles are calculated from the equations below.
.EQ
nu = {R + L} over 2
.EN
.EQ
delta = {R - L} over 2
.EN
From $nu$, the Mach number can be found from numerically solving the Prandtl-Meyer equation, and $mu$ can be found from the Mach number.
Lastly the coordinates of the new point are found from the intersection of the of the Mach lines at each point.
At the centerline, the flow direction is known to be horizontal, and the other values follow from there. 

.PP
At the top wall, the shape has to cancel out the incoming characteristic curves.
Since the curves bouncing off the wall will intersect with the final curve of the design Mach number with horizontal flow:
.EQ
R sub {roman wall} = L sub {roman final}
.EN
From this, the full state is known and all other values can be calculated.
.NH 2
Nozzle Contour
.NH 3
Convergence
.PP
The error in the calculation can be simply measured by comparing the end height to the theoretical area ratio for a nozzle of this design Mach number.
.EQ
A over {A sup \(**} =
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
The top plot showed the points in their $x$ and $y$ coordinates, bounded by the black line of the nozzle wall.
During the nozzle calculation, an ID was stored for the left and right curves at each point.
Then each point was connected in a line based on their left and right IDs.
The bottom plot is a contour map of the flow within the nozzle.
To plot a contour of the whole nozzle rather than solely the calculated points, two points were added at the top and bottom of the throat and one point at the end was added on the centerline.
The colors in both plots vary based on the Mach number at the point.
.NH 3
Centerline Distributions
.PP
After calculating the nozzle, the total pressure and total temperature ratios were calculated at each point.
Based on the IDs mentioned above, the pressure ratios, temperature ratios, and Mach numbers were selected to just the points on the centerline.
These values were then plotted against their $x$ position in one plot with separate y-axes to make the values of the ratios clearly visible compared to Mach number.
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
If the total pressure was lower, then the flow would be over-expanded, causing shocks to form at the exit, and the nozzle would not be started.
.NH 3
Mass Flow Rate
.PP
The mass flow rate is given by:
.EQ
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
.EQ
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
.EQ
T sub 0 =
T sub {0 0} left ( p sub 0 over p sub {0 0} right ) sup {{ gamma - 1 } over gamma }
.EN
So
.EQ
p sub 0 over T sub 0 =
{{p sub {0 0}} sup {{ gamma - 1 } smallover gamma } {p sub 0} sup { 1 smallover gamma }} over T sub {0 0}
.EN
.EQ
V over R left ( {d ( p sub 0 smallover T sub 0 )} over {dt} right ) =
V over R left ( {p sub {0 0}} sup {{ gamma - 1 } smallover gamma } over T sub {0 0} right ) ^ left ( {d ( {p sub 0} sup { 1 smallover gamma } )} over {dt} right ) =
V over R left ( {p sub {0 0}} sup {{ gamma - 1 } smallover gamma } over T sub {0 0} right ) ^ left ( 1 smallover gamma {p sub 0} sup {{ 1 - gamma } smallover gamma } right ) ^ left ( {d p sub 0} over {dt} right )
.EN
.EQ
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
.EQ
c =
{ gamma sqrt { gamma R T sub {0 0} }}
over V
{p sub {0 0}} sup {{ 1 - gamma } smallover gamma }
A left ( A sup \(** over A right ) ^
left ( 2 over { gamma + 1 } right ) sup 
{{ gamma + 1 } over { 2 ( gamma - 1 ) }}
.EN
.EQ
{d p sub 0} over {dt} =
-c ^
left ( {p sub 0} sup {{ 3 gamma - 1 } over gamma } right ) 
.EN
Separate and integrate.
Find $C sub 0$ from $p sub 0 ( 0 ) = p sub {0 0}$.
.EQ
left ( { 2 gamma } over { 1 - gamma } right ) 
{p sub 0} sup {{1 - gamma } over { 2 gamma }} =
-ct + C sub 0 = 
-ct + 
left ( { 2 gamma } over { 1 - gamma } right ) ^
{p sub {0 0}} sup {{1 - gamma } over { 2 gamma }}
.EN
Solving for $t$ at $p sub {0 roman min}$
.EQ
t =
left ( { 2 gamma } over { c ( gamma - 1 )} right ) ^
left ( {p sub {0 roman min}} sup {{1 - gamma } over { 2 gamma }} - {p sub {0 0}} sup {{1 - gamma } over { 2 gamma }} right )
.EN
Using this expression, the operation time $t$ was plotted against the operation total pressure $p sub 0$.
.NH 3
Unit Reynolds Number
.PP
The unit Reynolds number is given by the equation below, using Sutherland's law to approximate viscosity with respect to the temperature.
.EQ
Re prime =
{ rho U } over mu =
{ rho U } 
over {
mu sub {roman ref} left ( T over T sub {roman ref} right ) sup {3/2} ^
left ( { T sub {roman ref} + S } over {  T + S } right )
} =
{ ( p smallover p sub 0 ) p sub 0 U } 
over {
mu sub {roman ref} 
R ( T smallover T sub 0 ) T sub {0 roman min}
left ( T over T sub {roman ref} right ) sup {3/2} ^
left ( { T sub {roman ref} + S }
over
{  ( T smallover T sub 0 ) T sub {0 roman min} + S } right )
}
.EN
where $mu sub {roman ref} = 1.716 E -5 ~[ {roman Pa} fat cdot {roman s} ]$, $T sub {roman ref} = 273.15 roman K$, and $S = 110.4 roman K$.
A plot of unit Reynolds number versus operating total pressure was made using the above expression.
.NH 3
Tunnel Geometry
.PP
The length of the nozzle was found by multiplying the maximum $x$ value by the throat height.
The throat height was found by dividing the exit height by the area ratio, since the tunnel has a constant width.
.PP
The area ratio needed for a $M <= 0.05$ flow was found using the area ratio equation with that Mach number
Multiplying the area ratio by the width times the throat height gives the area needed.
.NH 3
Model Forces
.PP
The forces on the wedge were found first by calculating the areas of the sides using the given angle and the width.
Then normal vectors to the sides were defined based on the sines and cosines of the angle.
After that, the angle was used to find the angle of the shock wave on the top and bottom.
From there, the normal components of the flow to the oblique shocks were used to find the pressures on the top and bottom of the wedge with the normal shock relations.
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
.NH 3
Nozzle Contour
.PP
Below is the plot of the nozzle contour.
As discussed above, the final height is very close to the theoretical 12.5.
The nozzle is quite long, about 75 times the length of the throat height.
Because of this, the semi-circular portion of the nozzle at the throat is indistinguishable.
The Mach contour quickly increases along the x-axis, with the Mach number reaching 3 around 10% of the total length and reaching 5 around 20% of the length.
.PSPIC ./ps/mocNozzle_nozzle_M5_num512_rot.eps
.NH 3
Centerline Distributions
.PP
The centerline distributions of Mach number, total temperature and total pressure are shown in the plot below.
This plot shows how quickly the Mach number rises from 1 to 5, where it stays steady to the end of the length.
Since the total temperature and total pressure ratios are dependent on the Mach number, they in turn fall quickly from their initial values of .8333 and .5283 to .1667 and .00189, respectively.
As indicated above, $P sub {cr3} = 0.00189$.
From this chart, it is clear that the total pressure ratio drops more quickly than the total temperature ratio.
.PSPIC ./ps/mocNozzle_MpT.eps
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
The minimum total pressure found was 330 K; total temperatures less than this would result in liquefaction inside the nozzle.
.NH 3
Minimum Total Pressure
.PP
The minimum total pressure for the wind tunnel is 52.9090 [psia].
.NH 3
Mass Flow Rate
.PP
The mass flow rate for a nozzle 500 [mm] wide with an exit height of 350 [mm] operating at the minimum total temperature is plotted below.
The relationship is linear with the operating total pressure, as mandated by the driving expression, ranging from 5.62 to almost 53.12 $[ roman { kg smallover s } ]$.
.PSPIC ./ps/mocNozzle_mdot.eps
.NH 3
Operation Time
.PP
The possible operation time is plotted against the initial total pressure in the chart below.
The relationship is nonlinear, with increases in initial total pressure bringing diminishing returns in possible operation time.
The operation time given by this graph is less than it would be if the total temperature was assumed to remain constant.
The maximum operation time for this tunnel is 54.37 [s].
.PP
It should be noted that the total temperature of 300 K starts at, which is less than the minimum total temperature of 330 K, and decreases from there.
This indicates that these calculated results are inaccurate since they do not account for liquefaction.
.PSPIC ./ps/mocNozzle_time.eps
.NH 3
Unit Reynolds Number
.PP
The unit Reynolds number in the free-stream of the wind tunnel is shown in the plot below against the total pressure.
This plot is linear since the driving equation described above is linear with total pressure.
The unit Reynolds number in the free-stream of this tunnel ranges from 9029941 to 85334671 [1/m].
.PSPIC ./ps/mocNozzle_Re.eps
.NH 3
Tunnel Geometry
.PP
This nozzle is 1.0494 [m] long.
For the Mach number in the settling chamber to be \(<= 0.05, it would have to have a cross-sectional area of .08028 $[ roman m sup 2 ]$ or 80280 $[ roman {mm} sup 2 ]$.
.NH 3
Model Forces
.PP
The forces on the wedge would be 179.5354 [N] to the right and 810.6305 [N] downward if the tunnel was operating at full total pressure.

.NH
Summary and Conclusions
.PP


.SH
Appendix

