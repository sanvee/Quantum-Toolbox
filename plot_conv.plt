FILENAME ='converg2.txt'

set datafile separator "\t"
#set decimalsign locale "de_DE.UTF-8"
set decimalsign '.'
set autoscale x
set autoscale y
set autoscale z
set view equal xyz
set ticslevel 0
set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 1.5 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 1.5 # --- green

set mapping spherical
set angles radians

set parametric
#set hidden3d
set view 60
set isosamples 5,6
set xrange[-1 : 1]
set yrange[-1 : 1]
set zrange[-1 : 1]
splot [-pi:pi][-pi/2:pi/2] cos(u)*cos(v), sin(u)*cos(v), sin(v),\
FILENAME using 2:3:1:4 w p ls 1 lc palette