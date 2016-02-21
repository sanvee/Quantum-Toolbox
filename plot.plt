reset
unset multiplot
set terminal qt 0
#set terminal svg size 2024,2024 enhanced font 'Verdana,10'
#set output 'uebersicht.svg'
#set size square 1,1
set termoption enhanced

set xtics font 'Verdana,6'
set ytics font 'Verdana,6'
#set termoption dash

FILENAME ='test_tab.txt'
FILENAME2 ='converg.txt'

FILENAME_D ='suspekt.txt'

set datafile separator "\t"
#set decimalsign locale "de_DE.UTF-8"
set decimalsign '.'
set autoscale x
set autoscale y
set autoscale z
set view equal xyz
set ticslevel 0
#set xrange [0:1]
#set yrange [0:1]
#set zrange [0:1]
#set multiplot layout 4, 4

set style line 1 lc rgb '#8b1a0e' pt 1 ps 1 lt 1 lw 1.5 # --- red
set style line 2 lc rgb '#5e9c36' pt 6 ps 1 lt 1 lw 1.5 # --- green

#set style line 11 lc rgb '#808080' lt 1
#set border 3 back ls 11
#set tics nomirror
#set style line 12 lc rgb '#808080' lt 0 lw 1
#set grid back ls 12
 
#f(x)=1.0/3.0 
#stats FILENAME using($1)
#bin(x,width)=width*(floor(x/width)+0.5)
#bw = (STATS_max-STATS_min)/50
#set boxwidth bw*1.00

plot FILENAME2 using 3:4

#plot     FILENAME using (bin($1,bw)):(1.0)  smooth frequency with boxes ls 1,\
#         FILENAME using (bin($2,bw)):(1.0)  smooth frequency with boxes ls 2#,\
#        FILENAME using (bin($3,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME using (bin($4,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME using (bin($5,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME using (bin($6,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME using (bin($7,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME using (bin($8,bw)):(1.0)  smooth frequency with boxes
#     
#plot     FILENAME_D using (bin($1,bw)):(1.0)  smooth frequency with boxes ls 2
#plot     FILENAME_D using (bin($2,bw)):(1.0)  smooth frequency with boxes ls 1
#plot     FILENAME_D using (bin($3,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME_D using (bin($4,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME_D using (bin($5,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME_D using (bin($6,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME_D using (bin($7,bw)):(1.0)  smooth frequency with boxes
#plot     FILENAME_D using (bin($8,bw)):(1.0)  smooth frequency with boxes


#plot FILENAME using 2:1#, f(x)

#stats FILENAME using($1)
#bin(x,width)=width*(floor(x/width)+0.5)
#bw = (STATS_max-STATS_min)/50
#set boxwidth bw*1.00
#plot FILENAME using (bin($1,bw)):(1.0)  smooth frequency with boxes
#plot FILENAME using (bin($2,bw)):(1.0)  smooth frequency with boxes

unset multiplot