set term png
set output "test_strength_evol.png"

set grid back ls 12

set key box

set palette rgbformulae 33,13,10
set cbrange [0:100]
set cblabel "Time [kyr]"

set xlabel "Depth [km]"
set xrange [0:17.5]
set mxtics 4

set ylabel "Temperature [K]"
# set mytics 4

plot for [i=2:10000:100] './test_strength_evol.out' u 1:i w lines notitle \
     palette cb (i-1)*0.01

# set ylabel "Strength [MPa]"
# plot for [i=2:10000:100] './test_strength_evol.out' u 1:i w lines notitle \
#      palette cb (i-1)*0.01