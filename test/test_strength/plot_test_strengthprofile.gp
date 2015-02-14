ifile = "test_strengthprofile.out"
ofile = "test_strengthprofile.png"

set term png

set output ofile

set grid back ls 12
set xlabel "Depth [km]"
set ylabel "Strength [MPa]"
set yrange [0:300]
plot ifile using 1:2 with lines lw 2 title "Static Strength", \
     ifile using 1:3 with linespoints lw 2 title "Creep Strength", \
     ifile using 1:4 with lines lw 2 title "Strength", \
     ifile using 1:5 with lines lw 2 title "BZ1996 Creep Strength", \

print sprintf("Written to %s",ofile)