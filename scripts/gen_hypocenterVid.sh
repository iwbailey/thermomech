#!/bin/bash
#
# use gnuplot to plot the distribution of a set of numbers on the fault plane
#

pname=test_bz1996analog

ifile=./${pname}.out

# fault dimensions
nl=128;
xlen=70;
nd=32;
zlen=17.5;

dtp=4;

######################################################################
function plothypo()
{
    # plot one of the ifile columns before and after the earthquake

    ofile=$1 # output filename
    t1=$2;
    t2=$3;

    echo $t1 $t2

    cat $ifile | \
	gawk '{
if($2>=T1&&$2<T2)print($4,$6,S*$8-3.5,(65536 + 256 + 1)*int(255*(T2-$2)/(T2-T1)))}' \
	S=1 T1=$t1 T2=$t2 | sort -n -k 4 -r > tmp.dat

    #cat tmp.dat

    if [ -s tmp.dat ]; then
	gnuplot <<EOF
set terminal pngcairo
set output '$ofile'
set nokey
set xrange [0:${xlen}]
set yrange [${zlen}:0]
set title "t=$t2"
set size ratio -1
rgb(r,g,b) = 65536 * int(r) + 256 * int(g) + int(b)
plot "tmp.dat" using 1:2:3:4 \
     with circles lc rgb var fs transparent solid 0.3 noborder
# lc rgb "blue" fs transparent solid  noborder
EOF
    fi
}



######################################################################

minmax $ifile

t0=0.0; # start time
dt=.1 # plot interval
np=700; # number of plots

rm hypopngfiles/*

for ip in `seq 1 $np`; do
#for ip in 100; do

    t2=`echo "$t0 + $ip * $dt" | bc -l`
    t1=`echo "$t2 - $dtp" | bc -l`

    idx=`printf "%06i" $ip`

    # plot distribution
    plothypo ./hypopngfiles/${pname}.$idx.hypo.png $t1 $t2

done

# convert to movie file
vidfile=hypovid.avi
mencoder -mc 0 -noskip -skiplimit 0 -ovc lavc -lavcopts \
     vcodec=msmpeg4v2:vhq \
    "mf://hypopngfiles/${pname}.*.hypo.png" -mf \
    type=png:fps=15 -o $vidfile

mplayer $vidfile
