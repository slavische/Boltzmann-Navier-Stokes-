# set terminal and output
set terminal png size 900,400
set output "output2.png"

# set global properties
set style line 1 lw 2 lc rgb "red"
set style line 2 lw 2 lc rgb "blue"
set style line 3 lw 2 lc rgb "green"
set style line 4 lw 2 lc rgb "purple"
set style line 5 lw 2 lc rgb "orange"

# Multiplot to plot different figures on one canvas
set multiplot layout 1,3
unset key

set xtics 100 


set ylabel "Supersaturation"
set xlabel "Coordinate"
plot    "ss_ratio-2000.dat" using ($1*10):2 with lines lw 2 notitle

set ylabel "Number of droplet"
plot    "mom-2000.dat" using ($1*10):2 with lines lw 2 notitle

set ylabel "Liquid mass fraction"
plot    "mom-2000.dat" using ($1*10):5 with lines lw 2 notitle,

unset multiplot