# set terminal and output
set terminal png size 900,400
set output "output.png"

# set global properties
set style line 1 lw 2 lc rgb "red"
set style line 2 lw 2 lc rgb "blue"
set style line 3 lw 2 lc rgb "green"
set style line 4 lw 2 lc rgb "purple"
set style line 5 lw 2 lc rgb "orange"

# Multiplot to plot different figures on one canvas
set multiplot layout 1,1
unset key

set xtics 100 




set ylabel "Velocity"
unset xlabel
plot    "d_Ln.dat" using ($1*10):7 with lines lw 2 lc rgb "red" notitle,\
        "ev_con-n.dat" using ($1*10):5 with lines lw 2 lc rgb "red" notitle, \
        "d_Rn.dat" using ($1*10):7 with lines lw 2 lc rgb "red" notitle, \
        "./with_cond/20000/d_Ln.dat" using ($1*10):7 with lines lw 2 lc rgb "blue" notitle, \
        "./with_cond/20000/ev_con-n.dat" using ($1*10):5 with lines lw 2 lc rgb "blue" notitle, \
        "./with_cond/20000/d_Rn.dat" using ($1*10):7 with lines lw 2 lc rgb "blue" notitle,\
        "./without_cond/20000/d_Ln.dat" using ($1*10):7 with lines lw 2 lc rgb "black" notitle, \
        "./without_cond/20000/ev_con-n.dat" using ($1*10):5 with lines lw 2 lc rgb "black" notitle, \
        "./without_cond/20000/d_Rn.dat" using ($1*10):7 with lines lw 2 lc rgb "black" notitle

unset multiplot