set term postscript eps monochrome "Times,18"
set style data linespoints
set key right top

set output "l2error.eps"
set xlabel "degrees of freedom"
set logscale x
set logscale y
set grid
set ylabel "L_2 Error" 
plot \
"uniform_P1.dat" u 2:3 w l lw 4 title "uniform P1", \
"uniform_P2.dat" u 2:3 w l lw 4 title "uniform P2", \
"adaptive_P1.dat" u 2:3 w l lw 4 title "adaptive P1", \
"adaptive_P2.dat" u 2:3 w l lw 4 title "adaptive P2", \
1.0/x lw 4 title "1/N"
! epstopdf l2error.eps


