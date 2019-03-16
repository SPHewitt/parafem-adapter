#! /bin/bash
gnuplot -p << EOF
set grid
set title 'Displacement of the Flap Tip'
set xlabel 'Time [s]'
set ylabel 'Y-Displacement [m]'
set xrange [2:]
plot "precice-ParaFEM-watchpoint-point1.log" using 1:9 with lines title ""

EOF

