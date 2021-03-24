#!/usr/bin/env gnuplot

# set the output as pdf
set term pdfcairo
# Input file contains comma-separated values fields
set output "data/numerov1.pdf"
set title "Title"
set xlabel "x"
set ylabel "\Psi"
# set log y
# set xrange [-1:1]
set yrange [-2:2]
set grid
plot "data/numerov1.dat" using 1:2 with lines title "True" \
   ,"data/numerov1.dat" u 1:3 w lines title "Numerov"
# write to file
set output
