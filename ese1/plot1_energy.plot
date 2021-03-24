#!/usr/bin/env gnuplot

# set the output as pdf
set term pdfcairo
# Input file contains comma-separated values fields
set output "data/energy1.pdf"
set title "even eigenvalues"
set xlabel "E"
set ylabel "y_xmax"
# set log y
# set xrange [-1:1]
set yrange [-1000:1000]
set grid
plot "data/energy1.dat" using 1:2 with lines title "even" \
    ,"data/energy1.dat" using 1:3 with lines title "odd"
# write to file
set output
