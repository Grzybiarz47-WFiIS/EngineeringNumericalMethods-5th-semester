set terminal post enhanced colour solid font 20
set output "plots.eps"
set autoscale
set view map;
set colorbox border 2604

### global ###
set title 'Relaksacja globalna V(x,y) wG = 0.6'
set xlabel 'x'
set ylabel 'y'
splot 'out.dat' i 0 with pm3d 
set title 'Relaksacja globalna V(x,y) wG = 1.0'
splot 'out.dat' i 1 with pm3d 

set tics font "Helvetica,10"
### err ###
set title 'Relaksacja globalna error wG = 0.6'
splot 'err.dat' i 0 with pm3d 
set title 'Relaksacja globalna error wG = 1.0'
splot 'err.dat' i 1 with pm3d 

unset view;

set tics font "Helvetica,20"
### S ###
set xrange [1:100000]
set yrange [0:5000]
set logscale x
set title 'Relaksacja globalna S(it)'
set xlabel 'it'
set ylabel 'S'
plot 'S.dat' i 0 u 1:2 w lines lw 2 t 'wG = 0.6', \
     ''      i 1 u 1:2 w lines lw 2 t 'wG = 1.0'

### S ###
set title 'Relaksacja lokalna S(it)'
plot 'S.dat' i 2 u 1:2 w lines lw 2 t 'wL = 1.0', \
     ''      i 3 u 1:2 w lines lw 2 t 'wL = 1.4', \
     ''      i 4 u 1:2 w lines lw 2 t 'wL = 1.8', \
     ''      i 5 u 1:2 w lines lw 2 t 'wL = 1.9'