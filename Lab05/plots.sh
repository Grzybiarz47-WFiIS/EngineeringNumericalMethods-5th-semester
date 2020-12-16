set terminal post enhanced colour solid font 20
set output "plots.eps"
set view map;
set size square;
set colorbox border 2604
set xlabel 'x'
set ylabel 'y'
set xrange [0:26]
set yrange [0:26]
set palette defined (-1 'blue', 0 'yellow', 1 'red')

### V multigrid ###
set title 'Relaksacja wielosiatkowa V(x,y) k=16'
splot 'out.dat' i 0 with pm3d notitle
set title 'Relaksacja wielosiatkowa V(x,y) k=8'
splot 'out.dat' i 1 with pm3d notitle
set title 'Relaksacja wielosiatkowa V(x,y) k=4'
splot 'out.dat' i 2 with pm3d notitle
set title 'Relaksacja wielosiatkowa V(x,y) k=2'
splot 'out.dat' i 3 with pm3d notitle
set title 'Relaksacja wielosiatkowa V(x,y) k=1'
splot 'out.dat' i 4 with pm3d notitle

unset map;
unset size;
set autoscale;

### S stop ###
set title 'Relaksacja wielosiatkowa S(it)'
plot 'S.dat' i 0 u 1:2 w lines lw 5 t 'k=16', \
     ''      i 1 u 1:2 w lines lw 5 t 'k=8', \
     ''      i 2 u 1:2 w lines lw 5 t 'k=4', \
     ''      i 3 u 1:2 w lines lw 5 t 'k=2', \
     ''      i 4 u 1:2 w lines lw 5 t 'k=1'
