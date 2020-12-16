set terminal post enhanced colour solid font 20
set output "plots.eps"
set view map;
set size square;
set colorbox border 2604
set xlabel 'x'
set ylabel 'y'
set autoscale;
set palette defined (0 'blue', 20 'red', 40 'yellow')

### T ###
set title 'T(x,y) it = 100'
splot 'out.dat' i 0 with pm3d notitle
set title '{/Symbol \321}^{2}T(x,y) it = 100'
splot 'out.dat' i 1 with pm3d notitle
set title 'T(x,y) it = 200'
splot 'out.dat' i 2 with pm3d notitle
set title '{/Symbol \321}^{2}T(x,y) it = 200'
splot 'out.dat' i 3 with pm3d notitle
set title 'T(x,y) it = 500'
splot 'out.dat' i 4 with pm3d notitle
set title '{/Symbol \321}^{2}T(x,y) it = 500'
splot 'out.dat' i 5 with pm3d notitle
set title 'T(x,y) it = 1000'
splot 'out.dat' i 6 with pm3d notitle
set title '{/Symbol \321}^{2}T(x,y) it = 1000'
splot 'out.dat' i 7 with pm3d notitle
set title 'T(x,y) it = 2000'
splot 'out.dat' i 8 with pm3d notitle
set title '{/Symbol \321}^{2}T(x,y) it = 2000'
splot 'out.dat' i 9 with pm3d notitle
