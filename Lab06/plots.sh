set terminal post enhanced colour solid font 20
set output "plots.eps"
set view map;
set size square;
set colorbox border 2604
set xlabel 'x'
set ylabel 'y'
set autoscale;
set palette defined (-1 'blue', 0 'white', 1 'red')

### V ###
set title 'V(x,y) n_{x} = n_{y} = 50'
splot 'out.dat' i 1 with pm3d notitle
set title 'V(x,y) n_{x} = n_{y} = 100'
splot 'out.dat' i 2 with pm3d notitle
set title 'V(x,y) n_{x} = n_{y} = 200'
splot 'out.dat' i 3 with pm3d notitle
set cbrange [-0.8:0.8]
set title 'V(x,y) n_{x} = n_{y} = 100 {/Symbol e}_{1} = 1 {/Symbol e}_{2} = 1'
splot 'out.dat' i 4 with pm3d notitle
set title 'V(x,y) n_{x} = n_{y} = 100 {/Symbol e}_{1} = 1 {/Symbol e}_{2} = 2'
splot 'out.dat' i 5 with pm3d notitle
set title 'V(x,y) n_{x} = n_{y} = 100 {/Symbol e}_{1} = 1 {/Symbol e}_{2} = 10'
splot 'out.dat' i 6 with pm3d notitle