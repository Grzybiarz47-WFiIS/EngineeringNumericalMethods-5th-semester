set terminal post enhanced colour solid font 15
set output "plots.eps"
set xlabel 't'
set autoscale
set view map;
set colorbox border 2604
set palette defined (-1 'blue', 0 'white', 1 'red')

### E ###
set title 'E(t)'
set ylabel 'E(t)'
set grid
plot 'Eout.dat'  i 0 u 1:2 w lines lw 3 t "{/Symbol \a}=0.0 {/Symbol \b}=0.0", \
	 '' 		 i 1 u 1:2 w lines lw 3 t "{/Symbol \a}=0.0 {/Symbol \b}=0.1", \
     ''          i 2 u 1:2 w lines lw 3 t "{/Symbol \a}=0.0 {/Symbol \b}=1.0"
plot 'Eout.dat'  i 3 u 1:2 w lines lw 3 t "{/Symbol \a}=1.0 {/Symbol \b}=1.0"

### u(x,t) ###
set yrange [0:15]
set ylabel 'x'
set title 'u(x,t) {/Symbol a}=0.0 {/Symbol b}=0.0'
splot 'out.dat' i 0 with pm3d notitle
set title 'u(x,t) {/Symbol a}=0.0 {/Symbol b}=0.1'
splot 'out.dat' i 1 with pm3d notitle
set title 'u(x,t) {/Symbol a}=0.0 {/Symbol b}=1.0'
splot 'out.dat' i 2 with pm3d notitle
set title 'u(x,t) {/Symbol a}=1.0 {/Symbol b}=1.0'
splot 'out.dat' i 3 with pm3d notitle

