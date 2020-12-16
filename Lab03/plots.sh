set terminal post enhanced colour solid font 25  
set key vertical opaque outside top center
set grid
set autoscale 
set output "plots.eps"

### 1 ###
set xlabel "t"
set ylabel "v(t)"
set title "Metoda trapezow - v(t)"
plot 'out1.dat'  i 1 u 1:4 w lines lw 2 t "TOL=1e-5" , \
	 '' 		 i 0 u 1:4 w lines lw 2 t "TOL=1e-2"

set ylabel "x(t)"
set title "Metoda trapezow - x(t)"
plot 'out1.dat'  i 1 u 1:3 w lines lw 2 t "TOL=1e-5" , \
	 '' 		 i 0 u 1:3 w lines lw 2 t "TOL=1e-2"

set ylabel "deltaT(t)"
set title "Metoda trapezow - deltaT(t)"
plot 'out1.dat'  i 1 u 1:2 w lines lw 2 t "TOL=1e-5" , \
	 '' 		 i 0 u 1:2 w linespoints pt 7 ps 1 t "TOL=1e-2"

set xlabel "x"
set ylabel "v(x)"
set title "Metoda trapezow - v(x)"
plot 'out1.dat'  i 1 u 3:4 w lines lw 2 t "TOL=1e-5" , \
	 '' 		 i 0 u 3:4 w lines lw 2 t "TOL=1e-2"


### 2 ###
set xlabel "t"
set ylabel "v(t)"
set title "Metoda RK2 - v(t)"
plot 'out2.dat'  i 1 u 1:4 w lines lw 2 t "TOL=1e-5" , \
	 '' 		 i 0 u 1:4 w lines lw 2 t "TOL=1e-2"

set ylabel "x(t)"
set title "Metoda RK2 - x(t)"
plot 'out2.dat'  i 1 u 1:3 w lines lw 2 t "TOL=1e-5" , \
	 '' 		 i 0 u 1:3 w lines lw 2 t "TOL=1e-2"

set ylabel "deltaT(t)"
set title "Metoda RK2 - deltaT(t)"
plot 'out2.dat'  i 1 u 1:2 w lines lw 2 t "TOL=1e-5" , \
	 '' 		 i 0 u 1:2 w linespoints pt 7 ps 1 t "TOL=1e-2"

set xlabel "x"
set ylabel "v(x)"
set title "Metoda RK2 - v(x)"
plot 'out2.dat'  i 1 u 3:4 w lines lw 2 t "TOL=1e-5" , \
	 '' 		 i 0 u 3:4 w lines lw 2 t "TOL=1e-2"
	 




