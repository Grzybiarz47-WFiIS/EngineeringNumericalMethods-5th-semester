set terminal post enhanced colour solid font 25  
set key vertical opaque outside top center
set grid
set autoscale 
set output "plots.eps"

### 1 ###
set title "Metoda jawna Eulera"
set xlabel "t"
set ylabel "y(t)"
plot 'out1.dat'  i 0 u 1:2 w linespoints pt 7 ps 0.5 t "deltaT 0.01" , \
	 '' 		 i 1 u 1:2 w linespoints pt 7 ps 0.5 t "deltaT 0.1" , \
	 '' 		 i 2 u 1:2 w linespoints pt 7 ps 0.5 t "deltaT 1.0", \
	 			 exp(-x) w lines lw 2 t "e^-^x"

set title "Metoda jawna Eulera - bledy"
set ylabel "ERROR"		 
set logscale y
set format y "10^{%L}"
plot 'out1.dat' i 0 u 1:3 w lines t "deltaT 0.01", \
	 '' 		i 1 u 1:3 w linespoints pt 7 ps 0.5 t "deltaT 0.1", \
	 '' 		i 2 u 1:3 w linespoints pt 7 ps 0.5 t "deltaT 1.0"
unset logscale y
unset format y
	 
### 2 ###
set title "RK2"
set ylabel "y(t)"
plot 'out2.dat'  i 0 u 1:2 w linespoints pt 7 ps 0.5 t "deltaT 0.01" , \
	 '' 		 i 1 u 1:2 w linespoints pt 7 ps 0.5 t "deltaT 0.1" , \
	 '' 		 i 2 u 1:2 w linespoints pt 7 ps 0.5 t "deltaT 1.0", \
	 			 exp(-x) w lines t "e^-^x"

set title "RK2 - bledy"
set ylabel "ERROR"
set logscale y
set format y "10^{%L}"
plot 'out2.dat' i 0 u 1:3 w lines t "deltaT 0.01", \
	 '' 		i 1 u 1:3 w linespoints pt 7 ps 0.5 t "deltaT 0.1", \
	 '' 		i 2 u 1:3 w linespoints pt 7 ps 0.5 t "deltaT 1.0"
unset logscale y
unset format y

### 3 ###
set title "RK4"
set ylabel "y(t)"
plot 'out3.dat'  i 0 u 1:2 w linespoints pt 7 ps 0.5 t "deltaT 0.01" , \
	 '' 		 i 1 u 1:2 w linespoints pt 7 ps 0.5 t "deltaT 0.1" , \
	 '' 		 i 2 u 1:2 w linespoints pt 7 ps 0.5 t "deltaT 1.0", \
	 			 exp(-x) w lines t "e^-^x"

set title "RK4 - bledy"
set ylabel "ERROR" 
set logscale y
set format y "10^{%L}"
plot 'out3.dat' i 0 u 1:3 w lines t "deltaT 0.01", \
	 '' 		i 1 u 1:3 w linespoints pt 7 ps 0.5 t "deltaT 0.1", \
	 '' 		i 2 u 1:3 w linespoints pt 7 ps 0.5 t "deltaT 1.0"
unset logscale y
unset format y

### 4 ###
set title "RRZ2 - Q(t)"
set xlabel "t"
set ylabel "Q(t)"
plot 'out4.dat' i 0 u 1:2 w lines t "0.5w0", \
	 '' 		i 1 u 1:2 w lines t "0.8w0", \
	 '' 		i 2 u 1:2 w lines t "1.0w0", \
	 ''			i 3 u 1:2 w lines t "1.2w0"
	 
set title "RRZ2 - I(t)"
set ylabel "I(t)"
plot 'out4.dat' i 0 u 1:3 w lines t "0.5w0", \
	 '' 		i 1 u 1:3 w lines t "0.8w0", \
	 '' 		i 2 u 1:3 w lines t "1.0w0", \
	 ''			i 3 u 1:3 w lines t "1.2w0"
	 

