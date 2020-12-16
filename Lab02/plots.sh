set terminal post enhanced colour solid font 25  
set key vertical opaque outside top center
set grid
set autoscale 
set output "plots.eps"
set xlabel "t"
set ylabel "n"

### 1 ###
set title "Iteracja Picarda"
plot 'out1.dat'  u 1:2 w lines lw 5 t "u(t)" , \
	 '' 		 u 1:3 w lines lw 5 t "z(t)"
	 
### 2 ###
set title "Iteracja Newtona"
plot 'out2.dat'  u 1:2 w lines lw 5 t "u(t)" , \
	 '' 		 u 1:3 w lines lw 5 t "z(t)"

### 3 ###
set title "Niejawna metoda RK2"
plot 'out3.dat'  u 1:2 w lines lw 5 t "u(t)" , \
	 '' 		 u 1:3 w lines lw 5 t "z(t)"



