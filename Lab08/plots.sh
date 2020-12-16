set terminal post enhanced colour solid font 20
set output "plots.eps"
set view map
set xlabel 'x'
set ylabel 'y'
set palette defined (-2 'black', -1 'blue', 0 'red', 1 'yellow', 2 'white')
set tics font "Helvetica,15"
set size ratio 0.45
set autoscale

### vx vy ###
set title 'v_{x}(x,y)'
splot 'out1.dat' i 0 with pm3d notitle
set title 'v_{y}(x,y)'
splot 'out1.dat' i 1 with pm3d notitle

### u D = 0.0 ###
set cbrange[0:20]
set title 'u(x,y,t) D = 0.0 k = 1'
splot 'out2.dat' i 3 with pm3d notitle
set title 'u(x,y,t) D = 0.0 k = 2'
splot 'out2.dat' i 7 with pm3d notitle
set title 'u(x,y,t) D = 0.0 k = 3'
splot 'out2.dat' i 11 with pm3d notitle
set title 'u(x,y,t) D = 0.0 k = 4'
splot 'out2.dat' i 15 with pm3d notitle
set title 'u(x,y,t) D = 0.0 k = 5'
splot 'out2.dat' i 19 with pm3d notitle

### u D = 0.1 ###
set autoscale
set title 'u(x,y,t) D = 0.1 k = 1'
splot 'out2.dat' i 23 with pm3d notitle
set title 'u(x,y,t) D = 0.1 k = 2'
splot 'out2.dat' i 27 with pm3d notitle
set title 'u(x,y,t) D = 0.1 k = 3'
splot 'out2.dat' i 31 with pm3d notitle
set title 'u(x,y,t) D = 0.1 k = 4'
splot 'out2.dat' i 35 with pm3d notitle
set title 'u(x,y,t) D = 0.1 k = 5'
splot 'out2.dat' i 39 with pm3d notitle
unset view

### c x ###
set grid
set xlabel 't_{n}'
set ylabel 'c'
set title "c(t_{n})"
plot 'out3.dat'  i 0 u 1:2 w lines lw 3 t "D = 0.0" , \
	 '' 		 i 1 u 1:2 w lines lw 3 t "D = 0.1"

set ylabel 'x_{sr}'
set title "x_{sr}(t_{n})"
plot 'out3.dat'  i 0 u 1:3 w lines lw 3 t "D = 0.0" , \
	 '' 		 i 1 u 1:3 w lines lw 3 t "D = 0.1"


##### gif1 #####
reset
set term gif size 1000, 450 animate delay 20
set view map 

### u D = 0.0 ###
n1 = 19
set cbrange[0:20]
set title 'u(x,y,t) D = 0.0'
set output "anim1.gif"
do for [j=0:n1] {
  splot 'out2.dat' i j with pm3d notitle
}

##### gif2 #####
reset
set term gif size 1000, 450 animate delay 20
set view map 

### u D = 0.1 ###
n2 = 39
set autoscale
set title 'u(x,y,t) D = 0.1'
set output "anim2.gif"
do for [j=n1+1:n2] {
  splot 'out2.dat' i j with pm3d notitle
}
