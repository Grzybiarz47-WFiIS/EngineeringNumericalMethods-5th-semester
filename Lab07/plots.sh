set terminal post enhanced colour solid font 20
set output "plots.eps"
set view map
set xlabel 'x'
set ylabel 'y'
set palette defined (-1 'blue', 0 'red', 1 'yellow')
set tics font "Helvetica,15"
unset surface 

### Psi Q = -1000 ###
set contour
set cbrange [-55:-50]
set cntrparam levels increment -55, 0.25, -50
set title '{/Symbol y}(x,y) Q = -1000'
splot 'out-1000.dat' i 0 u 1:2:3 w l lt -1 lw 2 palette notitle

### Dzeta Q = -1000 ###
set cbrange [-200:300]
set cntrparam levels increment -200, 20, 300
set title '{/Symbol z}(x,y) Q = -1000'
splot 'out-1000.dat' i 1 u 1:2:3 w l lt -1 lw 2 palette notitle
unset contour

### u v Q = -1000 ###
set autoscale
set surface
set title 'u(x,y) Q = -1000'
splot 'out-1000.dat' i 2 with pm3d notitle
set title 'v(x,y) Q = -1000'
splot 'out-1000.dat' i 3 with pm3d notitle
unset surface

### Psi Q = -4000 ###
set contour
set cbrange [-220:-200]
set cntrparam levels increment -220, 1, -200
set title '{/Symbol y}(x,y) Q = -4000'
splot 'out-4000.dat' i 0 u 1:2:3 w l lt -1 lw 2 palette notitle

### Dzeta Q = -4000 ###
set cbrange [-400:1000]
set cntrparam levels increment -400, 50, 1000
set title '{/Symbol z}(x,y) Q = -4000'
splot 'out-4000.dat' i 1 u 1:2:3 w l lt -1 lw 2 palette notitle
unset contour

### u v Q = -4000 ###
set autoscale
set surface
set title 'u(x,y) Q = -4000'
splot 'out-4000.dat' i 2 with pm3d notitle
set title 'v(x,y) Q = -4000'
splot 'out-4000.dat' i 3 with pm3d notitle
unset surface

### Psi Q = 4000 ###
set contour
set cbrange [200:220]
set cntrparam levels increment 200, 1, 220
set title '{/Symbol y}(x,y) Q = 4000'
splot 'out4000.dat' i 0 u 1:2:3 w l lt -1 lw 2 palette notitle

### Dzeta Q = 4000 ###
set cbrange [-1000:400]
set cntrparam levels increment -1000, 50, 400
set title '{/Symbol z}(x,y) Q = 4000'
splot 'out4000.dat' i 1 u 1:2:3 w l lt -1 lw 2 palette notitle
unset contour

### u v Q = 4000 ###
set autoscale
set surface
set title 'u(x,y) Q = 4000'
splot 'out4000.dat' i 2 with pm3d notitle
set title 'v(x,y) Q = 4000'
splot 'out4000.dat' i 3 with pm3d notitle
unset surface