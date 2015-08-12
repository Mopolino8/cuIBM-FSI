reset;
set terminal pdf enhanced color font 'Palatino, 11' size 15cm, 15cm;

set title 'Velocity along the vertical centerline (Re=100)'
set xlabel 'y-coordinate'
set ylabel 'Centerline u-velocity'
set output '/home/chris/src/cuIBM/validation/lidDrivenCavity/Re100/uRe100.pdf'
plot [0:1] [-0.7:1.3] \
'/home/chris/src/cuIBM/validation-data/cavity-GGS82.txt' u 1:2 w p pt 13 ps 1 lc rgb '#4B5ED7' title 'Ghia et al, 1982', \
'/home/chris/src/cuIBM/validation/lidDrivenCavity/Re100/u.txt' u 1:2 w l lw 5 lc rgb '#228A4C' title 'cuIBM (Taira and Colonius, 2005)'

reset;
set terminal pdf enhanced color font 'Palatino, 11' size 15cm, 15cm;

set title 'Velocity along the horizontal centerline (Re=100)'
set xlabel 'x-coordinate'
set ylabel 'Centerline v-velocity'
set output '/home/chris/src/cuIBM/validation/lidDrivenCavity/Re100/vRe100.pdf'
plot [0:1] [-0.7:1.3] \
'/home/chris/src/cuIBM/validation-data/cavity-GGS82.txt' u 6:7 w p pt 13 ps 1 lc rgb '#4B5ED7' title 'Ghia et al, 1982', \
'/home/chris/src/cuIBM/validation/lidDrivenCavity/Re100/v.txt' u 1:2 w l lw 5 lc rgb '#228A4C' title 'cuIBM (Taira and Colonius, 2005)'
