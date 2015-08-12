reset;
set terminal pdf enhanced color font 'Palatino, 11' size 20cm, 15cm;

set title 'Flow over a NACA 0008 airfoil at angle of attach 4^o'
set xlabel 'Non-dimensional time'
set ylabel 'Force Coefficient'
set output '/home/chris/src/cuIBM/validation/NACA0008/AoA04_0.004/AoA04_0.004.pdf'
plot [0:12] [0:0.5] \
'/home/chris/src/cuIBM/validation/NACA0008/AoA04_0.004/forces' u 1:(2*$3) w l lw 5 lc rgb '#4B5ED7' title 'Lift coefficient', \
'/home/chris/src/cuIBM/validation/NACA0008/AoA04_0.004/forces' u 1:(2*$2) w l lw 5 lc rgb '#228A4C' title 'Drag coefficient'
