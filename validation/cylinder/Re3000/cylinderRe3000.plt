reset;
set terminal pdf enhanced color font 'Palatino, 11' size 20cm, 15cm;

set title 'Flow over an impulsively started cylinder at Reynolds number 3000'
set xlabel 'Non-dimensional time'
set ylabel 'Drag Coefficient'
set output '/home/chris/src/cuIBM/validation/cylinder/Re3000/cylRe3000Drag.pdf'
plot [0:3] [0:2] \
'/home/chris/src/cuIBM/validation-data/cylinderRe3000-KL95.txt' u (0.5*$1):2 w p pt 13 ps 1 lc rgb '#4B5ED7' title 'Koumoutsakos and Leonard, 1995', \
'/home/chris/src/cuIBM/validation/cylinder/Re3000/forces' u 1:(2*$2) w l lw 5 lc rgb '#228A4C' title 'cuIBM (Taira and Colonius, 2005)'
