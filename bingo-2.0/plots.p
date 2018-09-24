set term postscript eps enhanced color 

se output "plots/phi.eps"
set ytics format "%g"
set size 1,0.8
set title 'Behaviour of {/Symbol f} as a function of N'
set xtics nomirror font "Helvetica,20" 
set ytics nomirror font "Helvetica,20"

set xlabel 'N'
set ylabel '{/Symbol f}[N]'
plot 'plots/PHI.dat' w l lw 2


se output "plots/dphi_dn.eps"
set title 'Behaviour of {/Symbol f_N} as a function of N'
set size 1,0.8
set xtics nomirror font "Helvetica,20" 

set ytics nomirror font "Helvetica,20"

set xlabel 'N'
set ylabel '{/Symbol f_N}[N]'

pl 'plots/DPHI_DN.dat' w l lc rgb 'blue' lw 2 


se output "plots/epsilon.eps"
set title 'Behaviour of {/Symbol e} as a function of N'
set size 1,0.8
set xtics nomirror font "Helvetica,20" 

set ytics nomirror font "Helvetica,20"

set xlabel 'N'
set ylabel '{/Symbol e}[N]'

pl 'plots/EPSILON.dat' w l lc rgb 'blue' lw 2 

se output "plots/powerspectra.eps"
set title 'Behaviour of P_S(k) as a function of k'
set size 1,0.8
set log
set xtics nomirror font "Helvetica,20" 

set ytics nomirror font "Helvetica,20"
set xrange [1e-5:1e-1]
set xlabel 'k'
set ylabel 'P_S[k]'

pl 'plots/PS.txt' w l lc rgb 'blue' lw 2 

se output "plots/f_nl.eps"
set title 'Behaviour of f_{NL}(k) as a function of k'
set size 1,0.8
unset log
set log x
set xtics nomirror font "Helvetica,20" 

set ytics nomirror font "Helvetica,20"

set xlabel 'k'
set ylabel 'f_{NL}[k]'

pl 'plots/F_NL-RE.txt' w l lc rgb 'blue' lw 2