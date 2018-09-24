se te po enhanced eps color  
se output "plots/F_nl_2d.eps"
set palette model RGB defined (0 "dark-blue", 0 "blue", 1 "cyan", 2 "yellow", 3 "red", 3 "dark-red")
set xlabel "k_2/k_1" font "Times-Roman,24"
set ylabel "k_3/k_1" font "Times-Roman,24"
set xtics 0,0.2,1 add (0.5) font "Times-Roman,20"
set ytics 0.5,0.1,1 font "Times-Roman,20"
set xrange [0:1]
set yrange [0.5:1]
f(x)=x
g(x)=-x+1
unset key

pl 'plots/F_nl_2d.txt' pt 11 ps 1 palette,f(x) w l lw 3 lc rgb "black" lt 1,g(x) w l lw 3 lc rgb "black" lt 1

#filename(n)=sprintf("plots/F_nl_%d.txt", n)

#pl for [i=1:40] filename(i) pt 11 ps 1 palette,f(x) w l lw 3 lc rgb "black" lt 1,g(x) w l lw 3 lc rgb "black" lt 1
