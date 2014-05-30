#l = 1 + Number of time steps(Lax-Friedrichs)
l=`awk -F'\t' 'NR == 1 { print NF; exit }' Density_LF.dat`-1
#m = 1 + Number of time steps(Roe)
m=`awk -F'\t' 'NR == 1 { print NF; exit }' Density_ROE.dat`-1
#s = 1 + Number of time steps(HLL)
s=`awk -F'\t' 'NR == 1 { print NF; exit }' Density_HLL.dat`-1
#r = 1 + Number of time steps(HLLC)
r=`awk -F'\t' 'NR == 1 { print NF; exit }' Density_HLLC.dat`-1

set term png size 1600,1200
set output "Density.png"
set title 'Density t=0.2 '
set xlabel 'x'
set ylabel 'rho(x,t)'
set xrange [0:1]
set yrange [0.0:1.2]
set grid x y
set style line 1 lc rgb 'red' pt 4   # square
set style line 2 lc rgb 'blue' pt 6   # circle
set style line 3 lc rgb 'green' pt 8   # triangle
set style line 4 lc rgb 'black' pt 12   # tilted square
set style line 5 lc rgb 'brown' pt 3   

plot "Density_LF.dat" using 1:l with points ls 1 title "Lax-Friedrichs",\
	"Density_ROE.dat" using 1:m with points ls 2 title "Roe flux",\
	"Density_HLL.dat" using 1:s with points ls 3 title "HLL flux",\
	"Density_HLLC.dat" using 1:r with points ls 4 title "HLLC flux",\
	"riemann.data" u 1:2 w l title "Exact Solution",\
	"shocktube1D.csv" using 21:11 with linespoints ls 5	 title "CIAO_HLLC"

set term png size 1600,1200
set output "Velocity.png"
set title 'Velocity t=0.2 '
set xlabel 'x'
set ylabel 'u(x,t)'
set xrange [0:1]
set yrange [0.0:1.2]
set grid x y
set style line 1 lc rgb 'red' pt 4   # square
set style line 2 lc rgb 'blue' pt 6   # circle
set style line 3 lc rgb 'green' pt 8   # triangle
set style line 4 lc rgb 'black' pt 12   # tilted square
set style line 5 lc rgb 'brown' pt 3   

plot "Velocity_LF.dat" using 1:l with points ls 1 title "Lax-Friedrichs",\
	"Velocity_ROE.dat" using 1:m with points ls 2 title "Roe flux",\
	"Velocity_HLL.dat" using 1:s with points ls 3 title "HLL flux",\
	"Velocity_HLLC.dat" using 1:r with points ls 4 title "HLLC flux",\
	"riemann.data" u 1:3 w l title "Exact Solution",\
	"shocktube1D.csv" using 21:1 with linespoints ls 5 title "CIAO_HLLC"

set term png size 1600,1200
set output "Pressure.png"
set title 'Pressure t=0.2 '
set xlabel 'x'
set ylabel 'P(x,t)'
set xrange [0:1]
set yrange [0.0:1.2]
set grid x y
set style line 1 lc rgb 'red' pt 4   # square
set style line 2 lc rgb 'blue' pt 6   # circle
set style line 3 lc rgb 'green' pt 8   # triangle
set style line 4 lc rgb 'black' pt 12   # tilted square
set style line 5 lc rgb 'brown' pt 3   

plot "Pressure_LF.dat" using 1:l with points ls 1 title "Lax-Friedrichs",\
	"Pressure_ROE.dat" using 1:m with points ls 2 title "Roe flux",\
	"Pressure_HLL.dat" using 1:s with points ls 3 title "HLL flux",\
	"Pressure_HLLC.dat" using 1:r with points ls 4 title "HLLC flux",\
	"riemann.data" u 1:4 w l title "Exact Solution",\
	"shocktube1D.csv" using 21:12 with linespoints ls 5 title "CIAO_HLLC"


set term png size 1600,1200
set output "Mach_Number.png"
set title 'Mach_number t=0.2 '
set xlabel 'x'
set ylabel 'M(x,t)'
set xrange [0:1]
set yrange [0.0:1.2]
set grid x y
set style line 1 lc rgb 'red' pt 4   # square
set style line 2 lc rgb 'blue' pt 6   # circle
set style line 3 lc rgb 'green' pt 8   # triangle
set style line 4 lc rgb 'black' pt 12   # tilted square
set style line 5 lc rgb 'brown' pt 3   

plot "Mach_Number_LF.dat" using 1:l with points ls 1 title "Lax-Friedrichs",\
	"Mach_Number_ROE.dat" using 1:m with points ls 2 title "Roe flux",\
	"Mach_Number_HLL.dat" using 1:s with points ls 3 title "HLL flux",\
	"Mach_Number_HLLC.dat" using 1:r with points ls 4 title "HLLC flux"	

set term png size 1600,1200
set output "TEnergy.png"
set title 'Total Energy t=0.2 '
set xlabel 'x'
set ylabel 'TE(x,t)'
set xrange [0:1]
set yrange [0.0:3.5]
set grid x y
set style line 1 lc rgb 'red' pt 4   # square
set style line 2 lc rgb 'blue' pt 6   # circle
set style line 3 lc rgb 'green' pt 8   # triangle
set style line 4 lc rgb 'black' pt 12   # tilted square
set style line 5 lc rgb 'brown' pt 3   

plot "TEnergy_LF.dat" using 1:l with points ls 1 title "Lax-Friedrichs",\
	"TEnergy_ROE.dat" using 1:m with points ls 2 title "Roe flux",\
	"TEnergy_HLL.dat" using 1:s with points ls 3 title "HLL flux",\
	"TEnergy_HLLC.dat" using 1:r with points ls 4 title "HLLC flux",\
	"riemann.data" u 1:6 w l title "Exact Solution",\
	"shocktube1D.csv" using 21:($18*$11) with linespoints ls 5 title "CIAO_HLLC"

set term png size 1600,1200
set output "IEnergy.png"
set title 'Internal Energy t=0.2 '
set xlabel 'x'
set ylabel 'IE(x,t)'
set xrange [0:1]
set yrange [0.0:3.5]
set grid x y
set style line 1 lc rgb 'red' pt 4   # square
set style line 2 lc rgb 'blue' pt 6   # circle
set style line 3 lc rgb 'green' pt 8   # triangle
set style line 4 lc rgb 'black' pt 12   # tilted square
set style line 5 lc rgb 'brown' pt 3   

plot "IEnergy_LF.dat" using 1:l with points ls 1 title "Lax-Friedrichs",\
	"IEnergy_ROE.dat" using 1:m with points ls 2 title "Roe flux",\
	"IEnergy_HLL.dat" using 1:s with points ls 3 title "HLL flux",\
	"IEnergy_HLLC.dat" using 1:r with points ls 4 title "HLLC flux",\
	"riemann.data" u 1:5 w l title "Exact Solution",\
	"shocktube1D.csv" using 21:($18-0.5*$1*$1) with linespoints ls 5 title "CIAO_HLLC"
