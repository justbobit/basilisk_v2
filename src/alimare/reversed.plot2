reset

set term png; set output "ordre_reversed.png";

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)

f(x)=a+b*x
fit f(x) 'log' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'log' u (log($1)):(log($2)) via a2,b2

f3(x)=a3+b3*x
fit f3(x) 'log' u (log($1)):(log($7)) via a3,b3
f4(x)=a4+b4*x
fit f4(x) 'log' u (log($1)):(log($5)) via a4,b4

set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set key bottom left
set logscale
set xrange [16:256]
set xtics 16,2,256
set grid ytics
set cbrange [1:1]

plot 'log' u 1:4 t 'VOF max (cart)'  , exp(f(log(x)))       t ftitle(a ,b ), \
     'log' u 1:2 t 'VOF norm1 (cart)', exp(f2(log(x))) lw 3 t ftitle(a2,b2), \
     'log' u 1:7 t 'LS max (cart)'   , exp(f3(log(x)))      t ftitle(a3,b3), \
     'log' u 1:5 t 'LS norm1 (cart)' , exp(f4(log(x))) lw 3 t ftitle(a4,b4)