reset

ftitle(a,b) = sprintf("%.0f/x^{%4.2f}", exp(a), -b)

f(x)=a+b*x
fit f(x) 'cout' u (log($1)):(log($4)) via a,b
f2(x)=a2+b2*x
fit f2(x) 'cout' u (log($1)):(log($2)) via a2,b2

fc(x)=ac+bc*x
fit fc(x) 'cout2' u (log($1)):(log($4)) via ac,bc
fc2(x)=ac2+bc2*x
fit fc2(x) 'cout2' u (log($1)):(log($2)) via ac2,bc2

set xlabel 'Maximum resolution'
set ylabel 'Maximum error'
set key bottom left
set logscale
set xrange [16:256]
set xtics 16,2,256
set grid ytics
set cbrange [1:1]
plot 'cout' u 1:4 t 'max (adaptive)', exp(f(log(x))) t ftitle(a,b), \
     'cout2' u 1:4 t 'max (constant)', exp(fc(log(x))) t ftitle(ac,bc), \
     'cout' u 1:2 t 'norm1 (adaptive)', exp(f2(log(x))) t ftitle(a2,b2), \
     'cout2' u 1:2 t 'norm1 (constant)', exp(fc2(log(x))) t ftitle(ac2,bc2)

if (1) set term png; set output "interface.png"; else pause -1;
reset
set size ratio -1
plot [-0.5:0.5][-0.5:0.5]'cout' w l t "adaptive", 'cout2' w l t "constant"
