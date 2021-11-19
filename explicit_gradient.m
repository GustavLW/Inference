function aik = explicit_gradient(xi,xj,dx,theta)
dx = dx/2;
xi_plus1  = xi + [dx;0];
xi_minus1 = xi - [dx;0];
xi_plus2  = xi + [0;dx]; 
xi_minus2 = xi - [0;dx]; 


a1 = 1/(2*dx)*(phi(xi_plus1,xj,theta)-phi(xi_minus1,xj,theta));
a2 = 1/(2*dx)*(phi(xi_plus2,xj,theta)-phi(xi_minus2,xj,theta));

aik = [a1; a2];