function Aik = explicit_hessian(xi,xj,dx,theta)
% phi is a function handle
% xi is me, xj is you. COLUMN VECTORS
% dx is the numerical differentation step

xi_plus1  = xi + [dx;0];
xi_minus1 = xi - [dx;0];
xi_plus2  = xi + [0;dx]; 
xi_minus2 = xi - [0;dx]; 
A11 = (1/dx^2)*(phi(xi_plus1,xj,theta) - 2*phi(xi,xj,theta)...
                + phi(xi_minus1,xj,theta));
            
A22 = (1/dx^2)*(phi(xi_plus2,xj,theta) - 2*phi(xi,xj,theta)...
                + phi(xi_minus2,xj,theta));
            
A12 = (1/dx^2)*(  phi((xi_plus1+xi_plus2)/2,xj,theta)...
                - phi((xi_plus1+xi_minus2)/2,xj,theta)...
                - phi((xi_minus1+xi_plus2)/2,xj,theta)...
                + phi((xi_minus1+xi_minus2)/2,xj,theta));
            
Aik = [A11 A12; A12 A22];
          

