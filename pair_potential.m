function f = pair_potential(xi,xj,U_param) 
% calculate the influence of cell j on cell i 
r    = norm(xi-xj);
if r < 3 % extremely arbitrary
    f = -(xi-xj)*u(r,U_param)/r;
else
    f = [0;0];
end
