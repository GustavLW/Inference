function f = pair_potential(xi,xj,U_param) 
% calculate the influence of cell j on cell i 
r    = norm(xi-xj);
disp('check 2')
if r < 3 % extremely arbitrary
    f = -(xi-xj)*u(r,U_param)/r;
    g = V(xi,xj,U_param);
    disp(norm(f-g))
else
    f = [0;0];
end
