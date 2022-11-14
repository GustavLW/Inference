function [mi,Si] = get_coeff(xi,x,theta,dt)

dx  = 10^(-7);
mU  = [0;0];
L1a = zeros(2);

for j = 1:size(x,2)
    xj = x(:,j);
    if (xj~=xi)
        U = @(y) V(y,xj,theta);
        mU = mU - U(xi);
        L1a = L1a - [U(xi+dx*[1;0]) - U(xi), U(xi+dx*[0;1]) - U(xi)]/dx;
    end
end
mi  = mU*dt;
Si1 = sqrt(dt)*(eye(2)+dt*L1a/2);
Si2 = dt^(3/2)*L1a/sqrt(12);
Si = Si1'*Si1 + Si2*Si2;
