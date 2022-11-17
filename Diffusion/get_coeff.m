function [mi,Si,Zi] = get_coeff(xi,x,theta,dt,nei)

dx  = 10^(-8);
mU  = [0;0];
L1a = zeros(2);
nei_tmp = find(nei>0);
if ~isempty(nei_tmp)
    for n = 1:length(nei_tmp)
        j  = nei_tmp(n);
        xj = x(:,j);
        U = @(y) pair_potential(y,xj,theta);
        if (xj~=xi)
            mU  = mU + U(xi);
            L1a = L1a + [U(xi+dx*[1;0]) - U(xi), U(xi+dx*[0;1]) - U(xi)]/dx;
        end
    end
end

mi  = mU*dt;
Si1 = sqrt(dt)*(eye(2)+dt*L1a/2);
Zi1 = sqrt(dt)*eye(2);
Si2 = dt^(3/2)*L1a/sqrt(12);
Si = Si1'*Si1 + Si2'*Si2;
Zi = Zi1'*Zi1 + Si2'*Si2;
