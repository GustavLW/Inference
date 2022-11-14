function [m_ik,S_ik] = monte_carlo(observed_cells,theta,sigma,i,k,S,dt)
% m is a set of mean vectors, S is a set of variance matrices
% dt can for example be a second or 15 seconds or something.
time_span = observed_cells{end}*60;
m_ik      = zeros(2,S);
S_ik      = zeros(2,2,S);

N     = length(observed_cells)-1;
x     = zeros(2,N);
F     = zeros(2,N);
for ii = 1:N
    x(:,ii) = observed_cells{ii}.location(:,k);
end

k = 0;
for s = 1:S
    while k < time_span
        xi  = x(:,i);
        for j = 1:N
            if j~=i
                xj = x(:,j);
                F(:,i) = F(:,i) + pair_potential(xi,xj,theta(1),theta(2));
            end
        end
        xi = xi + dt*F(:,i) + sqrt(dt)*sigma*randn(2,1);
        k = k + dt;
    end
    
    m_ik(:,s)   = xi + dt*F(:,i);
    S_ik(:,:,s) = sqrt(dt)*sigma*eye(2);
end
