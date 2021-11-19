function [alphaK,betaK] = diffusion_shell(observed_cells,theta,dx,alpha0,beta0)
N = length(observed_cells)-1;
K = length(observed_cells{1}.location);
alphaK = alpha0*ones(N,1);
betaK  = beta0*ones(N,1);
dt = observed_cells{end}/60; % hours between observations
for k = 1:K-1
    nei = neighbours_matrix(observed_cells,k);
    for i = 1:N
        xik  = observed_cells{i}.location(:,k);
        xin  = observed_cells{i}.location(:,k+1);
        aik  = aijk(i,k,nei,observed_cells,dx,theta);
        Aik  = AAijk(i,k,nei,observed_cells,dx,theta);
        Aik  = 0;
        mik  = xik + aik*dt;
        S1ik = sqrt(dt)*(eye(2)+dt/2*Aik);
        S2ik = dt*sqrt(dt/12)*Aik;
        Sik  = S1ik'*S1ik + S2ik'*S2ik;
     
        alphaK(i) = alphaK(i) + 1;
        betaK(i)  = betaK(i) + 0.5*(xin-mik)'*(Sik\(xin-mik));
    end
end