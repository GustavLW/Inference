function [alphaK,betaK] = diffusion_shell(observed_cells,theta,dx)

sample_cell = observed_cells{1};
K           = length(sample_cell.location);
clear sample_cell
N           = length(observed_cells)-1;
ECEAT       = NaN*ones(N,2,K);
for i = 1:N
    for k = 1:K
        if k >= observed_cells{i}.b_time
            ECEAT(i,:,k) = observed_cells{i}.location(:,k);
        end
    end
end
dt = observed_cells{end};
alpha0 = ones(N,1);
beta0  = zeros(N,1);



for k   = 1:K-1
    x   = squeeze(ECEAT(:,:,k))';
    axn = squeeze(ECEAT(:,:,k+1))';
    for i = 1:N
        xn         = squeeze(ECEAT(i,:,k+1))';
        xi         = x(:,i);
        [mi,Si]    = get_coeff(xi,x,theta,dt,dx);     
        mi         = xi + mi;
        alpha0(i)  = alpha0(i) + 1;
        beta0(i)   = beta0(i) + 0.5*dot(xn-mi,Si\(xn-mi));
   
    end
end

alphaK = alpha0;
betaK = beta0;