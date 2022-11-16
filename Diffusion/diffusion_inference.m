function [alpha,beta,betaMSD] = diffusion_inference(observed_cells,observed_neighbours,td)

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

    theta   = observed_cells{end};
    base_dt = theta(1);
    theta   = theta(2:end-1);
    alpha   = ones(N,length(td));
    beta    = zeros(N,length(td));
    betaMSD = beta;
for c = 1:length(td)
    skip = td(c);
    dt = base_dt*skip;
    for k   = skip:skip:K-skip
        nei = observed_neighbours{1};
        x   = squeeze(ECEAT(:,:,k))';
        for i = 1:N
            nei_tmp     = nei(i,:);
            xn          = squeeze(ECEAT(i,:,k+skip))';
            xi          = x(:,i);
            [mi,Si]     = get_coeff(xi,x,theta,dt,nei_tmp);
            mi          = xi + mi;
            alpha(i,c)  = alpha(i,c)   + 1;
            beta(i,c)   = beta(i,c)    + 0.5*dot(xn-mi,Si\(xn-mi));
            betaMSD(i,c)= betaMSD(i,c) + 0.5*dot(xn-xi,xn-xi)/dt;
        end
    end
end
