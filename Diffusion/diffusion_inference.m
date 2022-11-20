function [alpha,beta,betaMSD,betaIMP] = diffusion_inference(observed_cells,td,U_param)

sample_cell = observed_cells{1};
K           = length(sample_cell.location);
clear sample_cell
N           = length(observed_cells)-2;
ECEAT       = NaN*ones(N,2,K);
for i = 1:N
    for k = 1:K
        if and(k >= observed_cells{i}.b_time,k < observed_cells{i}.d_time)
              ECEAT(i,:,k) = observed_cells{i}.location(:,k);
        end
    end
end

base_dt = observed_cells{end-1}(1);
switch nargin
    case 2
        theta   = observed_cells{end-1}(2:end-3-N); 
    case 3
        theta = U_param;
end
alpha   = zeros(N,length(td));
beta    = zeros(N,length(td));
betaMSD = beta;
betaIMP = beta;
observed_neighbours = observed_cells{end};
for c = 1:length(td)
    skip = td(c);
    dt = base_dt*skip;
    for k   = skip:skip:K-skip
        nei = observed_neighbours{k};
        fei = observed_neighbours{k+skip};
        x   = squeeze(ECEAT(:,:,k))';      % some are NaN
        axn = squeeze(ECEAT(:,:,k+skip))'; % some are NaN
        for i = 1:N
            if and(observed_cells{i}.b_time <= k,observed_cells{i}.d_time > k)
                if ~isnan(axn(:,i))
                    nei_tmp      = nei(i,:);
                    fei_tmp      = fei(i,:);
                    xi           = x(:,i);
                    xn           = axn(:,i);
                    [mi,Si,Zi]   = get_coeff(xi,x,theta,dt,nei_tmp);
                    [ni,~,~]     = get_coeff(xn,axn,theta,dt,fei_tmp);
                    ni           = xi + 0.5*(mi+ni);
                    mi           = xi + mi;
                    alpha(i,c)   = alpha(i,c)   + 1;
                    beta(i,c)    = beta(i,c)    + 0.5*dot(xn-mi,Si\(xn-mi));
                    betaIMP(i,c) = betaIMP(i,c) + 0.5*dot(xn-ni,Zi\(xn-ni));
                    betaMSD(i,c) = betaMSD(i,c) + 0.5*dot(xn-xi,xn-xi)/dt;
                end
            end
        end
    end
end
