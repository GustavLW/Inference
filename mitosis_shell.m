function mitosis_shell(observed_cells)

% vi behöver hämta ut   beta_ik: födda från i vid tid k
%                       rho_ik:  densitet runt i vid tid k

N = length(observed_cells)-1;
K = length(observed_cells{1}.location);

rho  = zeros(N,K);
kernel_window = 2.140162717247335;
% räkna ut rho
for i = 1:N
    for k = 1:K
        if k >= observed_cells{i}.b_time
            for j = 1:N
                if j ~= i
                    if k >= observed_cells{j}.b_time
                        distance = norm(observed_cells{i}.location(:,k)-...
                            observed_cells{j}.location(:,k));
                        rho(i,k) = rho(i,k) + exp(-kernel_window*distance); % a bit sloppy
                    end
                end
            end
        end
    end
end
% räkna ut beta
beta = zeros(N,K);
for i = 1:N
   child_tmp = observed_cells{i}.children;
   for j = 1:length(child_tmp)
       k = observed_cells{child_tmp(j)}.b_time;
       if k < K-1
       beta(i,k+1) = beta(i,k+1) + 1; % 0 is 1, 1 is 2... fucking matlab
       end
   end
end

dt = observed_cells{end}/60;

% SQP, constraining to positive lambda0,lambda1

clc
M  = 2^10;
%l0 = rand;
%l1 = rand/max(max(rho));
%m0 = 0.1*rand;
%m1 = 0.1*rand;
l0 = 0.25;
l1 = 0.6;
m0 = 0;
m1 = 0;
alpha = 0.000001;
% time step in hours
%

for m = 1:M
    likelihood      = 0;
    grad            = [0;0];
    hess            = zeros(2);
    for i=1:N
        for k = 1:K
            if k >= observed_cells{i}.b_time
                likelihood = likelihood + beta(i,k)*log(l0*(1-l1*rho(i,k))*dt)...
                            + l0*(1-l1*rho(i,k))*dt;
                grad = grad...
                    + [beta(i,k)/l0-(1-l1*rho(i,k))*dt;...
                    -beta(i,k)*rho(i,k)/(1-l1*rho(i,k))+l0*rho(i,k)*dt];
                hess = hess...
                    + [-beta(i,k)/l0^2 rho(i,k)*dt; rho(i,k)*dt...
                    -beta(i,k)*rho(i,k)^2/(1-l1*rho(i,k))^2];
            end
        end
    end
    if m == 1
        lili = likelihood;
    end
    
    grad_likelihood = [grad-[m0;m1];m0*l0;m1*l1];
    hess_likelihood = [hess, -eye(2);m0 0 l0 0; 0 m1 0 l1];
    pm     = (hess_likelihood-400*eye(4))\grad_likelihood;
    l_next = [l0;l1;m0;m1] - alpha*pm;
    l0 = l_next(1);
    l1 = l_next(2);
    m0 = l_next(3);
    m1 = l_next(4);
    if mod(log2(m),2) == 0
        disp(['Generation: ' num2str(m)])
        disp(['Likelihood: ' num2str(likelihood) '. Evolution from beginning: ' num2str(likelihood-lili)])
        disp(['Current position: [' num2str(l_next') ']'])
        disp(['Search direction: [' num2str(pm') ']'])
    end
    alpha = 100*(atan(32*m/M)./m.^(6/7)).^(7/3);
end




