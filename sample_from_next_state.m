function pop_tmp = sample_from_next_state(pop_tmp,U_param,sigma,freq,K)

% input:    population, theta, sigma
% freq:     seconds between observations
% K:        number of simulation steps
nei_tmp     = calculate_neighbours(pop_tmp);
dt          = freq/K;
N           = length(pop_tmp);

for k = 1:K
    F = interactions(pop_tmp,nei_tmp,U_param);
    for i = 1:N
        if ~isnan(pop_tmp(1,i))
            pop_tmp(:,i) = pop_tmp(:,i) + dt*F(:,i) + sqrt(dt)*sigma(i)*[randn;randn];
        end
    end
end


