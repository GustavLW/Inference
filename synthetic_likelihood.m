function [ell,pop_forecast] = synthetic_likelihood(S,k,observed_cells,U_param,sigma,L)
% outer layer!
% S:        number of MCMC samples
% k:        what image are we propagating from?
% input:    population, theta, sigma
% freq:     seconds between observations
% L:        number of simulation steps
% N = length(observed_cells)-2;
ell       = 0;
dt        = observed_cells{end-1}(1); 
pop_facit = unpack_image(observed_cells,k+1);
pop_tmp   = unpack_image(observed_cells,k);
nei_tmp   = observed_cells{end}{k};
pop_forecast  = cell(S,1);

for s = 1:S
    pop_next = sample_from_next_state(pop_tmp,nei_tmp,U_param,sigma,dt,L);
    F = interactions(pop_tmp,(U_param),nei_tmp);
    for i = 1:size(pop_next,2)
        ell = ell + mvnpdf(pop_facit(:,i),pop_next(:,i) + (dt/L)*F(:,i),(dt*sigma(i).^2/L)*eye(2));
        pop_next(:,i) = pop_next(:,i) + (dt/L)*F(:,i) + sqrt(dt/L)*sigma(i)*[randn;randn];
    end
    pop_forecast{s} = pop_next;
end
ell = ell/S;
ell = log(ell);
