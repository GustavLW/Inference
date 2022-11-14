function ell = synthetic_likelihood(S,k,observed_cells,U_param,sigma,freq,K)

% S:        number of MCMC samples
% k:        what image are we propagating from?
% input:    population, theta, sigma
% freq:     seconds between observations
% K:        number of simulation steps
%N = length(observed_cells)-1;
ell = 0;
pop_facit = unpack_image(observed_cells,k+1);
for s = 1:S
    pop_tmp = unpack_image(observed_cells,k);
    pop_next = sample_from_next_state(pop_tmp,U_param,sigma,freq,K);
    for i = 1:size(pop_next,2)
        ell = ell + mvnpdf(pop_facit(:,i),pop_next(:,i),(freq*sigma(i).^2/K)*eye(2));
    end
end
ell = ell/S;
ell = log(ell);