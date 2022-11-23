function agent = score_agent(agent,t,observed_cells,S,L)
[alpha,beta,~,~] = diffusion_inference(observed_cells,1,exp(agent.U_param(t,:))); % works beautifully!
agent.sigma  = sqrt(beta./alpha)';

switch nargin 
    case 3
agent.fitness(t) = -compare_radial_distribution(observed_cells,exp(agent.U_param(t,:)));
    case 5
        ell = 0;
        for k = 1:length(observed_cells{end})-1
            ell = ell + synthetic_likelihood(S,k,observed_cells,exp(agent.U_param(t,:)),agent.sigma,L);
        end
        agent.fitness(t)  = ell;
end
agent.penalty(t)  = penalize_agent(agent.U_param(t,:),t);
if t > 1
    agent.Bpenalty(t) = penalize_agent(agent.U_best(t,:),t);
end





