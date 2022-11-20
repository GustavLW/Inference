function agent = score_agent(agent,t,observed_cells,S,L)
[alpha,beta,~,~] = diffusion_inference(observed_cells,1,agent.U_param(t,:)); % works beautifully!
agent.sigma  = sqrt(beta./alpha)';
ell = 0;

for k = 1:length(observed_cells{end})-1
    ell = ell + synthetic_likelihood(S,k,observed_cells,agent.U_param(t,:),agent.sigma,L);
end

agent.fitness  = ell  - penalize_agent(agent.U_param(t,:),t);
agent.Bfitness = agent.Bfitness - penalize_agent(agent.U_best,t)...
                                + penalize_agent(agent.U_best,max(t-1,1));






