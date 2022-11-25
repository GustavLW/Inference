function [agent,forecast_array] = score_agent(RD,agent,t,observed_cells,L,S)

forecast_array    = cell(S,1);
for s = 1:S
forecast_array{s} = observed_cells;
end

[alpha,beta,~,~] = diffusion_inference(observed_cells,1,exp(agent.U_param(t,:))); % works beautifully!
agent.sigma      = sqrt(beta./alpha)';
ell = 0;
for k = 1:length(observed_cells{end})-1
    [elt,pop_forecast] = synthetic_likelihood(S,k,observed_cells,exp(agent.U_param(t,:)),agent.sigma,L);
    ell = ell + elt;
    for s = 1:S
        tmp_array = forecast_array{s};
        for i = 1:length(observed_cells)-2
            tmp_array{i}.location(:,k+1) = pop_forecast{s}(:,i);
        end
        if k == length(observed_cells{end})-1
            tmp_array{end-1}(2:end-3-(length(observed_cells)-2)) = exp(agent.U_param(t,:));
            tmp_array{end-1}(end-(length(observed_cells)-2)+1:end) = agent.sigma;
            tmp_array{end-1}(1) = observed_cells{end-1}(1)/L;
        end
        tmp_array{end}    = length(observed_cells{end});
        forecast_array{s} = tmp_array;
    end
end
[penalty,dev] = penalize_agent(RD,forecast_array,agent.U_param(t,:),t);

agent.fitness(t)  = ell;
agent.penalty(t)  = penalty;
agent.dev(t)      = dev;
if t > 1
    agent.Bpenalty(t) = penalize_agent(RD,agent.Bdev(t),agent.U_best(t,:),t);
end





