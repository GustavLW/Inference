function [best_ever_location,best_ever_fitness,best_ever_sigma,agents] = particle_swarm_optimization(q,A,pot_type,T,facit,s1,RD,observed_cells,L,S)
agents   = cell(A,1);
sim_an1  = @(t) (1/3)*(1./(1+exp(2-7*t/T)));
sim_an2  = @(t) (1/3)*(1./(1+exp(3-11*t/T)));
inertia  = (1:A)./((1:A));

for a = 1:A
    agents{a}          = create_agent(pot_type,T);
    agents{a}.U_param  = s1(a,:).*ones(T,length(facit));
end
best_ever_fitness  = -Inf*ones(T,1);
best_ever_penalty  = zeros(T,1);
best_ever_location = zeros(size(agents{randi(a)}.U_best)).*ones(T,1); % arbitrary initiation
best_ever_sigma    = zeros(T,length(observed_cells)-2);
for t = 1:T-1
    % agent specific calculations
    parfor a = 1:A
        [agent,~] = score_agent(RD,agents{a},t,observed_cells,L,S);
        agents{a} = agent;
        if  agents{a}.fitness(t) - agents{a}.penalty(t) > agents{a}.Bfitness(t) - agents{a}.Bpenalty(t)
            agents{a}.Bfitness(t:T) = agents{a}.fitness(t);
            agents{a}.Bpenalty(t:T) = agents{a}.penalty(t);
            agents{a}.U_best(t:T,:) = agents{a}.U_param(t,:).*ones(T-t+1,1);
        end
    end
    % propagation
    best_tmp = zeros(1,A);
    for a = 1:A
        if t > 1
            agents{a}.U_velo(t,:)  = 2*inertia(a)*(1 - sim_an1(t) - sim_an2(t))*agents{a}.U_velo(t,:)...
                + 0.25*(1/5000 + sim_an1(t))*(agents{a}.U_best(t,:)   - agents{a}.U_param(t,:))...
                + 0.25*(1/5000 + sim_an2(t))*(best_ever_location(t,:) - agents{a}.U_param(t,:));
        end
        vmax = 0.15*length(agents{a}.U_velo(t,:));
        if norm(agents{a}.U_velo(t,:)) > vmax
            agents{a}.U_velo(t,:) = vmax*agents{a}.U_velo(t,:)/norm(agents{a}.U_velo(t,:));
        end
        agents{a}.U_param(t+1,:) = agents{a}.U_param(t,:) + agents{a}.U_velo(t,:);
        best_tmp(a)              = agents{a}.Bfitness(t)  - agents{a}.Bpenalty(t); % we also throw this calculation in here
    end
    [best_now_fitness,best_index] = max(best_tmp);
    if best_now_fitness > best_ever_fitness(t) - best_ever_penalty(t)
        best_ever_fitness(t+1)    = agents{best_index}.Bfitness(t);
        best_ever_penalty(t+1)    = agents{best_index}.Bpenalty(t);
        best_ever_location(t+1,:) = agents{best_index}.U_best(t,:);
        best_ever_sigma(t+1,:)    = agents{best_index}.sigma;
    else
        best_ever_fitness(t+1)    = best_ever_fitness(t);
        best_ever_location(t+1,:) = best_ever_location(t,:);
        best_ever_sigma(t+1,:)    = best_ever_sigma(t,:);
        best_ever_penalty(t+1)    = penalize_agent(RD,0,best_ever_location(t+1,:),t);
    end
    %disp(['Iteration ' num2str(t) '/' num2str(T) ' of trial ' num2str(q) '/' num2str(Q) '.'])
end
end