clc
clear all
close all

TopFolder  = fileparts(pwd);
SimFolder  = [TopFolder '\Simulation'];
DataFolder = [SimFolder '\Datasets'];
%dir(SimFolder)
addpath(([SimFolder, filesep]))
addpath(([pwd '\Diffusion']))
addpath(([pwd '\Allee']))
addpath(([pwd '\Interaction']))
df = dir(DataFolder);
df = df(3:end);

d = 1;
load([DataFolder '\' df(d).name])
K = length(observed_cells{end});
A = 4;  % number of optimization agents
T = 50;  % number of generations
L = 3;
S = 30;
agents   = cell(A,1);
pot_type = 1;
for a = 1:A
    agents{a} = create_agent(pot_type,T);
end
% create a CORRECT agent
agents{1}.sigma   = observed_cells{end-1}(end)*ones(1,(length(observed_cells)-2));
agents{1}.U_param = observed_cells{end-1}(2:end-3-(length(observed_cells)-2)).*ones(T,length(observed_cells{end-1}(2:end-3-(length(observed_cells)-2))));
%%
clc
sim_an1   = @(t) (1/3)*exp(-2/(t^(1/5)));
sim_an2   = @(t) (1/3)*exp(-2/(t^(1/5)));
best_ever_fitness = -Inf;
best_ever_location = agents{randi(a)}.U_best; % arbitrary initiation
tic
for t = 1:T-1
    parfor a = 1:A
        agents{a} = score_agent(agents{a},t,observed_cells,S,L);
        if agents{a}.fitness > agents{a}.Bfitness
            agents{a}.Bfitness =  agents{a}.fitness;
            agents{a}.U_best = agents{a}.U_param(t,:);
        end
        agents{a}.U_velo         = (1-sim_an1(t)-sim_an2(t))*agents{a}.U_velo...
                                    + sim_an1(t)*(agents{a}.U_best   - agents{a}.U_param(t,:))...
                                    + sim_an2(t)*(best_ever_location - agents{a}.U_param(t,:));
        agents{a}.U_param(t+1,:) = agents{a}.U_param(t,:) + agents{a}.U_velo;
    end
    best_tmp = zeros(1,A);
    for a = 1:A
        best_tmp(a) = agents{a}.Bfitness;
    end
    [best_now_fitness,best_index] = max(best_tmp);
    if best_now_fitness > best_ever_fitness
        best_ever_fitness = best_now_fitness;
        best_ever_location = agents{best_index}.U_param(t,:);
    end

    disp([num2str(t) '/' num2str(T)])
toc
end
%%
facit = observed_cells{end-1}(2:end-3-(length(observed_cells)-2));

for t = 1:T
    hold off
    current_params = zeros(A,4*pot_type - 2);
    for a = 1:A
        current_params(a,:) = agents{a}.U_param(t,:);
    end
    scatter(current_params(:,1),current_params(:,2),'r')
    hold on
    plot(facit(1),facit(2),'ko')
    axis([0 5*facit(1) 0 5*facit(2)])
    drawnow;
    pause(0.5)
end


%%
winner_score = ones(1,T);
winner_penal = winner_score;
agent = agents{1};
tic
for t = 1:T
    ell = 0;
    for k = 1:length(observed_cells{end})-1
        ell = ell + synthetic_likelihood(S,k,observed_cells,agent.U_param(t,:),agent.sigma,L);
    end
    winner_score(t) = ell;
    winner_penal(t) = penalize_agent(agent.U_param(t,:),t);
    disp([num2str(t) '/' num2str(T)])
    toc
end
