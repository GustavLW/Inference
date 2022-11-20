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
addpath(([pwd '\Interactions']))
df = dir(DataFolder);
df = df(3:end);

d = 1;
load([DataFolder '\' df(d).name])
K = length(observed_cells{end});
A = 2;  % number of optimization agents
T = 10;  % number of generations
L = 3;
S = 4;
agents = cell(A,1);

for a = 1:A
    agents{a} = create_agent(1,T);
end
agents{1}
%%
clc
sim_an1   = @(t) (1/3)*exp(-2/(t^(1/5)));
sim_an2   = @(t) (1/3)*exp(-2/(t^(1/5)));
best_ever_fitness = -Inf;
best_ever_location = agents{randi(a)}.U_best; % arbitrary initiation
for t = 1:T-1
    for a = 1:A
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
    for a = 1:A
        best_tmp(a) = agents{a}.Bfitness;
    end
    best_now = max(best_tmp);
    best_tmp = zeros(1,A);
end
agents{1}