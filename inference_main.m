clc
clear all
datafolder = 'C:\Users\guslindw\Documents\MATLAB\Forskning, höst 2021\Datasets\';
load(strcat(datafolder,'dataset_20221004_1307'))


A = 4;  % number of optimization agents
T = 45; % number of generations
agents = cell(A,1);
dt = 1/3600;
for a = 1:A
    agents{a} = create_agent(dt);
end
%%
% hyperparameters
dx    = 0.01;
theta = [0.0002 4]; % parameters for morse potential
% priors


[alphaK,betaK] = diffusion_shell(observed_cells,theta,dx);
sigma = sqrt(betaK./alphaK)'

histogram(log(sigma),100)

%%
pop_tmp = randn(2,10);
U_param = 1*[1 2 2  0 1 5];
sigma   = 0.001*ones(1,100);
freq = 1200;
K = 30;

S = 20;
k = 1;
ell = synthetic_likelihood(S,k,observed_cells,U_param,sigma,freq,K)

