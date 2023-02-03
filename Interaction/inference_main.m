clc
clear all
close all

TopFolder  = fileparts(fileparts(pwd));
SimFolder  = [TopFolder '\Simulation'];
DataFolder = [SimFolder '\Datasets'];
%dir(SimFolder)
addpath(([SimFolder, filesep]))
addpath(([fileparts(pwd) '\Diffusion']))
addpath(([fileparts(pwd) '\Allee']))
addpath(([fileparts(pwd) '\Interaction']))

df = dir(DataFolder);
df = df(3:end);

load([DataFolder '\' df(end).name])
K     = length(observed_cells{end});
freq  = observed_cells{end-1}(1);
%facit = observed_cells{end-1}(5:end-(length(observed_cells)-2));

A = 2*feature('numcores');  % number of optimization agents
Q = 4;
T = 45;  % number of generations
L = freq/12;
S = 12;
pot_type = 2;

RD = calculate_radial_distribution(observed_cells,1);

repeated_trials = cell(Q,4); % location evolution, fitness, sigma for best location, all agents for good measure

try
    facit = observed_cells{end-1}(5:end-(length(observed_cells)-2));
    s1 = zeros(A,Q,length(facit));
    for q = 1:Q
        for a = 1:A
            s1(a,q,:) = log(facit) +  rand_sphere_shell(length(facit),q/4);
        end
    end
catch
    facit = ones(1,2 + (pot_type-1)*4);
    para_ranges = log(facit') + 1*[-1 1; -1 1; -1 1; -1 1; -1 1; -1 1];
    s0 = hyper_cube_sampling(Q,para_ranges);
    s1 = zeros(A,Q,length(facit));
    for q = 1:Q
        for a = 1:A
            s1(a,q,:) = s0(q,:) + rand_sphere_shell(length(facit),0.5);
        end
    end 
end
%%
tic
for q = 1:Q
    [best_ever_location,best_ever_fitness,best_ever_sigma,agents] =...
    particle_swarm_optimization(q,A,pot_type,T,facit,squeeze(s1(:,q,:)),RD,observed_cells,L,S);
    repeated_trials{q,1} = best_ever_location;
    repeated_trials{q,2} = best_ever_fitness;
    repeated_trials{q,3} = best_ever_sigma;
    repeated_trials{q,4} = agents;
    disp(['Repeated trial ' num2str(q) '/' num2str(Q) ' completed.'])
    toc
end
%%
all_fit = zeros(1,Q);
for q = 1:Q
    all_fit(q) = repeated_trials{q,2}(end);
end
[bfit,bindex] = sort(all_fit);
bindex = bindex(Q-A+1:Q);
besties=zeros(A,2+(pot_type-1)*4);
for a = 1:A
besties(a,:) = repeated_trials{bindex(a),1}(end,:);
end
% todo: pick the A best repeated_trials-results and run PSO again with those
% jag har nu extraherat bästa postionerna från våra Q startgissningar.
% Nästa steg är att skapa EN partikelsvärm som börjar i exakt dessa
% positioner. A<Q, och i "slutspelet" har vi de A bästa partiklarna från
% Q-spelet som får vibe:a

%%
dF = 0;
dW = 0;
close all
for k = 1:K-1
    r  = RD{k};
    rF = RDf{k};
    rW = RDw{k};
    dF = dF + norm(r(:,2)-rF(:,2));
    dW = dW + norm(r(:,2)-rW(:,2));
    hold off
    plot(rF(:,1),rF(:,2),'r','LineWidth',1.5)
    hold on
    plot(rW(:,1),rW(:,2),'b','LineWidth',1.5)
    plot(r(:,1),r(:,2),'k--','LineWidth',1)
    axis([0 5 0 1.5*max(r(:,2))])
    legend('Facit','Proposed','Data')
    drawnow;
    pause(0.1)
end
disp([dF dW])




 function point = rand_sphere_shell(n,r_radii)
        rang = [pi*rand(1,n-2) 2*pi*rand];
        point    = ones(1,n);
        point(2) = sin(rang(1));
        point(end) = prod(sin(rang));
        if n > 3
            for i = 3:n-1
                point(i) = point(i-1).*sin(rang(i-1));
            end
        end
        point(1:end-1) = point(1:end-1).*(cos(rang(1:end)));
        point = r_radii*point;
    end


% %% check RDF deviance for extremely high-fidelity simulation
% clc
% close all
% RD    = calculate_radial_distribution(observed_cells,1);
% freq  = observed_cells{end-1}(1);
% L     = freq;
% forecast_cells = observed_cells;
% U_param          = facit;
% [alpha,beta,~,~] = diffusion_inference(observed_cells,1,U_param);
% N                = (length(observed_cells)-2);
% for k = 1:K-1
%     pop_tmp      = unpack_image(observed_cells,k);
%     nei_tmp      = observed_cells{end}{k};
%     sigma        = sqrt(beta./alpha)';
%     pop_forecast = sample_from_next_state(pop_tmp,nei_tmp,U_param,sigma,freq,L);
%     F = interactions(pop_forecast,U_param,nei_tmp);
%     for i = 1:N
%         forecast_cells{i}.location(:,k+1) = pop_forecast(:,i) + (freq/L)*F(:,i) + sqrt(freq/L)*sigma(i)*[randn;randn];
%     end
% end
% %forecast_cells{end-1}(2:end-3-(length(observed_cells)-2)) = U_param;
% %forecast_cells{end-1}(end-N+1:end)                        = sigma;
% RDf = calculate_radial_distribution(forecast_cells,1);
% 
% %%% PROPOSAL %%%
% q = 1;
% forecast_cells = observed_cells;
% B_param          = exp(repeated_trials{q,1});
% [alpha,beta,~,~] = diffusion_inference(observed_cells,1,B_param);
% sigmaB           = sqrt(beta./alpha)';
% N                = (length(observed_cells)-2);
% for k = 1:K-1
%     pop_tmp      = unpack_image(observed_cells,k);
%     nei_tmp      = observed_cells{end}{k};
%     
%     pop_forecast = sample_from_next_state(pop_tmp,nei_tmp,B_param,sigmaB,freq,L);
%     F = interactions(pop_forecast,B_param,nei_tmp);
%     for i = 1:N 
%         forecast_cells{i}.location(:,k+1) = pop_forecast(:,i) + (freq/L)*F(:,i) + sqrt(freq/L)*sigmaB(i)*[randn;randn];
%     end
% end
% RDw = calculate_radial_distribution(forecast_cells,1);








