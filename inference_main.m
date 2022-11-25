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

load([DataFolder '\' df(end).name])
K     = length(observed_cells{end});
freq  = observed_cells{end-1}(1);
facit = (observed_cells{end-1}(2:end-3-(length(observed_cells)-2)));

A = 2*feature('numcores');  % number of optimization agents
Q = 1;
T = 180;  % number of generations
L = freq/12;
S = 4;
agents   = cell(A,1);
pot_type = 1;

RD = calculate_radial_distribution(observed_cells,1);
sim_an1   = @(t) (1/3)*exp(-4./((12*t).^(1/4)));
sim_an2   = @(t) (1/3)*exp(-2./((9*t).^(1/6)));

repeated_trials = cell(Q,2);
%%
clc
tic
for q = 1:Q
    for a = 1:A
        agents{a} = create_agent(pot_type,T);
    end
    %agents{1}.U_param = log(facit).*ones(T,length(facit));
    best_ever_fitness  = -Inf*ones(T,1);
    best_ever_penalty  = zeros(T,1);
    best_ever_dev      = zeros(T,1);
    best_ever_location = zeros(size(agents{randi(a)}.U_best)).*ones(T,1); % arbitrary initiation
    best_ever_index    = 0;
    for t = 1:T-1
        % agent specific calculations
        parfor a = 1:A
            [agent,~] = score_agent(RD,agents{a},t,observed_cells,L,S);
            agents{a} = agent;
            if  agents{a}.fitness(t) - agents{a}.penalty(t) > agents{a}.Bfitness(t) - agents{a}.Bpenalty(t)
                agents{a}.Bfitness(t:T) = agents{a}.fitness(t);
                agents{a}.Bpenalty(t:T) = agents{a}.penalty(t);
                agents{a}.Bdev(t:T)     = agents{a}.dev(t);
                agents{a}.U_best(t:T,:) = agents{a}.U_param(t,:).*ones(T-t+1,1);
            end
        end
        % propagation
        inertia  = 1:A;
        best_tmp = zeros(1,A);
        for a = 1:A
            if t > 1
                agents{a}.U_velo(t,:)  = inertia(a)*(1 - sim_an1(t) - sim_an2(t))*agents{a}.U_velo(t,:)...
                    + 0.2*(4+rand)*(1/250 + sim_an1(t))*(agents{a}.U_best(t,:)   - agents{a}.U_param(t,:))...
                    + 0.2*(4+rand)*(1/250 + sim_an2(t))*(best_ever_location(t,:) - agents{a}.U_param(t,:));
            end
            vmax = 0.075;
            if norm(agents{a}.U_velo(t,:)) > vmax
                agents{a}.U_velo(t,:) = vmax*agents{a}.U_velo(t,:)/norm(agents{a}.U_velo(t,:));
            end
            agents{a}.U_param(t+1,:) = agents{a}.U_param(t,:) + agents{a}.U_velo(t,:);
            best_tmp(a)              = agents{a}.Bfitness(t)  - agents{a}.Bpenalty(t); % we also throw this calculation in here
        end
        [best_now_fitness,best_index] = max(best_tmp);
        if best_now_fitness > best_ever_fitness(t) - best_ever_penalty(t)
            best_ever_fitness(t+1)    = agents{best_index}.fitness(t);
            best_ever_penalty(t+1)    = agents{best_index}.penalty(t);
            best_ever_dev(t+1)        = agents{best_index}.dev(t);
            best_ever_location(t+1,:) = agents{best_index}.U_best(t,:);
            best_ever_index           = best_index;
            new_best                  = 0;
        else
            new_best                  = new_best + 1;
            best_ever_fitness(t+1)    = best_ever_fitness(t);
            best_ever_dev(t+1)        = best_ever_dev(t)*1.00^new_best; % in the case of spurious mishaps, these guys get less valuable as time goes on
            best_ever_location(t+1,:) = best_ever_location(t,:);
            best_ever_penalty(t+1)    = penalize_agent(RD,best_ever_dev(t+1),best_ever_location(t+1,:),t+1);
        end

        disp(['Iteration ' num2str(t) '/' num2str(T) ' of trial ' num2str(q) '/' num2str(Q) '.'])
        toc
    end
    repeated_trials{q,1} = exp(best_ever_location(end,:));
    repeated_trials{q,2} = agents{best_ever_index}.sigma;
end
%% check RDF deviance for extremely high-fidelity simulation
clc
close all
RD    = calculate_radial_distribution(observed_cells,1);
freq  = observed_cells{end-1}(1);
L     = freq;
forecast_cells = observed_cells;
U_param          = facit;
[alpha,beta,~,~] = diffusion_inference(observed_cells,1,U_param);
N                = (length(observed_cells)-2);
for k = 1:K-1
    pop_tmp      = unpack_image(observed_cells,k);
    nei_tmp      = observed_cells{end}{k};
    sigma        = sqrt(beta./alpha)';
    pop_forecast = sample_from_next_state(pop_tmp,nei_tmp,U_param,sigma,freq,L);
    F = interactions(pop_forecast,U_param,nei_tmp);
    for i = 1:N
        forecast_cells{i}.location(:,k+1) = pop_forecast(:,i) + (freq/L)*F(:,i) + sqrt(freq/L)*sigma(i)*[randn;randn];
    end
end
forecast_cells{end-1}(2:end-3-(length(observed_cells)-2)) = U_param;
forecast_cells{end-1}(end-N+1:end)                        = sigma;
RDf = calculate_radial_distribution(forecast_cells,1);

%%% PROPOSAL %%%

forecast_cells = observed_cells;
B_param          = repeated_trials{1,1};
[alpha,beta,~,~] = diffusion_inference(observed_cells,1,B_param);
sigmaB           = sqrt(beta./alpha)';
N                = (length(observed_cells)-2);
for k = 1:K-1
    pop_tmp      = unpack_image(observed_cells,k);
    nei_tmp      = observed_cells{end}{k};
    
    pop_forecast = sample_from_next_state(pop_tmp,nei_tmp,B_param,sigmaB,freq,L);
    F = interactions(pop_forecast,B_param,nei_tmp);
    for i = 1:N 
        forecast_cells{i}.location(:,k+1) = pop_forecast(:,i) + (freq/L)*F(:,i) + sqrt(freq/L)*sigmaB(i)*[randn;randn];
    end
end
RDw = calculate_radial_distribution(forecast_cells,1);
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
%% COMPARE WINNING POTENTIAL COMPARED TO UNDERLYING
for q = 1:Q
    hold off
    x       = linspace(0,5,1001);
    xs      = linspace(0,5,21);
    facit   = observed_cells{end-1}(2:end-3-(length(observed_cells)-2));
    B_param = (repeated_trials{q,1});
    plot(xs,U(xs,facit),'ro')
    hold on
    plot(xs+0.125,U(xs+0.125,B_param),'bd')
    plot(x,U(x,facit),'r')
    plot(x,U(x,B_param),'b')
    plot([0.5 4],[0 0],'k')
    grid on
    axis([0.5 4 facit(1)*[-2 5]])
    xlabel('Distance')
    ylabel('Potential energy')
    title('Winning potential compared to underlying')
    legend('Underlying potential','Propsed by inference')
    drawnow;
    pause(1)
end
%% GRAFIK FÖR MORSE

close all
for t = 1:1:T
    current_params = zeros(A,4*pot_type - 2);
    current_best   = current_params;
    for a = 1:A
        current_params(a,:) = agents{a}.U_param(t,:);
        current_best(a,:)   = agents{a}.U_best(t,:);
    end
    hold off
    scatter(current_params(:,1),current_params(:,2),'k*')
    hold on
    scatter(current_best(:,1),current_best(:,2),'b*')
    scatter(best_ever_location(t,1),best_ever_location(t,2),'ro')
    plot(log(facit(1)),log(facit(2)),'rd')
    axis([-15 -0 0 5])
    grid on
    drawnow;
    pause(0.01)
end



