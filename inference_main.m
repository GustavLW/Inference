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
K = length(observed_cells{end});
A = 6;  % number of optimization agents
T = 100;  % number of generations
L = 5;
S = 5;
agents   = cell(A,1);
pot_type = 1;
for a = 1:A
    agents{a} = create_agent(pot_type,T);
end
%create a CORRECT agent
%agents{1}.sigma   = observed_cells{end-1}(end)*ones(1,(length(observed_cells)-2));
%agents{1}.U_param = log(observed_cells{end-1}(2:end-3-(length(observed_cells)-2))).*ones(T,length(observed_cells{end-1}(2:end-3-(length(observed_cells)-2))));
%%
clc
sim_an1   = @(t) (1/3)*exp(-7./((2*t/3).^(1/3)));
sim_an2   = @(t) (1/3)*exp(-3./((2*t/3).^(1/7)));
best_ever_fitness = -Inf*ones(T,1);
best_ever_penalty = zeros(T,1);
best_ever_location = zeros(size(agents{randi(a)}.U_best)).*ones(T,1); % arbitrary initiation
tic
for t = 1:T-1
    % agent specific calculations
    parfor a = 1:A
        agents{a} = score_agent(agents{a},t,observed_cells);
        if  agents{a}.fitness(t) - agents{a}.penalty(t) > agents{a}.Bfitness(t) - agents{a}.Bpenalty(t)
            agents{a}.Bfitness(t:T) = agents{a}.fitness(t);
            agents{a}.Bpenalty(t:T) = agents{a}.penalty(t);
            agents{a}.U_best(t:T,:) = agents{a}.U_param(t,:).*ones(T-t+1,1);
        end 
    end
    % propagation
    best_tmp    = zeros(1,A);
    for a = 1:A
        if t > 1
            agents{a}.U_velo(t,:)  = (4/3-sim_an1(t)-sim_an2(t))*agents{a}.U_velo(t,:)...
                                        + sim_an1(t)*(0.5 + 0.0*rand)*(agents{a}.U_best(t,:)   - agents{a}.U_param(t,:))...
                                        + sim_an2(t)*(0.5 + 0.0*rand)*(best_ever_location(t,:) - agents{a}.U_param(t,:));
        end
        agents{a}.U_param(t+1,:) = agents{a}.U_param(t,:) + agents{a}.U_velo(t,:);
        % we also throw this calculation in here
        best_tmp(a) = agents{a}.Bfitness(t) - agents{a}.Bpenalty(t);
    end
    [best_now_fitness,best_index] = max(best_tmp);
    if best_now_fitness > best_ever_fitness(t) - best_ever_penalty(t)
        best_ever_fitness(t+1)    = agents{best_index}.fitness(t);
        best_ever_penalty(t+1)    = agents{best_index}.penalty(t);
        best_ever_location(t+1,:) = agents{best_index}.U_param(t,:);
    else
        best_ever_fitness(t+1)    = best_ever_fitness(t);
        best_ever_penalty(t+1)    = penalize_agent(best_ever_location(t,:),t+1);
        best_ever_location(t+1,:) = best_ever_location(t,:);
    end

    disp([num2str(t) '/' num2str(T)])
toc
end
%%
facit = log(observed_cells{end-1}(2:end-3-(length(observed_cells)-2)));
close all
for t = 1:T
    hold off
    current_params = zeros(A,4*pot_type - 2);
    for a = 1:A
        current_params(a,:) = agents{a}.U_param(t,:);
    end
    scatter(current_params(:,1),current_params(:,2),'r')
    hold on
    scatter(best_ever_location(t,1),best_ever_location(t,2),'b')
    plot(facit(1),facit(2),'ko')
    axis([-15 -5 -5 5])
    title(num2str(t))
    drawnow;
    pause(0.000001)
end

%%
subplot(2,1,1)
for a = 1:A
plot(agents{a}.Bfitness,'b')
hold on
plot(agents{a}.fitness,'m')
plot(best_ever_fitness,'r')
%axis([1 T -abs(max(best_ever_fitness)) 0.9*abs(max(best_ever_fitness))])
end
subplot(2,1,2)
for a = 1:A
plot((agents{a}.penalty(2:T-1)),'b')
hold on
plot((best_ever_penalty(2:T-1)),'r')
%axis([1 T -max(best_ever_penalty) 1.1*max(best_ever_penalty)])
end

%%
r  = linspace(0,4,4001);
th = exp(best_ever_location(T-1,:));
y  = U(r,th);
z  = U(r,exp(facit));
plot(r,y,'b')
axis([0 4 -2*exp(facit(1)) 4*exp(facit(1))])
hold on
plot(r,z,'r')
%%
clc
clf
RD = calculate_radial_distribution(observed_cells);

subplot(2,1,1)
plot(RD(:,1),RD(:,2)/mean(RD(401:end,2)),'b')
forecast_cells = observed_cells;
U_param          = exp(facit);
[alpha,beta,~,~] = diffusion_inference(observed_cells,1,U_param);
N                = (length(observed_cells)-2);
for k = 1:K-1
    pop_tmp      = unpack_image(observed_cells,k);
    nei_tmp      = observed_cells{end}{k};
    sigma        = sqrt(beta./alpha)';
    pop_forecast = sample_from_next_state(pop_tmp,nei_tmp,U_param,sigma,1,observed_cells{end-1}(1));
    for i = 1:N 
        forecast_cells{i}.location(:,k+1) = pop_forecast(:,i);
    end
end
forecast_cells{end-1}(2:end-3-(length(observed_cells)-2)) = U_param;
forecast_cells{end-1}(end-N+1:end)                        = sigma;
RDf = calculate_radial_distribution(forecast_cells);
hold on
plot(RDf(:,1),RDf(:,2)/mean(RDf(401:end,2)),'r')

sim_PCF = RDf(:,2)/mean(RDf(401:end,2));
tru_PCF = RD(:,2)/mean(RD(401:end,2));
norm(sim_PCF-tru_PCF)

subplot(2,1,2)
plot(RD(:,1),RD(:,2)/mean(RD(401:end,2)),'b')
forecast_cells = observed_cells;
hold on
U_param          = exp(best_ever_location(end,:));
[alpha,beta,~,~] = diffusion_inference(observed_cells,1,U_param);
N                = (length(observed_cells)-2);
for k = 1:K-1
    pop_tmp      = unpack_image(observed_cells,k);
    nei_tmp      = observed_cells{end}{k};
    sigma        = sqrt(beta./alpha)';
    pop_forecast = sample_from_next_state(pop_tmp,nei_tmp,U_param,sigma,1,observed_cells{end-1}(1));
    for i = 1:N 
        forecast_cells{i}.location(:,k+1) = pop_forecast(:,i);
    end
end
forecast_cells{end-1}(2:end-3-(length(observed_cells)-2)) = U_param;
forecast_cells{end-1}(end-N+1:end)                        = sigma;
RDf = calculate_radial_distribution(forecast_cells);
hold on
plot(RDf(:,1),RDf(:,2)/mean(RDf(401:end,2)),'r')

sim_PCF = RDf(:,2)/mean(RDf(401:end,2));
tru_PCF = RD(:,2)/mean(RD(401:end,2));
norm(sim_PCF-tru_PCF)

%%
for i = 1:size(pop_next,2)
    ell = ell + mvnpdf(pop_facit(:,i),pop_next(:,i),(dt*sigma(i).^2/L)*eye(2));
end

% function pop_tmp = sample_from_next_state(pop_tmp,nei_tmp,U_param,sigma,freq,L)
% % has passed stage 1 testing yee always wins
% % input:    population, neighbours (matrix), theta, sigma
% % freq:     seconds between observations
% % K:        number of simulation steps
% dt          = freq/L;
% N           = size(pop_tmp,2);
% 
% for k = 1:L-1
%     F = interactions(pop_tmp,U_param,nei_tmp); %!! interactions(population,U_param)
%     for i = 1:N
%         if ~isnan(pop_tmp(1,i))
%             pop_tmp(:,i) = pop_tmp(:,i) + dt*F(:,i) + sqrt(dt)*sigma(i)*[randn;randn];
%         end
%     end
% end