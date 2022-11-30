clc
clear all
close all

TopFolder  = fileparts(fileparts(pwd));
SimFolder  = [TopFolder '\Simulation'];
DataFolder = [SimFolder '\Datasets'];
%dir(SimFolder)
addpath(([SimFolder, filesep]))
df = dir(DataFolder);
df = df(3:end);

    o0  = 0;
    o1  = 0;

d = length(df);
load([DataFolder '\' df(d).name])

[ECAET,birth_event_list,death_event_list] = get_locations_and_birth_events(observed_cells,1);

D             = distances(ECAET);
kernel_window = 3.068891070360678;
% räkna ut rho
[rho,~,~]     = density(D,kernel_window,0);
% räkna ut beta
N   = size(D{1},1);
K   = length(D);
dbeta  = zeros(N,K);
for b = 1:size(birth_event_list,1)
    dbeta(birth_event_list(b,1),birth_event_list(b,2)) = dbeta(birth_event_list(b,1),birth_event_list(b,2)) + 1;
end
rho = min(rho,1);
dt  = observed_cells{end-1}(1);
N   = size(D{1},1);
K   = length(D);
M   = 2000;           % M is the number of MCMC iterations
Q   = 10;             % Q is the number of trials

% ECOOM = squeeze(ECAET(:,1,:));
% ECOOM = ECOOM(:);
% ECOOM = ~isnan(ECOOM);
% BCOOM = beta(ECOOM);
% RCOOM = rho(ECOOM);
% RCOOM = sort(RCOOM);
% Rho   = cumsum(RCOOM)/sum(RCOOM);
% beta  = (RCOOM./Rho).*BCOOM;
% rho   = RCOOM;

all_l0 = ones(M,Q);
all_l1 = ones(M,Q);
all_pi = -10^5*ones(M,Q);

for q = 1:Q
    l0  = -30 + 2*(q-1)*(randn);
    l1  = -30 + 2*(q-1)*(randn); % LOGARITHM OF THE PARAMETERS
    for m = 2:M
        a_star  = l0 + (0.075*(1+5/sqrt(m))*randn);
        c_star  = l1 + (0.075*(1+5/sqrt(m))*randn);
        pi_star = likelihood_handler(rho,dbeta,D,dt,2,exp([a_star c_star]));
        pi_old  = likelihood_handler(rho,dbeta,D,dt,2,exp([l0 l1]));
        acc     = exp(pi_star-pi_old);  % acceptance rate
        tmp = rand;
        if tmp < acc                    % acceptance check
            l0 = a_star;
            l1 = c_star;
        end
        all_l0(m,q)  = l0;
        all_l1(m,q)  = l1;
        all_pi(m,q) = pi_old;
    end
end
%%
o0K = o0 + sum(death_event_list(:,2));
o1K = o1 + sum(death_event_list(:,1)-death_event_list(:,2));

death_pdf = @(r) (1/beta(o0K,o1K))*(1-exp(-dt*r)).^(o0K-1).*exp(-(o1K-1)*dt*r);
r = linspace(0,0.00001,10001);
plot(r,death_pdf(r))


%%
facit = observed_cells{end-1}(2:4);
[~,I] = max(mean(all_pi(M/2:end,:)));
max_likelihood = max(all_pi(:,I));
[~,Indices] = sort(all_pi(:,I));
Indices = Indices(end-M/10:end);
subplot(2,1,1)
histogram(all_l0(Indices,I))
hold on
plot([log(facit(1)) log(facit(1))],[0 40],'k--')
xlim([-20 0])
subplot(2,1,2)
histogram(all_l1(Indices,I))
hold on
plot([log(facit(2)) log(facit(2))],[0 40],'k--')
xlim([-20 0])


%%
ECOOM = squeeze(ECAET(:,1,:));
ECOOM = ECOOM(:);
ECOOM = ~isnan(ECOOM);
BCOOM = dbeta(ECOOM);
RCOOM = rho(ECOOM);
RCOOM = sort(RCOOM);
Rho   = cumsum(RCOOM)/sum(RCOOM);
BCOOM = (RCOOM./Rho).*BCOOM;

%%

unifrho = linspace(0,1,length(Rho))';

subplot(2,1,1)
plot(RCOOM,Rho,'r')
hold on
plot(unifrho,unifrho,'k')
xlim([0.01 0.99])
subplot(2,1,2)
plot(RCOOM,RCOOM./Rho,'b')
xlim([0.01 0.99])


%%
close all
plot(rho,dbeta)





