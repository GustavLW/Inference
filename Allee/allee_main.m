clc
clear all
close all

TopFolder  = fileparts(fileparts(pwd));
SimFolder  = [TopFolder '\Simulation'];
DataFolder = [SimFolder '\Datasets'];
RealFolder = [fileparts(pwd) '\Microscopy Data'];
%dir(SimFolder)
addpath(([SimFolder, filesep]))
df = dir(DataFolder);
%df = dir(RealFolder);

df   = df(3:end);
o0   = 0;
o1   = 0;
save = 0;
MCMC = 0;
for d = length(df)-2:length(df)
    try
        load([DataFolder '\' df(d).name])
        [ECAET,birth_event_list,death_event_list] = get_locations_and_birth_events(observed_cells,1);
        facit = observed_cells{end-1}(2:4);
    catch ME
        load([RealFolder '\' df(d).name])
    end
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

    % MCMC part
    if MCMC == 1
        M   = 4000;           % M is the number of MCMC iterations
        Q   = 20;             % Q is the number of trials
        all_l0 = ones(M,Q);
        all_l1 = ones(M,Q);
        all_pi = -10^5*ones(M,Q);

        for q = 1:Q
            l0  = -15 + 10*rand;
            l1  = -15 + 10*rand; % LOGARITHM OF THE PARAMETERS
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
        if save == 1
            clc
            close all
            h = figure;

            %facit = [NaN NaN];
            [~,I] = max(mean(all_pi(M/2:end,:)));
            max_likelihood = max(all_pi(:,I));
            [~,Indices] = sort(all_pi(:,I));
            Indices = Indices(end-M/10:end);
            subplot(2,1,1)
            histogram(all_l0(Indices,I))
            hold on
            plot([log(facit(1)) log(facit(1))],[0 40],'k--')
            %xlim([-20 0])
            subplot(2,1,2)
            histogram(all_l1(Indices,I))
            hold on
            plot([log(facit(2)) log(facit(2))],[0 40],'k--')
            %xlim([-20 0])
            figname = ['FIG' num2str(d)];
            saveas(h,figname,'png');
        end
    elseif MCMC == 0
        [L0,L1] = meshgrid(-15:1:-5,-15:1:-5);
        Z = L0;
        for i = 1:size(L0,1)
            for j = 1:size(L1,1)
                Z(i,j)       = likelihood_handler(rho,dbeta,D,dt,2,exp([L0(i,j) L1(i,j)]));
            end
        end
        
       
    end
end

%%
clf
o0K = o0 + sum(death_event_list(:,2));
o1K = o1 + sum(death_event_list(:,1)-death_event_list(:,2));

death_pdf = @(r) (1/beta(o0K,o1K))*(1-exp(-dt*r)).^(o0K-1).*exp(-(o1K-1)*dt*r);
max_death = max(death_pdf(r));
r = linspace(0,3*facit(3),10001);
plot(r,death_pdf(r),'k','LineWidth',1)
hold on
plot(facit(3)*[1 1],1.5*[0 max_death],'r--','LineWidth',1)
xlim([0 max(r)])
ylim([0 1.25*max_death])
xlabel('\omega')
ylabel('p(\omega|X)')
grid on
%%
facit = observed_cells{end-1}(2:4);
for q = 1:Q
    plot3(exp(all_l0(2:end,q)),exp(all_l1(2:end,q)),-log(log(log(-all_pi(2:end,q)))))
    xlabel('\lambda_0')
    ylabel('\lambda_1')
    hold on
end
LAMBDA = likelihood_handler(rho,dbeta,D,dt,2,[facit(1) facit(2)]);
scatter3((facit(1)),(facit(2)),-log(log(log(-LAMBDA))),'red','filled')
xlim([0 2*facit(1)])
ylim([0 2*facit(2)])

%%

rl0 = round(log(facit(1)));
rl1 = round(log(facit(2)));
plot_range = 5;
[L0,L1] = meshgrid(linspace(rl0-plot_range,rl0+plot_range,201),linspace(rl1-plot_range,rl1+plot_range,201));
        Z = L0;
        for i = 1:size(L0,1)
            for j = 1:size(L1,1)
                Z(i,j)       = likelihood_handler(rho,dbeta,D,dt,2,exp([L0(i,j) L1(i,j)]));
            end
        end     
LAMBDA = likelihood_handler(rho,dbeta,D,dt,2,[facit(1) facit(2)]);
%%
for log_it = 0:1
    hold off
    surf(L0,L1,repeated_sqrt(Z,log_it),'EdgeColor','none')
    hold on
    scatter3(log(facit(1)),log(facit(2)),repeated_sqrt(LAMBDA,log_it),'red','filled')
    xlabel('\lambda_0')
    ylabel('\lambda_1')
    view(0,90)
    drawnow;
    pause(1)
end

function Z = repeated_sqrt(Z,iterations)
tmp = max(Z)-Z;
for i = 1:iterations
    tmp = log(1+(tmp).^(2^i));
end
Z = -tmp;
end