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
%%
df   = df(3:end);
o0   = 0;
o1   = 0;
save = 1;
MCMC = 0;
resolution = 401;
Z    = zeros(resolution);
priorize = 0;
for d = length(df)-3:length(df)
    clc
    close all
    h = figure('units','centimeters','position',[0 0 16.8 21]);
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
        if save == 1
            figname = ['MCMC_' num2str(d)];
            saveas(h,figname,'png');
        end
    elseif MCMC == 0
        rl0 = round(log(facit(1)));
        rl1 = round(log(facit(2)));
        plot_range = 3;
        l0linspace = linspace(rl0-plot_range,rl0+plot_range,resolution);
        l1linspace = linspace(rl1-plot_range,rl1+plot_range,resolution);
        [L0,L1] = meshgrid(l0linspace,l1linspace);
        %Z = L0;
        if priorize == 0
        Z   = zeros(resolution);
        o0  = 0;
        o1  = 0;
        end
        for i = 1:size(L0,1)
            for j = 1:size(L1,1)
                Z(i,j)       = Z(i,j) + likelihood_handler(rho,dbeta,D,dt,2,exp([L0(i,j) L1(i,j)]));
            end
        end
        [m1,i1] = max(Z);
        [m2,i2] = max(m1);

        best_loc_index = [i2 i1(i2)]; % x-koordinat y-koordinat
        ML = [l0linspace(best_loc_index(1)) l1linspace(best_loc_index(2))];
        LAMBDA = likelihood_handler(rho,dbeta,D,dt,2,[facit(1) facit(2)]);
        MAMBDA = likelihood_handler(rho,dbeta,D,dt,2,exp(ML));
        log_it = 1;
        subplot(3,2,1:4)
        hold off
        surf(L0,L1,Z,'EdgeColor','none')
        hold on
        scatter3(log(facit(1)),log(facit(2)),LAMBDA,'red','filled')
        scatter3(ML(1),ML(2),MAMBDA,'blue','filled')
        xlabel('log(\lambda_0)')
        ylabel('log(\lambda_1)')
        xlim([rl0-plot_range rl0+plot_range])
        ylim([rl1-plot_range rl1+plot_range])
        view(0,90)

        subplot(3,2,5:6)
        o0 = o0 + sum(death_event_list(:,2));
        o1 = o1 + sum(death_event_list(:,1)-death_event_list(:,2));
        death_pdf = @(r) (o0-1)*log((1-exp(-dt*r))) - (o1-1)*dt*r;
                r = linspace(facit(3)/9,3*facit(3),10001);
        max_death = min(max(death_pdf(r)),10^100);
        min_death = max(min(death_pdf(r)),-10^100);
        plot(r,death_pdf(r),'k','LineWidth',1)
        hold on
        plot(facit(3)*[1 1],[1.25*min_death 0.75*max_death],'r--','LineWidth',1)
        xlim([min(r) max(r)])
        ylim([1.05*min_death 0.95*max_death])
        xlabel('\omega')
        ylabel('p(\omega|X)')
        grid on

        if save == 1
            c = clock;
            c1 = num2str(c(1)); % år
            c2 = num2str(c(2)); % månad
            tmp = length(c2);
            if tmp == 1
                c2 = strcat('0',c2);
            end
            c3 = num2str(c(3)); % dag
            tmp = length(c3);
            if tmp == 1
                c3 = strcat('0',c3);
            end
            c4 = num2str(c(4)); % timme
            tmp = length(c4);
            if tmp == 1
                c4 = strcat('0',c4);
            end
            c5 = num2str(c(5)); % sekund
            tmp = length(c5);
            if tmp == 1
                c5 = strcat('0',c5);
            end
            figname  = strcat(c1,c2,c3,'_',c4,c5,'_likelihood');

            saveas(h,figname,'png');
        end
    end
end

%%

        death_pdf = @(r) (o0-1)*log((1-exp(-dt*r))) - (o1-1)*dt*r;
        r = linspace(facit(3)/9,3*facit(3),10001);
        max_death = min(max(death_pdf(r)),10^100);
        min_death = max(min(death_pdf(r)),-10^100);
        max(death_pdf(r))
        plot(r,death_pdf(r),'k','LineWidth',1)
        hold on
        plot(facit(3)*[1 1],[1.25*min_death 0.75*max_death],'r--','LineWidth',1)
        xlim([min(r) max(r)])
        ylim([1.05*min_death 0.95*max_death])
        xlabel('\omega')
        ylabel('p(\omega|X)')

function Z = repeated_sqrt(Z,iterations)
tmp = max(Z)-Z;
for i = 1:iterations
    tmp = log(1+(tmp).^(2^i))/4;
end
Z = -tmp;
end


