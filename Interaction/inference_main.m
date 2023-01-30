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


%% GRAFIK FÖR ANIMERA LÖSNINGAR
clc
close all
figure('units','centimeters','position',[0 0 25 25]);
filename = 'particle_swarm_result.gif';
c = linspace(1,10,A);
zoom = 3;
for t = 2:T
    current_params = zeros(A,4*pot_type - 2);
    current_best   = current_params;
    for a = 1:A
        current_params(a,:) = agents{a}.U_param(t,:);
        current_best(a,:)   = agents{a}.U_best(t,:);
    end
    if size(agents{a}.U_param(t,:)) < 4
        hold off
        scatter(current_params(:,1),current_params(:,2),[],c,'*')
        hold on
        scatter(current_best(:,1),current_best(:,2),[],c,'o')
        scatter(best_ever_location(t,1),best_ever_location(t,2),25,'red','o')
        scatter(log(facit(1)),log(facit(2)),25,'red','s')
        axis([-15 -0 0 5])
        grid on
    else
        subplot(2,2,1)
        hold off
        scatter(current_params(:,1),current_params(:,4),[],c,'*')
        hold on
        scatter(current_best(:,1),current_best(:,4),[],c,'o')
        plot(best_ever_location(2:t,1),best_ever_location(2:t,4),'r--o')
        scatter(log(facit(1)),log(facit(4)),100,'red','s')
        axis([log(facit(1))-zoom log(facit(1))+zoom log(facit(4))-zoom log(facit(4))+zoom])
        grid on
        xlabel('k_1')
        ylabel('k_2')
        subplot(2,2,2)
        hold off
        scatter(current_params(:,2),current_params(:,5),[],c,'*')
        hold on
        scatter(current_best(:,2),current_best(:,5),[],c,'o')
        plot(best_ever_location(2:t,2),best_ever_location(2:t,5),'r--o')
        scatter(log(facit(2)),log(facit(5)),100,'red','s')
        axis([log(facit(2))-zoom log(facit(2))+zoom log(facit(5))-zoom log(facit(5))+zoom])
        grid on
        xlabel('l_1')
        ylabel('l_2')
        subplot(2,2,3)
        hold off
        scatter(current_params(:,3),current_params(:,6),[],c,'*')
        hold on
        scatter(current_best(:,3),current_best(:,6),[],c,'o')
        plot(best_ever_location(2:t,3),best_ever_location(2:t,6),'r--o')
        scatter(log(facit(3)),log(facit(6)),100,'red','s')
        axis([log(facit(3))-zoom log(facit(3))+zoom log(facit(6))-zoom log(facit(6))+zoom])
        grid on
        xlabel('a_1')
        ylabel('a_2')
        subplot(2,2,4)
        hold off
        x       = linspace(0,5,1001);
        xs      = linspace(0,5,21);
        facit = observed_cells{end-1}(5:end-(length(observed_cells)-2));
        B_param = exp(repeated_trials{q,1}(t,:));
        plot(xs,U_pot(xs,facit),'ro')
        hold on
        plot(xs+0.125,U_pot(xs+0.125,B_param),'bd')
        plot(x,U_pot(x,facit),'r')
        plot(x,U_pot(x,B_param),'b')
        plot([0.5 4],[0 0],'k')
        grid on
        axis([0.5 4 min(U_pot(x,facit))*[2 -5]])
        xlabel('Distance')
        ylabel('Potential energy')
    end
    drawnow;
    frame = getframe(1);
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
    if t == 1
         imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');

    end
end

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








