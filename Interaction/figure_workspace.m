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


load 28jan23.mat
T = size(repeated_trials{1,1},1);
Q = size(repeated_trials,1);
A = length(repeated_trials{1,4});
%% COMPARE WINNING POTENTIAL COMPARED TO UNDERLYING
clc
loser_distance = zeros(Q+1,A)+10^(-16);
winner_distance = zeros(Q+1,1)+10^(-16);
h = figure('units','centimeters','position',[0 0 20*0.8 25*0.8]);
for q = 1:Q
    filename = ['particle_swarm_result_' num2str(q) '.png'];
    x       = linspace(0,5,1001);
    xs      = linspace(0,5,21);
    facit = observed_cells{end-1}(5:end-(length(observed_cells)-2));
    for t = T:T
        subplot(4,1,q)
        hold off
        B_param = exp(repeated_trials{q,1}(t,:));
        plot(xs,U_pot(xs,facit),'ro')
        hold on
        plot(xs+0.125,U_pot(xs+0.125,B_param),'bd')
        
%         for a = 1:A
%             A_param = exp(repeated_trials{q,4}{a}.U_best(t,:));
%             plot(x,U_pot(x,A_param),'magenta')
%         end
%         
        plot(x,U_pot(x,facit),'r','LineWidth',1.5)
        plot(x,U_pot(x,B_param),'b','LineWidth',1.5)
        plot([0.5 4],[0 0],'k')
        grid on
        axis([0.5 3 min(U_pot(x,facit))*[2 -5]])
        xlabel('r')
        ylabel('U(r;\theta)')
        title(['12 best propsed potentials for \Delta t=' num2str(10*(2^(q-1))) ' minutes.'])
        legend('Underlying potential','Propsed by inference','Other agents')
        drawnow;
        winner_distance(q+1) = norm(U_pot(x(126:end),facit)-U_pot(x(126:end),B_param));
    end
    %saveas(h,filename)
    for a = 1:A
        A_param = exp(repeated_trials{q,4}{a}.U_best(t,:));
        loser_distance(q+1,a) = norm(U_pot(x(126:end),facit)-U_pot(x(126:end),A_param));
    end
end
subplot(4,1,Q+1)
rt  = 0:1/3:1;
rt2 = [rt, fliplr(rt)];
m1  = (log(min(loser_distance')));
m2  = mean(log(maxk(loser_distance',3)));
inBetween = [m1, fliplr(m2)];
hold off
fill(rt2, inBetween,[0.9 0.4 1]);
hold on
plot(rt,m1,'r','LineWidth',1.5)
plot(rt,m2,'r','LineWidth',1.5)
plot(rt,log(winner_distance),'k','LineWidth',2)
ylim([-15 0.75*log(max(max(loser_distance)))])
xticks([0:1/3:1])
xticklabels([0 10 20 40])
xlabel('\Delta t')
ylabel('log(|u(r)-u(r;\theta|)')
grid on
%% PERFORMANCE OF FACIT POTENTIAL FOR DIFFERENT TIME DILATIONS.
clc
Q = 12;
obs_array = cell(Q,2);
for q = 1:Q
obs_cop = observed_cells;
        for i = 1:length(observed_cells)-2
            obs_cop{i}.location = observed_cells{i}.location(:,q:q:end);
        end
        obs_cop{end-1}(1) = observed_cells{end-1}(1)*q;
        obs_cop{end}      = observed_cells{end}(q:q:end);
        obs_array{q,1}    = obs_cop;
end
%%
scsc = zeros(Q,1);
tic
for q = 1:Q
    agent = agents{1};
    agent.U_param = log(facit);
    S = 30;
    L = obs_array{q,1}{end-1}(1)/12;
    [agent,~] = score_agent(RD(q:q:end),agent,1,obs_array{q,1},L,S);
    scsc(q) = agent.fitness(1);
    disp(q)
    toc
end
%%
K_dt = floor(144./(1:12))';
plot(K_dt,scsc./K_dt,'k')
xlabel('Number of observations')
ylabel('Score of underlying potential')
axis([min(K_dt) max(K_dt) 0.5 4])
grid on
save('facit_score','scsc')
%% GRAFIK FÖR ANIMERA LÖSNINGAR
clc
close all
q = 3;
agents = repeated_trials{q,4};
best_ever_location = repeated_trials{q,1};

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