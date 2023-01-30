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


load 24jan23.mat
T = size(repeated_trials{1,1},1);
Q = size(repeated_trials,1);
A = length(repeated_trials{1,4});
%% COMPARE WINNING POTENTIAL COMPARED TO UNDERLYING
loser_distance = zeros(Q+1,A)+10^(-16);
winner_distance = zeros(Q+1,1)+10^(-16);
h = figure('units','centimeters','position',[0 0 20*0.8 25*0.8]);
for q = 1:Q
    filename = ['particle_swarm_result_' num2str(q) '.png'];
    x       = linspace(0,5,1001);
    xs      = linspace(0,5,21);
    facit = observed_cells{end-1}(5:end-(length(observed_cells)-2));
    for t = T:T
        subplot(3,2,q)
        hold off
        B_param = exp(repeated_trials{q,1}(t,:));
        plot(xs,U_pot(xs,facit),'ro')
        hold on
        plot(xs+0.125,U_pot(xs+0.125,B_param),'bd')
        
        for a = 1:A
            A_param = exp(repeated_trials{q,4}{a}.U_best(t,:));
            plot(x,U_pot(x,A_param),'magenta')
        end
        
        plot(x,U_pot(x,facit),'r','LineWidth',1.5)
        plot(x,U_pot(x,B_param),'b','LineWidth',1.5)
        plot([0.5 4],[0 0],'k')
        grid on
        axis([0.5 4 min(U_pot(x,facit))*[2 -5]])
        xlabel('r')
        ylabel('U(r;\theta)')
        title(['12 best propsed potentials for r_\theta=' num2str(q/4) '.'])
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
subplot(3,2,5:6)
rt  = 0:0.25:1;
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
xticks([0:0.25:1])
xlabel('r_\theta')
ylabel('log(|u(r)-u(r;\theta|)')
grid on
