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

dataset = 1;
load([DataFolder '\' df(dataset).name])
load([DataFolder '\' df(dataset+1).name])

td      = [1 2 5 10 15 20 30 40];

base_dt = observed_cells{end}(1);
sig     = observed_cells{end}(end);

[alpha,beta,betaMSD] = diffusion_inference(observed_cells,observed_neighbours,td);

W1 = wasserstein_gamma(alpha,beta,sig);
W2 = wasserstein_gamma(alpha,betaMSD,sig);
%%
hold off
plot(base_dt/60*td,W1,'ro')
hold on
plot(base_dt/60*td,W2,'bo')



%%
sample_cell = observed_cells{1};
K           = length(sample_cell.location);
clear sample_cell
N           = length(observed_cells)-1;
for c = 1:length(td)
    clf
log_modes = log10(sqrt(beta./(alpha-1)));
MSD_modes = log10(sqrt(betaMSD./(alpha-1)));
test = log_modes(:,c);
test1 = MSD_modes(:,c);
kern = @(x,r) normpdf(r,x,1/(length(log_modes))^(3/4));
r = linspace(log10(sig)-2,log10(sig)+2,1001);
y = zeros(size(r));
y1 = y;
for i = 1:length(test)
   y  = y  + kern(test(i),r);
   y1 = y1 + kern(test1(i),r);
end
%subplot(2,2,c)
plot(r,y,'r')
xlabel('log_{10}')
hold on
ylabel('Smoothed histogram')
plot(r,y1,'b')
axis([r(1) r(end) -10 1.1*max(max(y),max(y1))])
%plot([-3 -3],[-10 810],'k--')
plot(log10([sig sig]),[-10 1.1*max(max(y),max(y1))],'k--')
legend('Our method','MSD','Location','northwest')
title([num2str(N) ' cells. Minutes between observations: ' num2str(base_dt/60*td(c))])
drawnow;
pause(1)
end


