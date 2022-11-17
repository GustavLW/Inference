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
td = [1 2 5 7 10 15 20 30];

doff = 0;
for dataset = 1:length(df)/2
    d = dataset + doff;
    if strcmp(df(2*d-1).name(1:13),df(2*d).name(1:13)) == 1
        load([DataFolder '\' df(2*dataset-1).name])
        load([DataFolder '\' df(2*dataset).name])

        base_dt = observed_cells{end}(1);
        sig     = observed_cells{end}(end);

        [alpha,beta,betaMSD,betaIMP] = diffusion_inference(observed_cells,observed_neighbours,td);
        h = figure('units','normalized','outerposition',[0 0 1 1],'visible','off');
        % MF posterior
        alpha_grand = sum(alpha,1);
        beta_grand = sum(beta,1);
        betaMSD_grand = sum(betaMSD,1);
        betaIMP_grand = sum(betaIMP,1);

        r = linspace((sig)^(-2)*0.5,(sig)^(-2)*5,1001);
        for c = 1:length(td)
            a  = alpha_grand(c);
            b  = 1./beta_grand(c);
            bM = 1./betaMSD_grand(c);
            bI = 1./betaIMP_grand(c);

            p = flip(1./r);
            y = gampdf(p,a,1/(a^2*b));
            z = gampdf(p,a,1/(a^2*bM));
            w = gampdf(p,a,1/(a^2*bI));
            M = [max(y) max(z) max(w)];
            subplot(5,2,c)
            plot(p,y,'r')
            hold on
            plot(p,z,'b')
            plot(p,w,'g')
            plot([1 1]*(sig)^(2),[-0.1 1.1]*max(M),'k--')
            axis([p(1) p(end) [-0.1 1.1]*max(M)])
            xlabel('\sigma^{2}')
            ylabel('Probability')
            title(['MF assumption. Minutes between observations: ' num2str(base_dt/60*td(c))])
        end
        W1 = zeros(size(td,2),1);
        W2 = zeros(size(td,2),1);
        W3 = zeros(size(td,2),1);
        min_track_length = 1;
        for c = 1:length(td)
            W1(c) = wasserstein_gamma(alpha_grand(c),beta_grand(c),sig,min_track_length);
            W2(c) = wasserstein_gamma(alpha_grand(c),betaMSD_grand(c),sig,min_track_length);
            W3(c) = wasserstein_gamma(alpha_grand(c),betaIMP_grand(c),sig,min_track_length);
        end
        subplot(5,2,9:10)
        stem(base_dt/60*td,W1,'r')
        hold on
        stem(base_dt/60*td,W2,'b')
        stem(base_dt/60*td,W3,'g')
        grid on
        M = [max(W1) max(W2) max(W3)];
        axis([0 3*td(c)+10 [-0.1 1.1]*max(M)])
        legend('Our method','MSD','Implicit','Location','northwest')
        xlabel('Minutes between observations')
        title('Wasserstein distance')
        figname = ['figure' num2str(dataset)];
        saveas(h,figname,'png');
    else
        disp('Something is off')
        doff = doff - 0.5;
    end
    disp(['Dataset' num2str(dataset) ' is done!'])
end
%% Laplace-Wasserstein 
clf
W1 = zeros(size(td,2),1);
W2 = zeros(size(td,2),1);
W3 = zeros(size(td,2),1);
min_track_length = 1;
for c = 1:size(td,2)
    W1(c) = wasserstein_gamma(alpha(:,c),beta(:,c),sig,min_track_length);
    W2(c) = wasserstein_gamma(alpha(:,c),betaMSD(:,c),sig,min_track_length);
    W3(c) = wasserstein_gamma(alpha(:,c),betaIMP(:,c),sig,min_track_length);
end

hold on
stem(base_dt/60*td,W1,'r')
stem(base_dt/60*td,W2,'b')
%stem(base_dt/60*td,W3,'g')
grid on
axis([0 3*td(c)+10 [-0.1 1.1]*max(max(W1),max(W2))])
legend('Our method','MSD''Implicit','Location','northwest')
xlabel('Minutes between observations')
%%
sample_cell = observed_cells{1};
K           = length(sample_cell.location);
clear sample_cell
N           = length(observed_cells)-1;
clf
for c = 1:length(td)

    log_modes = log10(sqrt(beta./(alpha)));
    MSD_modes = log10(sqrt(betaMSD./(alpha)));
    IMP_modes = log10(sqrt(betaIMP./(alpha)));
    test = log_modes(:,c);
    test = test(~isnan(test));
    test1 = MSD_modes(:,c);
    test1 = test1(~isnan(test1));
    test2 = IMP_modes(:,c);
    test2 = test2(~isnan(test2));
    kern = @(x,r) normpdf(r,x,1/(length(log_modes))^(5/8));
    r = linspace(log10(sig)-1.5,log10(sig)+1.5,1001);
    y = zeros(size(r));
    y1 = y;
    y2 = y;
    for i = 1:length(test)
        y  = y  + kern(test(i),r);
        y1 = y1 + kern(test1(i),r);
        y2 = y2 + kern(test2(i),r);
    end
    subplot(4,2,c)
    hold on
    xlabel('log_{10}')
    ylabel('Distribution of modes')
    h0 = histogram(test,'FaceColor',[0.75 0 0],'Normalization','probability','FaceAlpha',0.35);
    h1 = histogram(test1,'FaceColor',[0 0 0.75],'Normalization','probability','FaceAlpha',0.35);
    h2 = histogram(test2,'FaceColor',[0 0.75 0],'Normalization','probability','FaceAlpha',0.35);
    plot(r,max(h0.Values)*y/max(y),'r','LineWidth',1.5)
    plot(r,max(h1.Values)*y1/max(y1),'b','LineWidth',1.5)
    plot(r,max(h2.Values)*y2/max(y2),'g','LineWidth',1.5)
    M = [max(h0.Values) max(h1.Values) max(h1.Values)];
    axis([r(1) r(end) -0.1*max(M) 1.1*max(M)])
    plot(log10([sig sig]),[-0.1 1.1]*max(M),'k--','LineWidth',1.5)
    %legend('Our method','MSD','Location','northwest')
    title([num2str(N) ' cells. Minutes between observations: ' num2str(base_dt/60*td(c))])

end

%%
% clf
% 
% 
% 
% c = 8;
% min_track_length = 1;
% a = alpha(:,c);
% b = beta(:,c);
% bM = betaMSD(:,c);
% 
% I    = find(a<min_track_length);
% b(I) = [];
% bM(I) = [];
% a(I) = [];
% 
% tau = a./b;
% h   = 0.5*a./(tau).^2;
% tauM = a./bM;
% hM   = 0.5*a./(tauM).^2;
% 
% for i = 1:length(tau)
%     subplot(1,2,1)
%     x = linspace(0.0*(1/sig)^2,4*(1/sig)^2,501);
%     y = normpdf(x,tau(i),sqrt(1/h(i)));
%     z = gampdf(x,a(i),1/b(i));
%     hold off
%     plot(x,y,'r')
%     hold on
%     plot(x,z,'k')
%     plot((1/sig)^2*[1 1],1.1*[0 max(max(y),max(z))],'k--')
%     axis([x(1) x(end) 1.1*[0 max(max(y),max(z))]])
%     legend('Normal estimation','Underlying Gamma-distribution','True \tau')
%     title('Our method')
%     subplot(1,2,2)
%     y = normpdf(x,tauM(i),sqrt(1/hM(i)));
%     z = gampdf(x,a(i),1/bM(i));
%     hold off
%     plot(x,y,'b')
%     hold on
%     plot(x,z,'k')
%     plot((1/sig)^2*[1 1],[0 0.004],'k--')
%     axis([x(1) x(end) 0 0.004])
%     legend('Normal estimation','Underlying Gamma-distribution','True \tau')
%     title('MSD')
%     drawnow;
% end
tau = 100;
a = 50000;
b = tau/a;
r = linspace(tau/2,3*tau/2,1001);
y = gampdf(r,a,b);
subplot(2,1,1)
plot(r,y)
subplot(2,1,2)
p = flip(1./r)
z = gampdf(p,a,1/(a^2*b));
plot(p,z);

