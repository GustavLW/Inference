clc
clear all
clf
close all
load dataset_20221115_1720

theta = [0.000075 4];
%theta = [0 4];
sig   = 0.01;


sample_cell = observed_cells{1};
K           = length(sample_cell.location);
clear sample_cell
N           = length(observed_cells)-1;
ECEAT       = NaN*ones(N,2,K);
for i = 1:N
    for k = 1:K
        if k >= observed_cells{i}.b_time
            ECEAT(i,:,k) = observed_cells{i}.location(:,k);
        end
    end
end
    alpha = zeros(N,5);
    beta  = zeros(N,5);
    betaMSD = beta;
    td = 3*[1 4 9 16 25];
for c = 1:5
    skip = td(c);
    dt = observed_cells{end}*skip;
    for k   = skip:skip:K-skip
        x   = squeeze(ECEAT(:,:,k))';
        axn = squeeze(ECEAT(:,:,k+skip))';
        for i = 1:N
            xn        = squeeze(ECEAT(i,:,k+skip))';
            xi        = x(:,i);
            [mi,Si]   = get_coeff(xi,x,theta,dt);
            mi        = xi + mi;
            alpha(i,c)  = alpha(i,c) + 1;
            beta(i,c)   = beta(i,c) + 0.5*dot(xn-mi,Si\(xn-mi));
            betaMSD(i,c)= betaMSD(i,c) + 0.5*dot(xn-xi,xn-xi)/dt;
        end
    end
end

%%
W1 = zeros(N,5);
W2 = zeros(N,5);
for c = 1:5
    for i=1:N
        W1(i,c) = wasserstein_gamma(alpha(i,c),1/beta(i,c),1/sig^2);
        W2(i,c) = wasserstein_gamma(alpha(i,c),1/betaMSD(i,c),1/sig^2);
    end
end

clf
clc
plot(3*td,mean(log(W1)),'r')
hold on
plot(3*td,mean(log(W2)),'b')
legend('Our method','MSD','Location','northwest')
% for c = 1:5
%     subplot(5,1,c)
%     boxplot([log(W1(:,c)),log(W2(:,c))],'Notch','on','Labels',{'Our method','MSD'})
%     
% end
%legend('Our method','MSD','Location','northwest')
%%
clf
W1 = log(W1);
W2 = log(W2);
for c = 2:5 
    kern = @(x,r) normpdf(r,x,1/N^(17/32));
    r  = linspace(min(min(W1(:,c)),min(W2(:,c)))-2,max(max(W1(:,c)),max(W2(:,c)))+2,1001);
    y  = zeros(size(r));
    y1 = y;
    for i = 1:size(W1,1)
        y  = y  + kern(W1(i,c),r);
        y1 = y1 + kern(W2(i,c),r);
    end
    subplot(2,2,c-1)
    h = histogram(W1(:,c),'FaceColor','r','NumBins',20,'Normalization','count','FaceAlpha',0.2);
    xlabel('Logarithm of Wasserstein distances')
    ylabel('Smoothed histogram')
    hold on
    h1 = histogram(W2(:,c),'FaceColor','b','NumBins',20,'Normalization','count','FaceAlpha',0.2);
    n  = max(h.Values);
    n1 = max(h1.Values);
    y  = y*n/max(y); 
    y1 = y1*n1/max(y1);
    plot(r,y,'r','LineWidth',2)
    plot(r,y1,'b','LineWidth',2)
    grid on
    axis([r(1) r(end) -10/sum(y1) 1.1*max(max(y),max(y1))])
    legend('Our method','MSD','Location','northwest')
    title([num2str(3*td(c)) ' minutes between observations'])
end

%%
for c = 1:5
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
title([num2str(N) ' cells. Minutes between observations: ' num2str(3*3^(c-1))])
drawnow;
pause(1)
end

%%

