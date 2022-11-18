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
MF = 0; % do we want to assume that all sigma are equal when evaluating performance?
for dataset = length(df)/2:length(df)/2
    d = dataset + doff;
    if strcmp(df(2*d-1).name(1:13),df(2*d).name(1:13)) == 1
        load([DataFolder '\' df(2*dataset-1).name])
        load([DataFolder '\' df(2*dataset).name])
        
        base_dt = observed_cells{end}(1);
        sig     = observed_cells{end}(end-length(observed_cells)+2:end);
        
        [alpha,beta,betaMSD,betaIMP] = diffusion_inference(observed_cells,observed_neighbours,td);
        h = figure('units','centimeters','position',[0 0 16.8 21]);%,'visible','off');
        if MF == 1% MF posterior
            alpha_grand = sum(alpha,1);
            beta_grand = sum(beta,1);
            betaMSD_grand = sum(betaMSD,1);
            betaIMP_grand = sum(betaIMP,1);
            r = linspace((min(sig))^(-2)*0.5,(max(sig))^(-2)*5,1001);
            for c = 1:length(td)
                a  = alpha_grand(c);
                b  = 1/beta_grand(c);
                bM = 1/betaMSD_grand(c);
                bI = 1/betaIMP_grand(c);
                
                p = flip(1./r);
                y = gampdf(p,a,1/(a^2*b));
                z = gampdf(p,a,1/(a^2*bM));
                w = gampdf(p,a,1/(a^2*bI));
                M = [max(y) max(z) max(w)];
                subplot(ceil(length(td)/2) + 1,2,c)
                plot(p,y,'r','LineWidth',1.5)
                hold on
                plot(p,z,'b','LineWidth',1.5)
                %plot(p,w,'g')
                plot([1 1]*(min(sig))^(2),[-0.1 1.1]*max(M),'k--','LineWidth',1.5)
                plot([1 1]*(max(sig))^(2),[-0.1 1.1]*max(M),'k--','LineWidth',1.5)
                plot(([p(1) p(end)]),[0 0],'k','LineWidth',1.0)
                axis([p(1) p(end) [-0.1 1.1]*max(M)])
                xlh = xlabel('$\sigma^{2}$','Interpreter','latex');
                xlh.Position(2) = xlh.Position(2) + 0.1;
                ylabel('$\pi(\sigma^2|X)$','Interpreter','latex')
                yticks([])
                title(['Minutes between observations: ' num2str(base_dt/60*td(c))])
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
            subplot(ceil(length(td)/2) + 1,2,ceil(length(td))+1:ceil(length(td)) + 2)
            xticklabels(string(base_dt/60*td));
            
            stem(W1,'r')
            hold on
            stem(W2,'b')
            %stem(W3,'g')
            grid on
            xticks(1:length(W1));
            xticklabels(string(base_dt/60*td));
            axis([0.5 length(W1)+0.5 [-0.1 1.1]*max(max(W1),max(W2))])
            %legend('Our method','MSD','Implicit','Location','northwest')
            xlabel('Minutes between observations')
            title('Wasserstein distance')
            suptitle(['Posterior distributions for ' num2str(dataset)])
            figname = ['MFfigure' num2str(dataset)];
            saveas(h,strcat('Results\',figname),'png');
            
        else % Laplace-Wasserstein - but with delta:s lol
            W1 = zeros(size(td,2),1);
            W2 = zeros(size(td,2),1);
            W3 = zeros(size(td,2),1);
            min_track_length = 2;
            for c = 1:size(td,2)
                W1(c) = wasserstein_gamma(alpha(:,c),beta(:,c),sig,min_track_length);
                W2(c) = wasserstein_gamma(alpha(:,c),betaMSD(:,c),sig,min_track_length);
                W3(c) = wasserstein_gamma(alpha(:,c),betaIMP(:,c),sig,min_track_length);
            end
            for c = 1:length(td)
                
                log_modes = log(sqrt(beta./(alpha)));
                MSD_modes = log(sqrt(betaMSD./(alpha)));
                IMP_modes = log(sqrt(betaIMP./(alpha)));
                test = log_modes(:,c);
                test = test(~isnan(test));
                test1 = MSD_modes(:,c);
                test1 = test1(~isnan(test1));
                test2 = IMP_modes(:,c);
                test2 = test2(~isnan(test2));
                kern = @(x,r) normpdf(r,x,1/(length(log_modes))^(5/8));
                span = 5;
                r = linspace(log((min(sig))/span),log((max(sig))*span),1001);
                y = zeros(size(r));
                y1 = y;
                y2 = y;
                for i = 1:length(test)
                    y  = y  + kern(test(i),r);
                    y1 = y1 + kern(test1(i),r);
                    y2 = y2 + kern(test2(i),r);
                end
                subplot(ceil(length(td)/2) + 1,2,c)

                 h0 = histogram(test,'FaceColor',[0.75 0 0],'Normalization','probability','FaceAlpha',0.35);
                 hold on
                 h1 = histogram(test1,'FaceColor',[0 0 0.75],'Normalization','probability','FaceAlpha',0.35);
                 %h2 = histogram(test2,'FaceColor',[0 0.75 0],'Normalization','probability','FaceAlpha',0.35);
                 plot(r,max(h0.Values)*y/max(y),'r','LineWidth',1.5)
                 plot(r,max(h1.Values)*y1/max(y1),'b','LineWidth',1.5)
                 %plot(r,max(h2.Values)*y2/max(y2),'g','LineWidth',1.5)
                 M = [max(h0.Values) max(h1.Values) max(h1.Values)];
                 axis([r(1) r(end) -0.1*max(M) 1.4*max(M)])
                 plot(log([min(sig) min(sig)]),[-0.1 1.2]*max(M),'k--','LineWidth',1.5)
                 plot(log([max(sig) max(sig)]),[-0.1 1.2]*max(M),'k--','LineWidth',1.5)
                 plot(([r(1) r(end)]),[0 0],'k','LineWidth',1.0)
                 %legend('Our method','MSD','Location','northwest')
                 title( ['\fontsize{10} Minutes between observations: ' num2str(base_dt/60*td(c))])
                 grid on
                 xlh = xlabel('\fontsize{12} \sigma');
                 xlh.Position(2) = xlh.Position(2) + 0.1;
                 ylabel('Distribution of modes')
                 yticks([])
            end
            
            
            subplot(ceil(length(td)/2) + 1,2,(length(td))+1:(length(td))+2)
            xticklabels(string(base_dt/60*td));
            
            stem(W1,'r')
            hold on
            stem(W2,'b')
            %stem(W3,'g')
            grid on
            xticks(1:length(W1));
            xticklabels(string(base_dt/60*td));
            axis([0.5 length(W1)+0.5 [-0.1 1.1]*max(max(W1),max(W2))])
            %legend('Our method','MSD','Implicit','Location','northwest')
            xlabel('Minutes between observations')
            title('Sum of mode deviations')
            suptitle(['Distribution of posteriors modes for Experiment ' num2str(dataset)])
            figname = ['MODEfigure' num2str(dataset)];
            saveas(h,strcat('Results\',figname),'png');
            
        end
    else
        disp('Something is off')
        doff = doff - 0.5;
    end
    disp(['Dataset ' num2str(dataset) ' is done!'])
end

