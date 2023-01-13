clc
clear all

desired_cell_line = 'C:\Users\Gustav\Documents\MATLAB\forskning\Data to Gustav\Data to Gustav (original images)\558\3123\Analysis\CellData_221111_111214_extracted';
data_folder = dir(desired_cell_line);
data_folder = data_folder(3:end);
mkdir(desired_cell_line(96:99))
%%
tic
for experiment = 3:3:length(data_folder)
    disp(['Begin inference on dataset ' regexp(data_folder(experiment).name,'[^_]+', 'match','once') ';'])
    too_large = 0;
    load([desired_cell_line '\' data_folder(experiment-2).name])
    load([desired_cell_line '\' data_folder(experiment-1).name])
    load([desired_cell_line '\' data_folder(experiment).name])
    disp(['Detecting ' num2str(size(ECAET,1)) ' cells over ' num2str(size(ECAET,3)) ' images...'])
    disp(['Loaded files are ' data_folder(experiment).name ', ' data_folder(experiment-1).name ', ' data_folder(experiment-2).name '.'])
    o0  = 0;
    o1  = 0;
    if numel(ECAET)>4*10^6
        str_experiment = regexp(data_folder(experiment).name,'[^_]+', 'match','once');
        disp(['Experiment ' str_experiment ' automatically dismissed for being too large.'])
        too_large = 1;
    else
    try 
    D   = distances(ECAET);
    catch ME
        if (strfind(ME.identifier,'MATLAB:array:SizeLimitExceeded'))
            str_experiment = regexp(data_folder(experiment).name,'[^_]+', 'match','once');
            disp(['Experiment ' str_experiment ' is too large! Array too large'])
            too_large = 1;
        end
        if (strfind(ME.identifier,'MATLAB:nomem'))
            str_experiment = regexp(data_folder(experiment).name,'[^_]+', 'match','once');
            disp(['Experiment ' str_experiment ' is too large! Out of memory'])
            too_large = 1;
        end
    end
    end
    if too_large == 0
    kernel_window = 3.068891070360678;
    % räkna ut rho
    [rho,~,~] = density(D,kernel_window,0);
    % räkna ut beta
    N   = size(D{1},1);
    K   = length(D);
    beta  = zeros(N,K); 
    for b = 1:size(birth_event_list,1)
        beta(birth_event_list(b,1),birth_event_list(b,2)) = beta(birth_event_list(b,1),birth_event_list(b,2)) + 1;
    end
    rho = min(rho,1);
    dt  = 1200; %HARDCODED
    N   = size(D{1},1);
    K   = length(D);
    M   = 1000;           % M is the number of MCMC iterations
    Q   = 10;             % Q is the number of trials
    
    all_l0 = ones(M,Q);
    all_l1 = ones(M,Q);
    all_pi = -10^5*ones(M,Q);
    
    for q = 1:Q
        l0  = -10 + 5*(q-1)*(randn);
        l1  = -10 + 5*(q-1)*(randn); % LOGARITHM OF THE PARAMETERS
        for m = 2:M
            a_star  = l0 + (0.075*(1+5/sqrt(m))*randn);
            c_star  = l1 + (0.075*(1+5/sqrt(m))*randn);
            pi_star = likelihood_handler(rho,beta,D,dt,2,exp([a_star c_star]));
            pi_old  = likelihood_handler(rho,beta,D,dt,2,exp([l0 l1]));
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
    
    % save logic
    save_this = cell(3,1);
    save_this{1} = all_l0;
    save_this{2} = all_l1;
    save_this{3} = all_pi;
    str_experiment = regexp(data_folder(experiment).name,'[^_]+', 'match','once');
    save_str = [desired_cell_line(96:99) '_' str_experiment '.mat'];
    save([cd '\' desired_cell_line(96:99) '\' save_str],'save_this')
    disp(['Inference completed on dataset ' str_experiment])
    else
        save_this = cell(3,1);
    save_this{1} = NaN*ones(M,1);
    save_this{2} = NaN*ones(M,1);
    save_this{3} = NaN*ones(M,1);
    str_experiment = regexp(data_folder(experiment).name,'[^_]+', 'match','once');
    save_str = [desired_cell_line(96:99) '_' str_experiment '.mat'];
    save([cd '\' desired_cell_line(96:99) '\' save_str],'save_this')
    disp(['Inference failed on dataset ' str_experiment])
    end
    toc
    
end


%%


