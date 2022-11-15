% this code takes a cell line, with 48 plates, and plots histogram over the
% 100 best samples for lambda0, lambda1

clc
clear all
folder_of_inferences = [cd '\3123'];
inference_folder = dir(folder_of_inferences);
inference_folder = inference_folder(3:end);
panel = 1;
for l = 1:length(inference_folder)
    load([folder_of_inferences '\' inference_folder(l).name]);
    l0_tracks = save_this{1};
    l1_tracks = save_this{2};
    pi_tracks = save_this{3};
    [~,I] = max(mean(pi_tracks(500:end,:)));
    max_likelihood = max(pi_tracks(:,I));
    [~,Indices] = sort(pi_tracks(:,I));
    Indices = Indices(end-100:end);
    subplot(8,6,panel)
    histogram(l0_tracks(Indices,I))
    title(inference_folder(l).name(1:end-4),'Interpreter','none')
    xlim([-20 0])
    panel = panel + 1;
end