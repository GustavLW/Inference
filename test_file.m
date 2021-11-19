clc

theta = [0.25 3];
Vmin  = 0.25;           % depth of potential
alpha = 3;              % steepness of potential
xi    = [1 1]';
xj    = [0 0]';
dx    = 0.01;

Aik = explicit_hessian(xi,xj,dx,theta);
aik = explicit_gradient(xi,xj,dx,theta);

 testcell = cell(2,1);
 testcell{1} = cell(4,1);
 testcell{2} = [3 4 5];
 
 %%
 clc
 theta = [0.2 3];
 load dataset_20211119_1035.mat
 %

 dx = 10^(-3);
 alpha0 = 0;
 beta0 = 0;
 
 [alphaK,betaK] = diffusion_shell(observed_cells,theta,dx,alpha0,beta0);
 plot(sort(sqrt(betaK./(alphaK))))
 hold on
%
%%
ren_MSD = sort(sqrt(betaK./(alphaK)));
%%

histogram(min_metod,40)
