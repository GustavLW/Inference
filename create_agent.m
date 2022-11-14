function agent = create_agent(dt)
theta = [rand*20 rand*0.01/dt rand 2*rand rand*0.001/dt 5*rand];
agent = struct('sigma',NaN,'U_param',theta,'U_velo',theta/1000,'fitness',NaN,'U_best',theta,'Bfitness',-Inf);