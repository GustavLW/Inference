function agent = create_agent(type,T)
if type == 1
    theta = [15*rand-15 5*rand];
elseif type == 2
    k1 =    - 5 + 9*rand;
    k2 = k1 - 3 + 2*rand;
    a1 =    - 5 + 5*rand;
    a2 = a1 - 3 - 5*rand;
    l1 =    - 5 + 5*rand;
    l2 = l1 + 1 + 3*rand;
    theta = [k1 l1 a1 k2 l2 a2];
end

agent = struct('sigma',NaN,'U_param',theta.*ones(T,length(theta)),...
    'U_velo',randn(size(theta)).*ones(T,length(theta))/10,...
    'fitness',-Inf*ones(T,1),'penalty',zeros(T,1),...
    'U_best',theta.*ones(T,length(theta)),'Bfitness',-Inf*ones(T,1),'Bpenalty',zeros(T,1));