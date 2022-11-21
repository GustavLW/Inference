function agent = create_agent(type,T)
if type == 1
    theta = [(0.00005 + 0.0099*rand/2) 1+5*rand];
elseif type == 2
    k1 = 2      + rand*18;
    k2 = k1/20  + rand;
    a1 = (0.002 + 0.098*rand)*3;
    a2 = a1/15  + 0.005*rand;
    l1 = 0.1    + 0.9*rand;
    l2 = 1.5*l1 + 2*rand;
    theta = [k1 l1 a1 k2 a2 l2];
end

agent = struct('sigma',NaN,'U_param',theta.*ones(T,length(theta)),'U_velo',(1 + 0.2*abs(randn(size(theta)))).*theta/100,'fitness',NaN,'U_best',theta,'Bfitness',-Inf);