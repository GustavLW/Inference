function likel = evaluate_likelihood(rho,beta,D,dt,l0,l1)
% takes in one lambda0 and one lambda1 and return that pair's likelihood
if size(rho,2) > 1
    N = size(D{1},1);
    K = length(D);
    likel = 0;
    for k = 1:K-1
        for i = 1:N
            hik = dt*(1 - rho(i,k))*(l0+l1*rho(i,k));
            hik = max(eps,hik);
            likel = likel + beta(i,k)*log(hik) - hik;
        end
    end
else
    NK     = length(rho);
    likel  = 0;
    for ik = 1:NK
        hik = dt*(1 - rho(ik))*(l0+l1*rho(ik));
        hik = max(eps,hik);
        likel = likel + beta(ik)*log(hik) - hik;
    end
end