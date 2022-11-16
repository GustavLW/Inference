function W = wasserstein_gamma(a,b,sig)

N  = size(a,1);

W = zeros(size(a,2),1);
for c = 1:size(a,2)
    tau  = a(:,c)./b(:,c);
    h    = 0.5*a(:,c)./(tau).^2;
    W(c) = norm((1/sig)^2-tau)^2 + sum(h);
end
W = sqrt(W)*sig^2/N;

%F1 = @(x) abs(gamcdf(x,a,b) - heaviside(x-sigma0));
%W = integral(F1,0,inf)/sigma0;
