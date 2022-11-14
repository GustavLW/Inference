function W = wasserstein_gamma(a,b,sigma0)

F1 = @(x) abs(gamcdf(x,a,b) - heaviside(x-sigma0));
W = integral(F1,0,inf);
