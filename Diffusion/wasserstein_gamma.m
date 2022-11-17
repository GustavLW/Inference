function W = wasserstein_gamma(a,b,sig,min_track_length)
if size(a,1) > 1 %does a laplace approximation of the multivariate distr.
                 % and reports the W-distance for that gaussian to true
                 % value
I    = find(a<min_track_length);
b(I) = [];
a(I) = [];

tau = a./b;
h   = 0.5*a./(tau).^2;
W   = norm((1/sig)^2-tau)^2 + sum(1./h);
W   = sqrt(W);

elseif size(a,1) == 1 % we can use gamma directly under this assumption
W   = integral(@(x) abs(gamcdf(x,a,b/(a^2))-heaviside(x-(sig)^2)),0,2*sig^2);

end
