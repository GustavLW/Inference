function W = wasserstein_gamma(a,b,sig,min_track_length)
if size(a,1) > 1 %does a laplace approximation of the multivariate distr.
                 % and reports the W-distance for that gaussian to true
                 % value
I    = find(a<min_track_length);
b(I) = [];
a(I) = [];

tau = b./a;
h   = 0.5*a./(tau).^2;
W   = norm((sig)-(sqrt(tau)));% + sum(1./h);
%W   = sqrt(W);

elseif size(a,1) == 1 % we can use gamma directly under this assumption
W   = integral(@(x) abs(gamcdf(x,a,b/(a^2))-heaviside(x-(mean(sig))^2)),0,2*mean(sig)^2);

end
