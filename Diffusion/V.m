function v = V(xi,xj,theta)
type = 1;

if type == 0
    Vmin = theta(1);
    alpha = theta(2);
    
    r = max(0.8,norm(xi-xj));

    
    De     = -Vmin/(2^(-2)-2^(-1));
    a   = 2^(-1/alpha);
    phi = 1/r;
    phi0 = 1;
    v = 2*a*De*((phi/phi0).^(2*a)-2*(phi/phi0).^a);

elseif type == 1
    v = zeros(2,1);
    r = max(0.001,norm(xi-xj));
    De   = theta(1);
    a    = theta(2);
    phi  = exp(-r);
    phi0 = exp(-1);
    %v = De*((phi/phi0).^(2*a)-2*(phi/phi0).^a);
     vamp = 2*a*De*(exp(-a*(r-1)) - exp(-2*a*(r-1)));
     v(1) = (xi(1)-xj(1))/r*vamp;
     v(2) = (xi(2)-xj(2))/r*vamp;
elseif type == 2
    
    
end