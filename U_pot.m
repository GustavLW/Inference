function v = U_pot(r,theta)
type = 1; % if we want LJ, set to 0. otherwise ignore. this is the PRIMITIVE of the interaction
% the POTENTIAL ENERGY
if length(theta) == 6
   type = 2; 
end
if type == 0
    Vmin   = theta(1) + 0*theta(2);
    alpha    = theta(2);
% LJ variables
D     = -Vmin/(2^(-2)-2^(-1));
rho   = 2^(-1/alpha);

v =  D*((rho./r).^(2*alpha)-(rho./r).^alpha); % L-J potential
elseif type == 1
    Vmin  = theta(1);
    alpha = theta(2);
    v     = Vmin*(exp(-2*alpha*(r-1))-2*exp(-alpha*(r-1))); % Morse potential
elseif type == 2
p1   = 1./(1+exp(-theta(1)*(r.^2-theta(2)^2)));
p2   = 1./(1+exp(-theta(4)*(r.^2-theta(5)^2)));
v = -theta(3)*(p1-1) + theta(6)*(p2-p1);
end