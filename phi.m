function v = phi(xi,xj,theta)
% this is a morse potential but it'll have to do for now
% you have to be aware of what potential is currenty in phi, and what
% parameters to feed it
Vmin  = theta(1);
alpha = theta(2);
r     = sqrt((xi(1)-xj(1)).^2+(xi(2)-xj(2)).^2);
v     = Vmin*(exp(-2*alpha*(r-1))-exp(-alpha*(r-1))); 
end