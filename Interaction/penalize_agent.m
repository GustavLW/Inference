function [penalty,dev] = penalize_agent(RD,forecast_array,U_param,t)


if iscell(forecast_array)
dev = 0;

for s = 1:length(forecast_array)
    ftmp    = forecast_array{s};
    RDf     = calculate_radial_distribution(ftmp,1);
    for k = 1:ftmp{end}
        Rf      = RDf{k};
        R       = RD{k};
        sim_PCF = Rf(:,2);
        tru_PCF = R(:,2);
        dev     = dev + norm(sim_PCF-tru_PCF); 
    end
end
dev = dev/(length(forecast_array)); % dev as in deviance
else
    dev = forecast_array;   % technically not an array now, just a number. 
                            % BIG slop :)
end


if length(U_param) > 4
    g1 = [exp(U_param(4)) - exp(U_param(1));
         exp(U_param(6)) - exp(U_param(3));
         exp(U_param(2)) - exp(U_param(5))];
else
    g1 = -exp(U_param(2));
end
Re       = 4;
epsilon = 10^(-2);
u  = @(y) U_pot(y,exp(U_param));
v  = @(y) V(y,[0;0],exp(U_param));
v1 = @(a,y)a(y);
g2 = [                  1*(abs(u(Re))-epsilon*(abs(u(1)) - epsilon));...
                                                   epsilon.^2-u(0);...
                                                    v1(v([0;0]),1);...
                                                              u(1);...
                                                   -v1(v([1;0]),1);...
                                                    v1(v([1;0]),1);...
     -(v1(v([1+epsilon/2;0]),1)-v1(v([1-epsilon/2;0]),1))/epsilon];
penalty = 2^(sqrt(t))*(sum((max(g1,0)).^2)+sum((max(g2,0)).^2));
