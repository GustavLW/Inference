function y = likelihood_handler(rho,beta,D,dt,model,x)
% rho, beta, D and dt are observations
% we insert "observed cells" in order to get some timings right and such
% model, 1 through 4, is the parameters considered
% they're stored in "x"
% evaluate likelihood currently does not support model 2 and 4
% it is evident that unconstrained optimization fucks up, so i'll attempt
% to implement a SQP adaptation.

if model == 1
    % corresponds to only a base division and logistic growth
    y = evaluate_likelihood(rho,beta,D,dt,x(1),0);
elseif model == 2
    % corresponds to allee effect
    y = evaluate_likelihood(rho,beta,D,dt,x(1),x(2));   
end
end