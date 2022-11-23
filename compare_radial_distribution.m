function dev = compare_radial_distribution(observed_cells,U_param)

RD = calculate_radial_distribution(observed_cells);
%plot(RD(:,1),RD(:,2)/mean(RD(401:end,2)),'b')
forecast_cells = observed_cells;
[alpha,beta,~,~] = diffusion_inference(observed_cells,1,U_param);
N                = length(observed_cells)-2;
K                = length(observed_cells{end});
for k = 1:K-1
    pop_tmp      = unpack_image(observed_cells,k);
    nei_tmp      = observed_cells{end}{k};
    sigma        = sqrt(beta./alpha)';
    pop_forecast = sample_from_next_state(pop_tmp,nei_tmp,U_param,sigma,1,observed_cells{end-1}(1)/6);
    for i = 1:N 
        forecast_cells{i}.location(:,k+1) = pop_forecast(:,i);
    end
end
forecast_cells{end-1}(2:end-3-(length(observed_cells)-2)) = U_param;
forecast_cells{end-1}(end-N+1:end)                        = sigma;
RDf = calculate_radial_distribution(forecast_cells);
%hold on
%plot(RDf(:,1),RDf(:,2)/mean(RDf(401:end,2)),'r')

sim_PCF = RDf(:,2)/mean(RDf(401:end,2));
tru_PCF = RD(:,2)/mean(RD(401:end,2));
dev     = norm(sim_PCF-tru_PCF);