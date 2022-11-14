function theta = interaction_shell(observed_cells,theta,dt,dx)
% dt here is the time step to the next 
for k = 1:K-1
    nei = neighbours_matrix(observed_cells,k);
    for i = 1:N
        xik  = observed_cells{i}.location(:,k);
        xin  = observed_cells{i}.location(:,k+1);
        aik  = aijk(i,k,nei,observed_cells,dx,theta);
        Aik  = AAijk(i,k,nei,observed_cells,dx,theta);
        mik  = xik + aik*dt;
        S1ik = sqrt(dt)*(eye(2)+dt/2*Aik);
        S2ik = dt*sqrt(dt/12)*Aik;
        Sik  = S1ik'*S1ik + S2ik'*S2ik;
     

    end
end