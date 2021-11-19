function nei = neighbours_matrix(observed_cells,k)
N   = length(observed_cells)-1;
nei = zeros(N);
for i = 1:N
    for j = 1:N
        xi = observed_cells{i}.location(:,k);
        xj = observed_cells{j}.location(:,k);
        if norm(xi-xj) < 3 && j ~= i
            if (observed_cells{i}.b_time <= k) && (observed_cells{i}.d_time > k)
                if (observed_cells{j}.b_time <= k) && (observed_cells{j}.d_time > k)
                    nei(i,j) = 1;
                end
            end
        end
    end
end

% dead cells are technically 


