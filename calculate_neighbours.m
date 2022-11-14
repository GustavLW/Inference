function nei_tmp = calculate_neighbours(pop_tmp)

N      = length(pop_tmp);
nei_tmp = zeros(N);
% this code calculate which cells are our neighbours at time k
for i = 1:N
    for j = 1:N
        if j~=i
            if norm(pop_tmp(:,i)-pop_tmp(:,j)) < 3
                nei_tmp(i,j) = 1;
            end
        end
    end
end


