function F = interactions(pop_tmp,nei_tmp,U_param)
F = zeros(size(pop_tmp));
N = length(pop_tmp);
for i = 1:N
    xi  = pop_tmp(:,i);
    for j = 1:N
        if nei_tmp(i,j) == 1
            xj = pop_tmp(:,j);
            F(:,i) = F(:,i) + pair_potential(xi,xj,U_param);
        end
    end
end

