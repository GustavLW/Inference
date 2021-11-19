function Aik = AAijk(i,k,nei,observed_cells,dx,theta)
xik = observed_cells{i}.location(:,k);
Aik = zeros(2);
N = size(nei,1);

for j = 1:N
    if nei(i,j) == 1
        xjk = observed_cells{j}.location(:,k);
        Aik = Aik + explicit_hessian(xik,xjk,dx,theta);
    end
end