function aik = aijk(i,k,nei,observed_cells,dx,theta)

xik = observed_cells{i}.location(:,k);
aik = [0;0];
N = size(nei,1);

for j = 1:N
    if nei(i,j) == 1
        xjk = observed_cells{j}.location(:,k);
        aik = aik + explicit_gradient(xik,xjk,dx,theta);
    end
end