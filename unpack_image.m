function pop_tmp = unpack_image(observed_cells,k)

N = length(observed_cells) - 2;
pop_tmp = NaN*ones(2,N);

for i = 1:N
   pop_tmp(:,i) = observed_cells{i}.location(:,k); 
end

