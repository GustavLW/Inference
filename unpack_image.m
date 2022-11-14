function pop_tmp = unpack_image(observed_cells,k)

N = length(observed_cells) - 1;
pop_tmp = zeros(2,N);

for i = 1:N
   pop_tmp(:,i) = observed_cells{i}.location(:,k); 
end

