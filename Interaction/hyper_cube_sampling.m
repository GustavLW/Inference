function s1 = hyper_cube_sampling(Q,para_ranges)
% now featuring no hypercube!
 % variable ranges go here
Y = cell(2,2,2,2,2,2);
for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
                    for i6 = 1:2
                        Y{i1,i2,i3,i4,i5,i6} = [para_ranges(1,i1) para_ranges(2,i2) para_ranges(3,i3) para_ranges(4,i4) para_ranges(5,i5) para_ranges(6,i6)];
                    end
                end
            end
        end
    end
end

samples = randperm(64);

s1 = zeros(Q,6);
for q = 1:Q
    s1(q,:) = Y{samples(q)};
end