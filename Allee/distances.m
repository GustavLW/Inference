function D = distances(ECEAT)
N = size(ECEAT,1);
if N ==2
    N = size(ECEAT,2);
    K = size(ECEAT,3);
    D = cell(K,1);
    for k = 1:K
        tmp = NaN*ones(N);
        for i = 1:N
            for j = 1:N
                if j ~= i
                    if ~isnan(ECEAT(1,j,k))
                        if norm(ECEAT(:,i,k) - ECEAT(:,j,k)) < 9999
                            distance = norm(ECEAT(:,i,k) - ECEAT(:,j,k));
                            tmp(i,j) = distance;
                        end
                    end
                end
            end
        end
            D{k} = tmp;
    end

else
    K = size(ECEAT,3);
    D = cell(K,1);
    for k = 1:K
        tmp = NaN*ones(N);
        for i = 1:N
            for j = 1:N
                if j ~= i
                    if ~isnan(ECEAT(j,1,k))
                        if norm(ECEAT(i,:,k) - ECEAT(j,:,k)) < 9999
                            distance = norm(ECEAT(i,:,k) - ECEAT(j,:,k));
                            tmp(i,j) = distance;
                        end
                    end
                end
            end
        end
        D{k} = tmp;
    end
end
end
