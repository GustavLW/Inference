function [output1,output2,output3] = density(D,r,g)
N = size(D{1},1);
K = length(D);
if g == 0
    output1 = zeros(N,K);
    for k = 1:K
        Dtmp = D{k};
        for i = 1:N
            tmp = Dtmp(i,:);
            tmp = tmp(~isnan(tmp));
            tmp = tmp(tmp>0);
            if ~isempty(tmp)
                output1(i,k) = sum(r*exp(-r*tmp));
            else
                output1(i,k) = 0;
            end
        end
    end
    output1 = output1;
    output2 = [];
    output3 = [];
elseif g == 1
    output1 = zeros(N,K);
    output2 = zeros(N,K);
    output3 = zeros(N,K);
    for k = 1:K
        for i = 1:N
            tmp = squeeze(D(i,:,k));
            output1(i,k) = sum(exp(-r*tmp));
            output2(i,k) = sum(tmp.*exp(-r*tmp));
            output3(i,k) = sum((tmp.^2).*exp(-r*tmp));
        end
    end
else
    disp('Du kodar som en apa!!!')
end
