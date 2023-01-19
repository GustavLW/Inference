function [X,Y,Z] = penalty_check(plot1,plot2,vary,varyn)

k1v = -8:4;      % true value:  2.3026
k2v = -9:3;      % true value:  1.3863
l1v = -10:2;     % true value: -0.5978
l2v = -10:2;     % true value:  0.1823
a1v = -16:-4;    % true value: -5.5215
a2v = -20:-8;    % true value: -9.7212

paramatrix = [k1v; l1v; a1v; k2v; l2v; a2v];

ulim = ones(1,6);
ulim(plot1) = 13;
ulim(plot2) = 13;
ulim(vary)  = varyn;
blim = ones(1,6);
blim(vary)  = varyn;
pen_surf = @(theta) penalize_agent(0,0,theta,1);
Z = zeros(13);
for i1 =  blim(1):ulim(1)
    for i2 = blim(2):ulim(2)
        for i3 = blim(3):ulim(3)
            for i4 = blim(4):ulim(4)
                for i5 = blim(5):ulim(5)
                    for i6 = blim(6):ulim(6)
                        th       = [k1v(i1) l1v(i3) a1v(i5) k2v(i2)  l2v(i4)  a2v(i6)];
                        Z(i1,i2) = pen_surf(th);
                    end
                end
            end
        end
    end
end

[X,Y] = meshgrid(paramatrix(plot1,:),paramatrix(plot2,:));
