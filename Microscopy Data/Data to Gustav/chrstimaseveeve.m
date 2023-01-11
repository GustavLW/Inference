clc 
clear all
tmp_dir = dir([cd '\Tracks\558']);
tmp_dir = tmp_dir(3:end);

dstring = tmp_dir(1).name;
data = readtable([cd '\Tracks\558\' dstring]);

clean = 1;

celldata = table2array(data);
N = max(data.particle);
K = max(data.frame);
ECAET = nan*ones(2,N,232);
for ik = 1:length(celldata)
    tmp_row = celldata(ik,:);
    if and(tmp_row(9) < 298, tmp_row(9) > 67)
        k            = tmp_row(9)-67;
        i            = tmp_row(10);
        ECAET(1,i,k) = tmp_row(2); 
        ECAET(2,i,k) = tmp_row(1);
    end
end

if clean == 1
    % remove shorties
    min_length = 2;
    for i = 1:size(ECAET,2)
        tmp_exist = squeeze(~isnan(ECAET(1,i,:)));
        tmp_index = find(tmp_exist==1);
        if size(tmp_index)>0
            obs_span = [tmp_index(1) tmp_index(end)];
            if (obs_span(2) - obs_span(1)) < min_length % is the track short?
                if obs_span(1) + min_length > size(ECAET,3) % and is it so for thisreason?
                elseif obs_span(2) - min_length < 1 % or this reason?
                else
                    ECAET(:,i,:) = nan;
                end
            end
        end
    end
    % remove shit
    i = 1;
    while i <= size(ECAET,2)
        if sum(isnan(ECAET(1,i,:))) == size(ECAET,3)
            ECAET(:,i,:) = [];
            %disp(['Row removed. Current cell count is ' num2str(size(ECAET,2))])
        else
            i = i + 1;
            %disp('No cell removed! We are moving forward lads')
        end
    end
    % fill gap
    clc
    for i = 1:size(ECAET,2)
        A = squeeze(ECAET(1,i,:));
        B = find(~isnan(A));
        if max(diff(B)) > 1
            k0 = B(1);
            k1 = B(end);
        else
            k0 = max(1,min(B));
            k1 = min(max(B),size(ECAET,3));
        end
        for k = k0:k1-1
            now  = squeeze(ECAET(1,i,k));
            next = squeeze(ECAET(1,i,k+1));

            if and(~isnan(now),isnan(next))
                ECAET(1,i,k+1) = ECAET(1,i,k);
                ECAET(2,i,k+1) = ECAET(2,i,k);
            end
        end
    end
end
%%
N = size(ECAET,2);
K = size(ECAET,3);
rdf_data = [];
for i = 1:N
    for k = 1:K
        if ~isnan(squeeze(ECAET(1,i,k)))
            rdf_data = [rdf_data; ECAET(1,i,k) ECAET(2,i,k) k];
        end
    end
end
[result,t] = gr2D(rdf_data,0,1,0,0,250,2);


%%
r = result(:,1);
g = smoothdata(result(:,2));
plot(r,g/mean(g(120:end)),'k')
xlabel('r (pixels)')
ylabel('g(r)')
grid on

[~,I] = max(g);
r0 = r(I);
ECAET = ECAET/r0;

save("extracted_tracks.mat",'ECAET')
%%
clc
N = size(ECAET,2);
K = size(ECAET,3);
h = figure('units','centimeters','position',[0 0 2*14.08 2*10.40]);
bmfolder = dir([cd '\Images\558\' dstring(1:2)  '\binary_mask']);
bmages   = bmfolder(71:length(bmfolder)-70);
imfolder = dir([cd '\Images\558\' dstring(1:2)]);
images   = imfolder(71:length(imfolder)-71);
for k = 1:K
bm_tmp = imread([cd '\Images\558\' dstring(1:2) '\binary_mask\' bmages(k).name]);
bm_tmp = im2double(bm_tmp);
im_tmp = imread([cd '\Images\558\' dstring(1:2) '\' images(k).name]);
im_tmp = im2double(im_tmp);
cm_tmp         = im_tmp;
em_tmp         = find(im_tmp==0);
cm_tmp(em_tmp) = max(max(cm_tmp));
fm_tmp         = (imgaussfilt(bm_tmp,5)/max(max(imgaussfilt(bm_tmp,5))) + imgaussfilt(bm_tmp,15)/max(max(imgaussfilt(bm_tmp,15))) + imgaussfilt(bm_tmp,35)/max(max(imgaussfilt(bm_tmp,35))))/3;
fm             = fm_tmp.*imgaussfilt(im_tmp./cm_tmp,4);
filtered_image = uint8(fm*255);
hold off
imshow(filtered_image)
hold on
scatter(squeeze(ECAET(1,:,k)),squeeze(ECAET(2,:,k)),'red','filled')
%     for i = 1:N
%         k0 = max(k-0,1);
%         plot(squeeze(ECAET(1,i,k0:k)),squeeze(ECAET(2,i,k0:k)),'y-o')
%     end
     drawnow;
end
function [result,t] = gr2D(data,dataType,dt,interactive,rmin,rmax,deltar)
% pair distribution function
%
% MODIFICATION HISTORY:
%     12-10-01   Eric R. Weeks
% see:  http://www.physics.emory.edu/~weeks/idl/gofr0.html
%
%     10-1-2018  Transmuted to Matlab by Rui Zhang, Dpe. of Math, Villanova Univ.
%                 Email: rzhang5@villanova.edu
%
% formatted for tabstop=4
%
%  usage:
%
% function ericgr2d,data,rmin=rmin,rmax=rmax,deltar=deltar,track=track
% assumes pretrack data unless /track is used
%
% IDL> tr=read_gdf('tracked.data.set')
% IDL> gr=ericgr2d(tr,/track,rmin=0.0,rmax=10.0)
% IDL> plot,gr(0,*),gr(1,*),xtitle='r',ytitle='g(r)'
%
% dataType: 1 for trackfile. The data output by track()
%                x   y   t  id
%           0 for orginal data. without using track(). 
%                x   y   t
%tic
if nargin < 4, interactive = 0; end % set 1 if you want to watch movie
if nargin < 5, rmin = 0; end
if nargin < 6, rmax = 10; end
if nargin < 7, deltar = 0.01; end
if and(length(data(1,:)) == 3,length(data(:,1))>0)
    %disp('WE HERE')
    %disp(data)
    nel = length(data(1,:)); % how many columns
    npts = length(data(:,1));% how many particles
    if dataType==0
        tel=nel;
    elseif dataType==1
        tel = nel - 1;
    else
        warning('only input 0 or 1');
        result=[];
        return
    end
    tmin = min(data(:,tel));
    tmax = max(data(:,tel));
    nr = (rmax - rmin)/deltar+1;
    rvec = (0:nr-1)'*deltar+rmin;
    rsqr = rvec.^2;
    result=zeros(nr,2);
    result(:,1) = rvec;
    rmin2=rmin^2;
    rmax2=rmax^2;
    x0 = min(data(:,1));
    x1 = max(data(:,1));
    y0 = min(data(:,2));
    y1 = max(data(:,2));
    density=npts/(tmax-tmin+1)/((x1-x0)*(y1-y0));
    %disp(strcat('number density = ',num2str(density)));
    for t=tmin:tmax
        temp = zeros(nr,1);
        w = find(data(:,tel) == t);
        nw=length(w);
        if nw > 0
            one = ones(nw,1);
            wtemp1 = find(data(w,1)>x0+rmax & data(w,1)<x1-rmax);
            wtemp2 = find(data(w,2)>y0+rmax & data(w,2)<y1-rmax);
            w4 = intersect(wtemp1,wtemp2);
            nw4 = length(w4);
            flag = one;
            if nw4>0
                flag(w4)=0;
                for i = 1:nw
                    pos0 = data(w(i),1:2);   %reference point
                    dd = (one*pos0)-data(w,1:2);
                    dis = sum(dd.^2,2);      % distance squared
                    dis(i) = 9e9;
                    w2 = find(dis>rmin2 & dis<rmax2);
                    nw2=length(w2);
                    if nw2>0
                        newdis = sqrt(dis(w2));
                        thehisto = histcounts(newdis,(rmin-deltar:deltar:rmax),'Normalization','probability');
                        thehisto = thehisto';
                        if (flag(i)<0.5)   %it's far from the edges
                            theta = 2*pi;
                        else
                            % now need to set correction factor based on location of
                            % pos0 -- ifd it's near corners. Check all four quadrants
                            tx=0;ty=0;tx2=0;ty2=0;
                            %checkquadrant,pos0,xref,yref,hix=hix,hiy=hiy,nr,rvec,rsqr,rmax,thetax=thetax,thetay=thetay
                            theta1 = checkquadrant(pos0,x0,y0,0,0,nr,rvec,rsqr,rmax,tx,ty);
                            theta2 = checkquadrant(pos0,x1,y0,1,0,nr,rvec,rsqr,rmax,tx2,ty);
                            theta3 = checkquadrant(pos0,x1,y1,1,1,nr,rvec,rsqr,rmax,tx2,ty2);
                            theta4 = checkquadrant(pos0,x0,y1,0,1,nr,rvec,rsqr,rmax,tx,ty2);
                            theta = theta1+theta2+theta3+theta4;
                        end
                        area = theta.*rvec*deltar; %area of each ring
                        w3 = find(area<1e-9);
                        nw3=length(w3);
                        if nw3>0
                            area(w3)=9e9;    % avoid divide-by-zero
                            temp = temp +thehisto./area;
                        end
                    end
                end
                temp = temp/nw/density;
            end
        end
        result(:,2)=result(:,2)+temp;
        if interactive == 1
            if mod(t,5)==0
                scatter(result(:,1),result(:,2)/(dt*(t-tmin+1)),'.')
                pause(0.001)
            end
        end
    end
    tmax=dt*tmax;tmin=dt*tmin;
    result(:,2) = result(:,2)/(density*(tmax-tmin+1));
    %plot(result(:,1),result(:,2))
    %toc
else
    result = [0 0];
end
end
%--------------------------------------------------------------------------
% local function
function theta = checkquadrant(pos0,xref,yref,hix,hiy,nr,rvec,rsqr,rmax,thetax,thetay)
% This function calculates the angle(in radians) of the arc that falls
% within this quadrant; used later for normalization.
%
% Clearly, this subroutine is the slowest part of the program. I think it's
% optimized but worth double-checking. Could possibly speed up the acos
% calls with a lookup table?
%
% pos0 is the test-point
% xref & yref are the edges that may be closest
% hix and hiy indicate the high-value edges x1 & y1 rather than the
%             low-value edges x0 & y0
% nr,rvec,rsqr,ramax all same as regular program
% thetax & thetay are used to pass variables which will be used for other
%                 quadrants and thus save time recalculating them
ind1 = rvec<=0.001;
rvec2=rvec;
rvec2(ind1,:)=0.001;% to avoid divide-by-zero
if thetax == 0
    if hix == 1
        if (xref - rmax) >= pos0(1)
            thetax = zeros(nr,1);
        else
            xprime = min(abs(zeros(nr,1)+xref-pos0(1)),rvec);
            thetax = acos(xprime./rvec2);
        end
    else  
        if (xref + rmax) < pos0(1)
            thetax = zeros(nr,1);
        else 
            xprime = min(abs(zeros(nr,1)+pos0(1)-xref),rvec);
            thetax = acos(xprime./rvec2);
        end
    end
end
if thetay == 0
    if hiy == 1
        if (yref-rmax) > pos0(1,2)
            thetay = zeros(nr,1);
        else 
            yprime = min(abs(zeros(nr,1)+yref-pos0(1,2)),rvec);
            thetay = acos(yprime./rvec2);
        end
    else
        if (yref+rmax)<pos0(1,2)
            thetay=zeros(nr,1);
        else
            yprime = min(abs(zeros(nr,1) + pos0(1,2) - yref),rvec);
            thetay = acos(yprime./rvec2);
        end
    end
end
theta = (zeros(nr,1) + 0.5*pi) - thetax - thetay;
dcorner = pos0 - [xref,yref];
cornerdist = sum(dcorner.^2);
w = find(rsqr>=cornerdist);
nw=length(w);
if nw > 0 
    theta(w) = 0;
end
end




