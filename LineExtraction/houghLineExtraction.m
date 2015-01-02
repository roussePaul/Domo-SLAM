 %function [r1,t1,phi] = lineExtraction(X,n)
 % X vector of points, n minimum number of points for a line

function [Z,C] = houghLineExtraction(X,n)

minx = min(X(:,1));
maxx = max(X(:,1));
miny = min(X(:,2));
maxy = max(X(:,2));

N = 100;

f = max([(maxx-minx),(maxy-miny)]);
f = (N-1)/f;

idx = 1+floor((X(:,1)-minx)*f);
idy = 1+floor((X(:,2)-miny)*f);
BW = zeros(100);
for i=1:size(X,1)
    BW(idx(i),idy(i)) = 1;
end
[H, theta, rho] = hough(BW);
peaks = houghpeaks(H, 1);


I = -5:5;
np = size(peaks,1);
tr = zeros(np,1);
Z=zeros(np,2);
C=zeros(2,2,np);
mask = zeros(np,1);
for i=1:np
    
    r=peaks(i,1);
    t=peaks(i,2);   
    
    
    
    ctot=0;
    for k=I
        for l=I
            g = [k,l];
            if 0<r+k && r+k<=size(H,1) && 0<t+l && t+l<size(H,2)
                C(:,:,i) = C(:,:,i) + H(r+k,t+l)*(g')*g;
                ctot = ctot+H(r+k,t+l);
            end
        end
    end
    
    mask(i)=(ctot>=n);
    
    
    C(:,:,i) = C(:,:,i)/ctot;
    a = theta(t)/180*pi;
    a = (pi/2-a);
    d = rho(r)/f+cos(a)*minx+sin(a)*miny;
    if d<0;
        d = abs(d);
        a = a+pi;
    end
    Z(i,:) = [d,a];
    tr(i) = norm(eig(C(:,:,i)));
end

% 
% figure
% imagesc(H)
% 
% 
% figure
% plot(idx,idy,'.');
% hold on
% x = linspace(0,100,10);
% disp('Rt')
% r = rho(peaks(1,1))
% t = theta(peaks(1,2))/180*pi
% plot(x,(r-x*cos(t))/sin(t))
% 
% f
% 
% figure
% 
% plot(X(:,1),X(:,2),'.');
% hold on
% x = linspace(0,1,10);
% disp('Rt')
% t = theta(peaks(1,2))/180*pi
% r = rho(peaks(1,1))/f+cos(t)*minx+sin(t)*miny
% plot(x,(r-x*cos(t))/sin(t))
% 
% d
% a
% 
[m,I] = find(mask==1);

Z = Z(I,:);
C = C(:,:,I);
tr = tr(I);
[tr,I] = sort(tr,'ascend');

Z = Z(I,:);
C = C(:,:,I);
