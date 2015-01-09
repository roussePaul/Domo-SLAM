 %function [Z,C,cl] = lineExtraction(X,S,d)

function [Z,C,cl] = lineExtraction(X,S,d,alpha,s,verbose)

%% 
if nargin<6
    verbose=0;
end

k_min=10;

%% Nearest classification

[class,type] = dbscan2(X,k_min,d);

mask = (type>=0);

[m,I_type] = find(mask);

cl = class*0;
ClassVal = unique(class(I_type));
N = size(ClassVal,2);
Z=[];
C=[];
d_noise_eig = norm(eig(S));

if verbose == 1
    disp('Max dist clustering:');
    d
    disp('Number of class detected:');
    N
    figure
    gscatter(X(:,1),X(:,2),class);
end

j=1;
for i=1:N
    I = find( (class==ClassVal(i)) .* mask);
    data = X(I,:);
    d_eig = norm(eig(cov(data)));
    
    if d_eig>50*d_noise_eig
        [Z1,C1] = covarianceLineExtraction(data,1,alpha,s);
        Z(j,:) = Z1;
        C(:,:,j) = C1;
        j=j+1;
        cl(I) = i;
    end
    
end

