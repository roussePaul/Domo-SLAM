 %function [r1,t1,phi] = lineExtraction(X,n)

function [Z,C,cl] = lineExtraction(X,S,d)

%% 
verbose = 0;

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

for i=1:N
    [m,I] = find( (class==ClassVal(i)) .* mask);
    data = X(I,:);
    d_eig = norm(eig(cov(data)));
    
    if d_eig>2*d_noise_eig
        [Z1,C1] = houghLineExtraction(data,10);
        Z = [Z;Z1];
        C = [C;C1];
        cl(I) = i;
    end
    
end

