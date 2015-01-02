function [Z,C] = covarianceLineExtraction(X,n,alpha,s)

r = @(t) mean(cos(t)*X(:,1)+sin(t)*X(:,2));
c = @(t) cov(cos(t)*X(:,1)+sin(t)*X(:,2));

t1 = alpha;
t2 = alpha+pi/2;
t=t1;
if c(t1)>c(t2)
    t=t2;
end

Z = [r(t),t];
C = zeros(2);


