function [mu,sigma,tablOutliers] = addFeatures(mu,sigma,Q,z,angleMeasure,tablOutliers,outlier,verbose)

if nargin<8
    verbose=0;      
end

np=5;
bufferSize = 100;


indOutliers = find(outlier~=0);
pts = getPosFromScan(mu(1:3)',z(indOutliers),angleMeasure(indOutliers));

tablOutliers = [tablOutliers;pts];
sTO = size(tablOutliers,1);


if verbose
    disp(['Size of the outlier table: ',num2str(sTO)])
end

if sTO>np
    
    [Z,C,I] = lineExtraction(tablOutliers,0.0002*ones(2),0.2,mu(4),sigma(4,4)+0.01,verbose);

    if size(Z,1)>0
        for i=1:size(Z,1)
            [mu,sigma] = addFeature(mu,sigma,Z(i,:),C(:,:,i));
            disp(['Feature added, r=',num2str(Z(i,1)),' theta=',num2str(Z(i,2)),' alpha=',num2str(mu(4))]);
        end
        ind = find(I==0);
        
        tablOutliers = tablOutliers(ind,:);
        
    end
    
    n = size(tablOutliers,1); 
    if n>bufferSize
        tablOutliers = tablOutliers((n-bufferSize):end,:);
    end

end
end