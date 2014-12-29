%[nbrOfFeatures,sizeOfState,indexFeatures,nbrOfParamPerFeature,sizeOfMeasures] = defSizes(mu)
function [nbrOfFeatures,sizeOfState,indexFeatures,nbrOfParamPerFeature,sizeOfMeasures] = defSizes(mu)

indexFeatures = 5; % index of last state definition
nbrOfParamPerFeature = 2; % number of parameters per features
sizeOfMeasures = 1;

sizeOfState = size(mu,1);
nbrOfFeatures = (sizeOfState-indexFeatures+1)/nbrOfParamPerFeature;