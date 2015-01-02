function [muO,sigmaO] = addFeature(mu,sigma,Z,C)

Ne = size(mu,1);

muO = [mu;Z'];
sigma_ = 1e-4 * ones(size(muO,1));
sigma_((Ne+1):end,(Ne+1):end) = 1e-1;
sigma_(1:Ne,1:Ne) = sigma;
sigmaO = sigma_;

end