function opt = initOption(opt)
if nargin<1
    opt = struct();
end

opt = setIfEmpty(opt,'showEstimate',1);
opt = setIfEmpty(opt,'showTrue',1);
opt = setIfEmpty(opt,'showOdometry',0);
opt = setIfEmpty(opt,'verbose',2);