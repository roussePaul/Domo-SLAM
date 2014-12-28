function opt = setIfEmpty(opt, field, value)

if isfield(opt,field)~=1
    opt = setfield(opt, field, value);
end
