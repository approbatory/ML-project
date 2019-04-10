function v = normalize(v, param)
if nargin == 1
    error('supply second argument as "norm"');
end
if ~strcmp(param, 'norm')
    error('unsupported param for normalize: %s', param);
end

v = v ./ norm(v);
end