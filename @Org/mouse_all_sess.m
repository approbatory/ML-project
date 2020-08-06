function res_list = mouse_all_sess(o, varname, mouse_name)

if isempty(mouse_name)
    filt_sess = true(size(o.mouse));
elseif ischar(mouse_name) && strcmp(mouse_name, 'restrict')
    filt_sess = o.default_filt(true);
elseif iscell(mouse_name) && strcmp(mouse_name{2}, 'restrict')
    filt_sess = strcmp(o.mouse, mouse_name{1}) & o.default_filt(true);
else
    filt_sess = strcmp(o.mouse, mouse_name);
end

n_sess = sum(filt_sess);
assert(n_sess > 1, 'no sessions, or just 1');

%sess_indices = find(filt_sess);

%res_list = arrayfun(@(sess_idx) o.sess_by_bins(varname, sess_idx),...
%    sess_indices, 'UniformOutput', false);
V = o.fetch(varname);
V = V(filt_sess);
for sess_idx = 1:numel(V)
    if ~iscell(V{sess_idx})
        continue;
    end
    num_dims = cellfun(@numel,V{sess_idx});
    num_dims(num_dims==0) = [];
    n = min(num_dims);
    V{sess_idx} = cellfun(@(x)sel_(x,n), V{sess_idx},...
        'UniformOutput', false);
    V{sess_idx} = cell2mat(V{sess_idx});
end
res_list = V;


sizes = cellfun(@(x)size(x,1), res_list);
max_size = max(sizes);

res_list = cellfun(@(x)padarray(x,max_size-size(x,1),'post'), res_list,...
    'UniformOutput', false);
%assert(isequal(sizes{:}),...
%    'unequal sizes, consider trimming');
res_list = cat(3, res_list{:});

end



function r = sel_(x, n)

if isempty(x)
    r = nan(n,1);
else
    r = x(1:n);
    r = r(:);
end
end