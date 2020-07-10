function [summary, summary_err] = report(o, name, subname, varargin)
p = inputParser;
p.addRequired('name', @ischar);
p.addRequired('subname', @ischar);
p.addParameter('reps', 'mean');
p.addParameter('bins', 'median');
p.addParameter('sessions', 'mean');
p.addParameter('permouse', false);
p.parse(name, subname, varargin{:});
r = p.Results;

L = load([name '.mat']);

%start within each session

for s_i = 1:o.num_sess
    
    my_res = L.res{s_i};
    my_res = cellfun(@(s)index(s,subname), my_res, 'UniformOutput', false);
    dims = L.dim_names{s_i};
    
    %reps
    [my_res, my_res_err] = aggregate('reps', dims, my_res, []);
    
    
    %i = find(strcmp('reps', dims),1);
    %if ~isempty(i)
    %    [my_res, my_res_err] = aggregate(i, my_res);
    %end
end


end

function x = index(s, subname)
    if isemtpy(s)
        x = [];
        return;
    end
    x = s.(subname);
end

function [my_res_agg, my_res_agg_err] = aggregate(dim_name, dims, my_res, my_res_err)
i = find(strcmp(dim_name, dims),1);
if isempty(i)
    my_res_agg = my_res;
    my_res_agg_err = [];
    return;
end

[my_res_agg, my_res_agg_err] = cell_mean_sem(my_res, i);

end