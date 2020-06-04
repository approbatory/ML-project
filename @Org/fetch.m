function res = fetch(o, varname)

if isfield(o.vars, varname)
    res = o.vars.(varname);
    fprintf('Fetched a saved variable: %s\n', varname);
elseif isfield(o.derived, varname)
    each = @(f,c)cellfun(f,c,'UniformOutput',false);
    
    varlist = o.derived.(varname).v;
    func = o.derived.(varname).f;
    my_vars = each(@(v)o.fetch(v), varlist);
    res = cell(1, o.total_sessions);
    for sess_idx = 1:o.total_sessions
        res{sess_idx} = cell(1,40);
        for b_idx = 1:40
            lengths = cellfun(@(x)numel(x{sess_idx}), my_vars);
            if any(lengths < b_idx & lengths ~= 1)
                continue;
            end
            v = each(@(x)lax_index(x,sess_idx,b_idx), my_vars);
            if any(cellfun(@isempty, v))
                continue;
            end
            res{sess_idx}{b_idx} = func(v{:});
        end
    end
    fprintf('Fetched a derived variable: %s\n', varname);
else
    fprintf('The variable %s does not exist: %s\n', varname);
    res = [];
    error('no such var');
end

end

function r = lax_index(x, s, b)
r = x{s};
if iscell(r)
    r = r{b};
end
end