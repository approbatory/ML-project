function data = make_xlsx(data, name)

F = fields(data);
n_fields = numel(F);

num_elems = zeros(1, n_fields);
for i = 1:n_fields
    fname = F{i};
    num_elems(i) = numel(data.(fname));
end
max_elems = max(num_elems);

for i = 1:n_fields
    fname = F{i};
    x = data.(fname);
    x = x(:);
    c_ = iscell(x);
    if c_
        x = string(x);
    end
    y = nan(max_elems, 1);
    if c_
        y = string(y);
    end
    y(1:numel(x)) = x;
    data.(fname) = y;
end

data = struct2table(data);

save_directory = fileparts(which(mfilename));
xlsx_filename = sprintf('%s.xlsx', name);
writetable(data, fullfile(save_directory, xlsx_filename));