function lean_extract_output(fname)
%mat_struct = load(fname);
%mat_struct.output.info = [];
%mat_struct.output.spatial_weights = sparsify_spatial_weights(mat_struct.output.spatial_weights);

bak = [fname '.bak_23381949'];
copyfile(fname, bak);

m = matfile(fname);
m.Properties.Writable = true;

try_fields = {'output', 'sortoutput'};
existing_fields = who('-file', fname);


for i = 1:numel(try_fields)
    f = try_fields{i};
    if ismember(f, existing_fields)
        fld = m.(f);
        if isfield(fld, 'info')
            fld.info = [];
        end
        fld.spatial_weights = sparsify_spatial_weights(fld.spatial_weights);
        m.(f) = fld;
    end
end

[directory, base, ext] = fileparts(fname);
fname_out = fullfile(directory, [base '_lean' ext]);
clear m

movefile(fname, fname_out);
movefile(bak, fname);

%save(fname_out, '-struct', 'mat_struct');
end
