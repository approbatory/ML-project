function [ds, X, ks, errf] = load_day(d, varargin)
p = inputParser;

p.addRequired('d', @isstruct);
p.addParameter('ds', {'deprobe', 'nocells'}, @iscell);
p.addParameter('data', {}, @iscell);
p.parse(d, varargin{:});
d = p.Results.d;
ds_opts = p.Results.ds;
data_opts = p.Results.data;

ds = quick_ds(fullfile(d.directory, d.day), ds_opts{:});
[X, ks, errf] = ds_dataset(ds, data_opts{:});
end